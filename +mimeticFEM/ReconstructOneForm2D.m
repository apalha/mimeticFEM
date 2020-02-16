function [reconstructedX, reconstructedY, xiRefinedGrid, etaRefinedGrid] = ReconstructOneForm2D(discreteOneForm,gridType,nReconstruction,mesh)
%ReconstructOneForm2D computes the reconstruction of a discretized 1-form
%   in 2 dimensions. It performs the reconstruction for all the elements
%   given.
%
%   Returns two matrices containing the values of the reconstructed x and y
%   components of the 1-form at the reconstructed nodes, following the same
%   structure as for the input of discreteOneForm.
%
%   USAGE
%   -----
%       [reconstructedX, reconstructedY, xiRefinedGrid, etaRefinedGrid] = 
%                   ReconstructOneForm2D(discreteOneForm,gridType,nReconstruction,
%                                        mesh)
%
%   INPUTS
%   ------
%       discreteOneForm   :: The 1-form discretized in each element.
%                            discreteOneForm(i,j) is the coefficient
%                            associated to the ith basis function of the jth
%                            element. Therefore columns contain the
%                            coefficients in each element. Since we are in
%                            2D, for spectral elements of order p, we have
%                            p(p+1) basis for each component, therefore there
%                            are 2p(p+1) basis functions (coefficients).
%                            First the xi component and then the eat
%                            component.
%                            (type: float64, size: [(p+1)(p+1),n])
%       gridType :: The type of grid. Can be:
%                'Lobatto', 'EGauss' or 'Gauss'
%                (type: string, size: single value)
%       nReconstruction :: The number of points to use in the x and y
%                          direction for the plotting. As an alternative
%                          a list of points xi\in[-1,1] can be given.
%                          (type: int64, shape: single value, [1,N] or [N,1])
%       mesh :: (type: mimeticFEM2.Mesh, size: single value)
%               A mesh object over which the 2-form is defined.
%
%
%   OUTPUTS
%   -------  
%       reconstructedX :: The reconstructed x-component of the 1-form at the reconstruction
%                        points given locally by xiRefinedGrid and
%                        etaRefinedGrid. n is the number of elements.
%                        (type: float64, size: [nReconstruction*nReconstruction, n])
%       reconstructedY :: The reconstructed y-component of the 1-form at the reconstruction
%                        points given locally by xiRefinedGrid and
%                        etaRefinedGrid. n is the number of elements.
%                        (type: float64, size: [nReconstruction*nReconstruction, n])
%       xiRefinedGrid :: The xi coordinates of the points of the reconstruction
%                        grid in the computational domain. 
%                        xiRefinedGrid \in [-1,1].
%                        (type: float64, size: [nReconstruction*nReconstruction,1])
%       etaRefinedGrid :: The eta coordinates of the points of the reconstruction
%                        grid in the computational domain. 
%                        etaRefinedGrid \in [-1,1].
%                        (type: float64, size: [nReconstruction*nReconstruction,1])
%   
%   Copyright 2009 Artur Palha

%   Revisions:  2009-11-04 (apalha) First implementation.
%               2009-11-05 (apalha) Created the function for reconstruction 
%                                   of one element.
%               2009-11-05 (apalha) Added the possibility to reconstruct at 
%                                   a given set of parametric points.
%               2009-11-06 (apalha) Added the possibility to reconstruct for
%                                   more than one element.
    
    % the number of elements
    nElements = size(discreteOneForm,2);
    
    % compute the number of degrees of freedom for xi and eta components
    dofOneFormXi = 0.5*size(discreteOneForm,1);
    dofOneFormEta = 0.5*size(discreteOneForm,1);
    
    % compute the order of the discretization
    p = 0.5*(-1 + sqrt(1+2*size(discreteOneForm,1)));
    
    if length(nReconstruction)>1
        xiRefined = nReconstruction;
        etaRefined = xiRefined;
        nReconstruction = length(nReconstruction);
    else
        % compute the local nodes where to compute the refinement
        xiRefined = (2/(nReconstruction-1))*(0:(nReconstruction-1))-1;
        etaRefined = xiRefined;
    end
    
    [xiRefinedGrid, etaRefinedGrid] = meshgrid(xiRefined,etaRefined);
    xiRefinedGrid = xiRefinedGrid(:);
    etaRefinedGrid = etaRefinedGrid(:);
    
    % compute the basis functions at the refined nodes
    
    % xi basis
    xiBasisRefinedXi = mimeticFEM.EdgePoly(xiRefined,p,gridType);
    etaBasisRefinedXi = mimeticFEM.LobattoPoly(etaRefined,p);
    
    xietaBasisXiKron = kron(xiBasisRefinedXi,etaBasisRefinedXi)';
    
    % eta basis
    xiBasisRefinedEta = etaBasisRefinedXi;
    etaBasisRefinedEta = xiBasisRefinedXi;
    
    xietaBasisEtaKron = kron(xiBasisRefinedEta,etaBasisRefinedEta)';
    
    % allocate memory space for reconstructed result
    reconstructedX = zeros([nReconstruction*nReconstruction nElements]);
    reconstructedY = zeros([nReconstruction*nReconstruction nElements]);
    
    if nElements>1
        for element=1:nElements
            % compute the metric at the refined points
            gEvaluated = mesh.g(element,xiRefinedGrid(:),etaRefinedGrid(:));
            
            % compute dPhiXdXi
            dXdXiEvaluatedMatrix = spdiags(mesh.dXdXi(element,xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);
            
            % compute dPhiXdEta
            dXdEtaEvaluatedMatrix = spdiags(mesh.dXdEta(element,xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);
            
            % compute dPhiYdXi
            dYdXiEvaluatedMatrix = spdiags(mesh.dYdXi(element,xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);
            
            % compute dPhiYdEta
            dYdEtaEvaluatedMatrix = spdiags(mesh.dYdEta(element,xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);
            
            % the combined basis in 2d already scaled with the metric
            xietaBasisXiX = dYdEtaEvaluatedMatrix*xietaBasisXiKron;
            xietaBasisXiY = -dXdEtaEvaluatedMatrix*xietaBasisXiKron;
            xietaBasisEtaX = -dYdXiEvaluatedMatrix*xietaBasisEtaKron;
            xietaBasisEtaY = dXdXiEvaluatedMatrix*xietaBasisEtaKron;

            % reconstruct
            reconstructedX(:,element) = xietaBasisXiX*discreteOneForm(1:dofOneFormXi,element) +...
                                        xietaBasisEtaX*discreteOneForm((dofOneFormXi+1):end,element);
            reconstructedY(:,element) = xietaBasisXiY*discreteOneForm(1:dofOneFormXi,element) +...
                                        xietaBasisEtaY*discreteOneForm((dofOneFormXi+1):end,element);
            
        end
    else
        % compute the metric at the refined points
        gEvaluated = mesh.g(1,xiRefinedGrid(:),etaRefinedGrid(:));

        % compute dPhiXdXi
        dXdXiEvaluatedMatrix = spdiags(mesh.dXdXi(1,xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);

        % compute dPhiXdEta
        dXdEtaEvaluatedMatrix = spdiags(mesh.dXdEta(1,xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);

        % compute dPhiYdXi
        dYdXiEvaluatedMatrix = spdiags(mesh.dYdXi(1,xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);

        % compute dPhiYdEta
        dYdEtaEvaluatedMatrix = spdiags(mesh.dYdEta(1,xiRefinedGrid(:),etaRefinedGrid(:))./gEvaluated,0,nReconstruction*nReconstruction,nReconstruction*nReconstruction);

        % the combined basis in 2d already scaled with the metric
        xietaBasisXiX = dYdEtaEvaluatedMatrix*xietaBasisXiKron;
        xietaBasisXiY = -dXdEtaEvaluatedMatrix*xietaBasisXiKron;
        xietaBasisEtaX = -dYdXiEvaluatedMatrix*xietaBasisEtaKron;
        xietaBasisEtaY = dXdXiEvaluatedMatrix*xietaBasisEtaKron;

        % reconstruct
        reconstructedX(:) = xietaBasisXiX*discreteOneForm(1:dofOneFormXi) +...
                                    xietaBasisEtaX*discreteOneForm((dofOneFormXi+1):end);
        reconstructedY(:) = xietaBasisXiY*discreteOneForm(1:dofOneFormXi) +...
                                    xietaBasisEtaY*discreteOneForm((dofOneFormXi+1):end);
    end
        
end
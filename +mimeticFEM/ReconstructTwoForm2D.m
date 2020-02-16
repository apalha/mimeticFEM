function [reconstructed,xiRefinedGrid,etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm,gridType,mesh,nReconstruction)
%ReconstructTwoForm2D computes the reconstruction of a discretized 2-form
%   in 2 dimensions. It performs the reconstruction for all the elements
%   given.
%
%   Returns a matrix containing the values of the
%   reconstructed 2-form at the reconstructed nodes, following the same
%   structure as for the input of discreteZeroForm.
%
%   USAGE
%   -----
%       [reconsctructed,xiRefinedGrid,etaRefinedGrid] = 
%                       ReconstructTwoForm2D(discreteTwoForm,gridType,mesh,nReconstruction)
%
%   INPUTS
%   ------
%       discreteTwoForm :: The 2-form discretized in each element.
%                          discreteZeroForm(i,j) is the coefficient
%                          associated to the ith basis function of the jth
%                          element. Therefore columns contain the
%                          coefficients in each element. Since we are in
%                          2D, for spectral elements of order p, we have
%                          p basis in each direction, therefore there
%                          are p*p basis functions (coefficients).
%                          (type: float64, size: [p*p,n])
%       gridType :: The type of grid. Can be:
%                'Lobatto', 'EGauss' or 'Gauss'
%                (type: string, size: single value)
%       mesh :: (type: mimeticFEM2.Mesh, size: single value)
%               A mesh object over which the 2-form is defined.
%       nReconstruction :: The number of points to use in the x and y
%                          direction for the plotting. As an alternative
%                          a list of points xi\in[-1,1] can be given.
%                          (type: int64, shape: single value, [1,N] or [N,1])
%
%   OUTPUTS
%   -------  
%       reconstructed :: The reconstructed 2-form at the reconstruction
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
%   Copyright 2009-2015 Artur Palha

%   Revisions:  2017-22-02 (apalha) Inputs changed. Mesh object is now an
%                                   input, instead of g.
%               2015-06-09 (apalha) Updated documentation.
%               2009-11-06 (apalha) Added the possibility to reconstruct for
%                                   more than one element.
%               2009-11-05 (apalha) Added the possibility to reconstruct at 
%                                   a given set of parametric points.
%               2009-11-05 (apalha) Created the function for reconstruction 
%                                   of one element.
%               2009-11-04 (apalha) First implementation.              
%               
    
    % the number of elements
    nElements = mesh.numElements;

    % compute the order of the discretization
    p = sqrt(size(discreteTwoForm,1));
    
    if iscell(nReconstruction)
        xiRefinedGrid = nReconstruction{1}(:);
        etaRefinedGrid = nReconstruction{2}(:);
        nReconstructed = length(xiRefinedGrid);
        
        % compute the basis functions at the refined nodes
        xiBasisRefined = mimeticFEM.EdgePoly(xiRefinedGrid,p,gridType);
        etaBasisRefined = mimeticFEM.EdgePoly(etaRefinedGrid,p,gridType);
    
        xietaBasisKron = zeros(p*p,nReconstructed);
        for kxi=1:p
            for keta=1:p
                xietaBasisKron((kxi-1)*p + keta,:) = xiBasisRefined(kxi,:).*etaBasisRefined(keta,:);
            end
        end
        xietaBasisKron = xietaBasisKron';
        
    elseif isscalar(nReconstruction)
        % compute the local nodes where to compute the refinement
        xiRefined = (2/(nReconstruction-1))*(0:(nReconstruction-1))-1;
        etaRefined = xiRefined;
        nReconstructed = nReconstruction*nReconstruction;
        [xiRefinedGrid,etaRefinedGrid] = meshgrid(xiRefined,etaRefined);
        xiRefinedGrid = xiRefinedGrid(:);
        etaRefinedGrid = etaRefinedGrid(:);
    
        % compute the basis functions at the refined nodes
        xiBasisRefined = mimeticFEM.EdgePoly(xiRefined,p,gridType);
        etaBasisRefined = xiBasisRefined;
    
        xietaBasisKron = kron(xiBasisRefined,etaBasisRefined)';
        
    elseif isvector(nReconstruction)
        xiRefined = nReconstruction;
        etaRefined = xiRefined;
        nReconstruction = length(nReconstruction);
        nReconstructed = nReconstruction*nReconstruction;
        [xiRefinedGrid,etaRefinedGrid] = meshgrid(xiRefined,etaRefined);
        xiRefinedGrid = xiRefinedGrid(:);
        etaRefinedGrid = etaRefinedGrid(:);
    
        % compute the basis functions at the refined nodes
        xiBasisRefined = mimeticFEM.EdgePoly(xiRefined,p,gridType);
        etaBasisRefined = xiBasisRefined;
    
        xietaBasisKron = kron(xiBasisRefined,etaBasisRefined)';
    end
    
    
    % allocate memory space for reconstructed result
    reconstructed = zeros([nReconstructed nElements]);
    
    if nElements>1
        for element=1:nElements
            % compute the metric at the refined points
            gEvaluatedMatrix = spdiags(1./mesh.g(element,xiRefinedGrid(:),etaRefinedGrid(:)),0,nReconstructed,nReconstructed);

            % the combined basis in 2d already scaled with the metric
            xietaBasis = gEvaluatedMatrix*xietaBasisKron;

            % reconstruct
            reconstructed(:,element) = xietaBasis*discreteTwoForm(:,element);
        end
    else
        % compute the metric at the refined points
        gEvaluatedMatrix = spdiags(1./mesh.g(1,xiRefinedGrid(:),etaRefinedGrid(:)),0,nReconstructed,nReconstructed);

        % the combined basis in 2d already scaled with the metric
        xietaBasis = gEvaluatedMatrix*xietaBasisKron;

        % reconstruct
        reconstructed(:) = xietaBasis*discreteTwoForm(:);
    end
        
end
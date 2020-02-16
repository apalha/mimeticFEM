function [reconstructed, xiRefinedGrid, etaRefinedGrid] = ReconstructZeroForm2D(discreteZeroForm,gridType,nReconstruction)
%ReconstructZeroForm2D computes the reconstruction of a discretized 0-form
%   in 2 dimensions. It performs the reconstruction for all the elements
%   given.
%
%   Returns a matrix containing the values of the
%   reconstructed 0-form at the reconstructed nodes, following the same
%   structure as for the input of discreteZeroForm.
%
%   USAGE
%   -----
%       [reconsctructed,xiRefinedGrid,etaRefinedGrid] = 
%                       ReconstructZeroForm2D(discreteZeroForm,gridType,nReconstruction)
%
%   INPUTS
%   ------
%       discreteZeroForm :: The 0-form discretized in each element.
%                           discreteZeroForm(i,j) is the coefficient
%                           associated to the ith basis function of the jth
%                           element. Therefore columns contain the
%                           coefficients in each element. Since we are in
%                           2D, for spectral elements of order p, we have
%                           (p+1) basis in each direction, therefore there
%                           are (p+1)(p+1) basis functions (coefficients).
%                           (type: float64, size: [(p+1)(p+1),n])
%       gridType :: The type of grid. Can be:
%                'Lobatto', 'EGauss' or 'Gauss'
%                (type: string, size: single value)
%       nReconstruction :: The number of points to use in the x and y
%                          direction for the plotting. As an alternative
%                          a list of points xi\in[-1,1] can be given.
%                          (type: int64, shape: single value, [1,N] or [N,1])
%
%   OUTPUTS
%   -------  
%       reconstructed :: The reconstructed 0-form at the reconstruction
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

    
    % % the number of elements
    % n = size(discreteZeroForm,2);

    % compute the order of the discretization
    p = sqrt(size(discreteZeroForm,1))-1;
    
    if iscell(nReconstruction)
        xiRefinedGrid = nReconstruction{1}(:);
        etaRefinedGrid = nReconstruction{2}(:);
        nReconstructed = length(xiRefinedGrid);
            
        % compute the basis functions at the refined nodes
        if strcmp(gridType,'Lobatto')
            xiBasisRefined = mimeticFEM.LobattoPoly(xiRefinedGrid,p);
            etaBasisRefined = mimeticFEM.LobattoPoly(etaRefinedGrid,p);
        elseif strcmp(gridType,'Gauss')
            xiBasisRefined = mimeticFEM.GaussPoly(xiRefinedGrid,p);
            etaBasisRefined = mimeticFEM.GaussPoly(etaRefinedGrid,p);
        end
        
        xietaBasisKron = zeros((p+1)*(p+1),nReconstructed);
        for kxi=1:(p+1)
            for keta=1:(p+1)
                xietaBasisKron((kxi-1)*(p+1) + keta,:) = xiBasisRefined(kxi,:).*etaBasisRefined(keta,:);
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
        if strcmp(gridType,'Lobatto')
            xiBasisRefined = mimeticFEM.LobattoPoly(xiRefined,p);
        elseif strcmp(gridType,'Gauss')
            xiBasisRefined = mimeticFEM.GaussPoly(xiRefined,p);
        end
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
        if strcmp(gridType,'Lobatto')
            xiBasisRefined = mimeticFEM.LobattoPoly(xiRefined,p);
        elseif strcmp(gridType,'Gauss')
            xiBasisRefined = mimeticFEM.GaussPoly(xiRefined,p);
        end
        etaBasisRefined = xiBasisRefined;
    
        xietaBasisKron = kron(xiBasisRefined,etaBasisRefined)';
    end
    
    
    
    % allocate memory space for reconstructed result
    %reconstructed = zeros([nReconstruction*nReconstruction n]);
    
    % reconstruct
    reconstructed = xietaBasisKron*discreteZeroForm;
        
end
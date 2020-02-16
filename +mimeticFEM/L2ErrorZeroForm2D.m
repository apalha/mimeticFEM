function [globalError, localError, nodesX, nodesY] = L2ErrorZeroForm2D(discreteZeroForm,analyticalZeroForm,mesh,nReconstruction,pErrorInt,gridType,figureNumber,varargin)
%L2ErrorZeroForm2D Computes the error of a 0-form.
%
%   L2ErrorZeroForm2D(discreteZeroForm,analyticalTwoForm,...
%                    mesh,nReconstruction,gridType,figureNumber)
%   Where:
%
%   INPUTS:
%       discreteZeroForm   :: the 0-form discretized as a vector: f_{1,1},
%                            f_{1,2}, f_{1,3}, ..., f_{2,1},... on all
%                            elements
%       analyticalZeroform :: the analytical 0-form
%       mesh :: (type: mimeticFEM2.Mesh, size: single value)
%               A mesh object over which the 2-form is defined.
%       nPointsRefinement :: the number of points to use in the x and y
%                            direction for the refinement to plot the
%                            error.
%       pErrorInt         :: the order of the quadrature used to compute
%                            the global error.
%       gridType          :: the type of grid.
%       figureNumber      :: plot or not the local error if
%                            figureNumber==0, the plot is not done,
%                            otherwise the plot is made to figure numbered
%                            figureNumber.
%
%
%   OUTPUTS:
%
%       globalError     :: the total L2 error, the integral
%       localError      :: the L2 error at each point
%       nodesX          :: the x coordinates of the nodes where the L2 error
%                          is computed
%       nodesY          :: the y coordinates of the nodes where the L2 error
%                          is computed
%
%   Copyright 2017 Artur Palha
%   $Revision: 1.0 $  $Date: 2017/05/23 $
    
    % do the optional tests
    if any(strcmp('PlotType',varargin))
        optIndex = find(strcmp('PlotType',varargin));
        PlotType = varargin{optIndex+1};
    else
        PlotType = 'surf';
    end
    
    nElements = mesh.numElements;
    
    %% Compute the global error
    
    % compute the integration nodes of Gauss quadrature
    [intErrorNodes,quadWeights] = mimeticFEM.GaussQuad(pErrorInt);
    quadWeightsGrid = kron(quadWeights,quadWeights');
    
    % Reconstruct the 0-form
    %[reconstructedForm,xiRefinedGrid,etaRefinedGrid] = mimeticFEM.ReconstructTwoForm2D(discreteTwoForm,gridType,g,intErrorNodes);
    [xiIntErrorNodes,etaIntErrorNodes] = meshgrid(intErrorNodes,intErrorNodes);
    [reconstructedForm,xiRefinedGrid,etaRefinedGrid] = mimeticFEM2.ReconstructZeroForm2D(discreteZeroForm,gridType,{xiIntErrorNodes(:),etaIntErrorNodes(:)});
    
    globalError = 0;
    
    % loop over the elements
    for element = 1:nElements
        if mesh.dim(2) == 2
            % transform the local coordinates into physical coordinates
            [xGrid,yGrid] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);

            % compute the analytical 1-form at the physical coordinates
            zeroFormAnalytical = analyticalZeroForm(xGrid,yGrid);
        elseif mesh.dim(2) == 3
            % transform the local coordinates into physical coordinates
            [xGrid,yGrid,zGrid] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);

            % compute the analytical 1-form at the physical coordinates
            zeroFormAnalytical = analyticalZeroForm(xGrid,yGrid,zGrid);
        end
        
        evaluatedg = abs(mesh.g(element,xiRefinedGrid,etaRefinedGrid));
        
        % compute the global error at the element
        globalError = globalError + quadWeightsGrid(:)' * ((((reconstructedForm(:,element)-zeroFormAnalytical).^2)).*evaluatedg(:));
    end
    
    globalError = sqrt(globalError);
    
    %% Compute the local error
    
    % Reconstruct the 1-form
    [reconstructedForm,xiRefinedGrid,etaRefinedGrid] = mimeticFEM2.ReconstructZeroForm2D(discreteZeroForm,gridType,nReconstruction);
        
    localError = zeros(size(reconstructedForm));
    
    nodesX = zeros(size(reconstructedForm));
    nodesY = zeros(size(reconstructedForm));
    if mesh.dim(2) == 3
        nodesZ = zeros(size(reconstructedForm));
    end
    
    if(figureNumber>0)
        figure(figureNumber)
        axis([mesh.xBounds mesh.yBounds]);
        hold on
    end
    
    % loop over the elements
    for element = 1:nElements
        if mesh.dim(2) == 2
            % transform the local coordinates into physical coordinates
            [nodesX(:,element),nodesY(:,element)] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);
            % compute the analytical 1-form at the physical coordinates
            zeroFormAnalytical = analyticalZeroForm(nodesX(:,element),nodesY(:,element));
        elseif mesh.dim(2) == 3
            % transform the local coordinates into physical coordinates
            [nodesX(:,element),nodesY(:,element), nodesZ(:, element)] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);
            % compute the analytical 1-form at the physical coordinates
            zeroFormAnalytical = analyticalZeroForm(nodesX(:,element),nodesY(:,element),nodesZ(:,element));
        end
        
        

        % compute the global error at the element
        localError(:,element) = abs(reconstructedForm(:,element)-zeroFormAnalytical);
        
        if(figureNumber>0)
            if strcmp('surf',PlotType)
                if mesh.dim(2) == 2
                    surf(reshape(nodesX(:,element),nReconstruction,nReconstruction),reshape(nodesY(:,element),nReconstruction,nReconstruction),log10(reshape(localError(:,element),nReconstruction,nReconstruction)),'EdgeColor','None')
                elseif mesh.dim(2) == 3
                    surf(reshape(nodesX(:,element),nReconstruction,nReconstruction),reshape(nodesY(:,element),nReconstruction,nReconstruction),reshape(nodesZ(:,element),nReconstruction,nReconstruction),log10(reshape(localError(:,element),nReconstruction,nReconstruction)),'EdgeColor','None')
                end
            elseif strcmp('pcolor',PlotType)
                pcolor(reshape(nodesX(:,element),nReconstruction,nReconstruction),reshape(nodesY(:,element),nReconstruction,nReconstruction),log10(reshape(localError(:,element),nReconstruction,nReconstruction)))
            end
        end
    end
    
    if(figureNumber>0)
        shading interp
        axis square
        box on
    end

end
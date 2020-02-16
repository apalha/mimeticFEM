function [globalError, localError, nodesX, nodesY] = L2ErrorOneForm2D(discreteOneForm,analyticalOneForm,mesh,nReconstruction,pErrorInt,gridType,n,figureNumber,Hodge)
%L2Error2D Plots a discrete OneForm.
%
%   [globalError, localError, nodesX, nodesY] = L2ErrorOneForm2D(discreteOneForm,analyticalOneForm,mesh,nReconstruction,pErrorInt,gridType,n,figureNumber,Hodge)
%   Where:
%
%   INPUTS:
%       discreteOneForm   :: the 1-form discretized as a vector: f_{1,1},
%                            f_{1,2}, f_{1,3}, ..., f_{2,1},... on all
%                            elements
%       analyticalOneform :: the analytical 1-form as flexoutput
%       mesh :: (type: mimeticFEM2.Mesh, size: single value)
%               A mesh object over which the 2-form is defined.
%       xBounds           :: the x bounds of the physical domain
%       yBounds           :: the y bounds of the physical domain
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
%       Hodge             :: Flag that specifies if the error is computed
%                            with the current form (false) of with the
%                            Hodge (true)
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
%   Copyright 2011 Artur Palha
%   $Revision: 1.0 $  $Date: 2011/12/12 $

    nElements = n(1)*n(2);
    
    %% Compute the global error
    
    % compute the integration nodes of Gauss quadrature
    [intErrorNodes, quadWeights] = mimeticFEM.GaussQuad(pErrorInt);
    quadWeightsGrid = kron(quadWeights,quadWeights');
    
    % Reconstruct the 1-form
    [reconstructedX, reconstructedY, xiRefinedGrid, etaRefinedGrid] = mimeticFEM2.ReconstructOneForm2D(discreteOneForm,gridType,intErrorNodes,mesh);
    globalError = 0;
    
    if Hodge % compute the Hodge-* of the reconstructed form
        reconstructedTemp = - reconstructedX;
        reconstructedX = reconstructedY;
        reconstructedY = reconstructedTemp;
    end
    
    % loop over the elements
    for element = 1:nElements
        % transform the local coordinates into physical coordinates
        [xGrid, yGrid] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);
        
        % compute the analytical 1-form at the physical coordinates
        [xAnalytical, yAnalytical] = analyticalOneForm(xGrid,yGrid);
        
        evaluatedg = mesh.g(element,xiRefinedGrid,etaRefinedGrid);
        
        % compute the global error at the element
        globalError = globalError + quadWeightsGrid(:)' * ((((reconstructedX(:,element)-xAnalytical).^2)+((reconstructedY(:,element)-yAnalytical).^2)).*evaluatedg(:));
    end
    
    globalError = sqrt(globalError);
    
    %% Compute the local error
    
    if(figureNumber>0)
        % Reconstruct the 1-form
        [reconstructedX, reconstructedY, xiRefinedGrid, etaRefinedGrid] = mimeticFEM.ReconstructOneForm2D(discreteOneForm,gridType,nReconstruction,mesh);

        if Hodge % compute the Hodge-* of the reconstructed form
            reconstructedTemp = - reconstructedX;
            reconstructedX = reconstructedY;
            reconstructedY = reconstructedTemp;
        end

        localError = zeros(size(reconstructedX));

        nodesX = zeros(size(reconstructedX));
        nodesY = zeros(size(reconstructedX));
    
        figure(figureNumber)
        axis([xBounds yBounds]);
        hold on
    
        % loop over the elements
        for element = 1:nElements
            % transform the local coordinates into physical coordinates
            [nodesX(:,element), nodesY(:,element)] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);

            % compute the analytical 1-form at the physical coordinates
            [xAnalytical, yAnalytical] = analyticalOneForm(nodesX(:,element),nodesY(:,element));

            % compute the global error at the element
            localError(:,element) = ((((reconstructedX(:,element)-xAnalytical).^2)+((reconstructedY(:,element)-yAnalytical).^2)));

            if(figureNumber>0)
                surf(reshape(nodesX(:,element),nReconstruction,nReconstruction),reshape(nodesY(:,element),nReconstruction,nReconstruction),log10(reshape(localError(:,element),nReconstruction,nReconstruction)),'EdgeColor','None')
            end
        end

        shading interp
        axis square
        box on
    
    end
end
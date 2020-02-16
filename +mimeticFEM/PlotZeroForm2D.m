function PlotZeroForm2D(discreteZeroForm,mesh,nReconstruction,gridType,figureNumber,varargin)
%PlotZeroForm2D Plots a discrete ZeroForm.
%
%   USAGE
%   -----
%       PlotZeroForm2D(discreteZeroForm,mesh,
%                      nReconstruction,gridType,figureNumber)
%
%       PlotZeroForm2D(discreteZeroForm,mesh,
%                      nReconstruction,gridType,figureNumber,
%                      'IsSubPlot',subPlotNumber,
%                      'NormalizePlot',normalizedMaxValue)
%
%   INPUTS
%   ------
%
%       discreteZeroForm :: The 0-form discretized in each element.
%                           discreteZeroForm(i,j) is the coefficient
%                           associated to the ith basis function of the jth
%                           element. Therefore columns contain the
%                           coefficients in each element. Since we are in
%                           2D, for spectral elements of order p, we have
%                           (p+1) basis in each direction, therefore there
%                           are (p+1)(p+1) basis functions (coefficients).
%                           (type: float64, size: [(p+1)(p+1),n])
%       mesh :: The mesh over which the TwoForm is defined.
%               (type: mimeticFEM.Mesh, size: single value)
%       nReconstruction :: The number of points to use in the x and y
%                          direction for the plotting. As an alternative
%                          a list of points xi\in[-1,1] can be given.
%                          (type: int64, shape: single value, [1,N] or [N,1])
%       gridType :: The type of grid. Can be:
%                'Lobatto', 'EGauss' or 'Gauss'
%                (type: string, size: single value)
%       figureNumber :: The number for the figure where to plot. If
%                       the figure with number figureNumber does not
%                       exist, a new figure is generated.
%                       (type: int32, size: single value)
%       subPlotNumber :: If the figure is to be generated in a subplot, use
%       -------------    this value as the subplot number. Default is [].
%                        (type: int32, size: single value)
%       normalizedMaxValue :: Use this value to normalize the maximum value
%       ------------------    of the plot. Default is [].
%                             (type: float64, size: single value)
%
%   OUTPUTS
%   -------  
%
%
%   Copyright 2011 Artur Palha

%   Revisions:  2011-11-25 (apalha) First implementation.
    
    % the discrete zero form contains the degrees of freedom within an
    % element in the columns, therefore there are as many elements as
    % columns in discreteZeroForm
    nElements = mesh.numElements;
    
    figure(figureNumber)
    
    % do the optional tests
    if any(strcmp('IsSubPlot',varargin))
        optIndex = find(strcmp('IsSubPlot',varargin));
        subplot(varargin{optIndex+1}(1),varargin{optIndex+1}(2),varargin{optIndex+1}(3))
    else
        %clf(figureNumber)
    end
    
    if any(strcmp('NormalizePlot',varargin))
        optIndex = find(strcmp('NormalizePlot',varargin));
        normalizedMaxValue = varargin{optIndex+1};
    else
        normalizedMaxValue = [];
    end
    
    if mesh.dim(2) == 3
        axis([mesh.xBounds mesh.yBounds]);
        hold on

    %    if nElements>1
    %        if length(discreteTwoForm) == size(discreteTwoForm,1)
    %            discreteTwoForm = discreteTwoForm';
    %        end
    %    end

        if nElements == 1
            for element = 1:nElements
                [reconstructed, xiRefinedGrid, etaRefinedGrid] = mimeticFEM.ReconstructZeroForm2D(discreteZeroForm(:,element),gridType,nReconstruction);
                [xGrid, yGrid, zGrid] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);
                surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(zGrid,nReconstruction,nReconstruction),reshape(reconstructed,nReconstruction,nReconstruction),'EdgeColor','None')
            end
        else
            [reconstructed, xiRefinedGrid, etaRefinedGrid] = mimeticFEM.ReconstructZeroForm2D(discreteZeroForm,gridType,nReconstruction);
            maxValue = max(reconstructed(:));
            minValue = min(reconstructed(:));
            if ~isempty(normalizedMaxValue)
                reconstructed = reconstructed * (normalizedMaxValue/max(abs([minValue maxValue])));
            end

            for element = 1:nElements
                %[reconstructed xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm(:,element),g{element},nReconstruction,gridType);
                [xGrid, yGrid, zGrid] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);
                surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(zGrid,nReconstruction,nReconstruction),reshape(reconstructed(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            end
        end
        view(37.5,70)
    else
        axis([mesh.xBounds mesh.yBounds]);
        hold on

    %    if nElements>1
    %        if length(discreteTwoForm) == size(discreteTwoForm,1)
    %            discreteTwoForm = discreteTwoForm';
    %        end
    %    end

        if nElements == 1
            for element = 1:nElements
                [reconstructed, xiRefinedGrid, etaRefinedGrid] = mimeticFEM.ReconstructZeroForm2D(discreteZeroForm(:,element),gridType,nReconstruction);
                [xGrid, yGrid] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);
                surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructed,nReconstruction,nReconstruction),'EdgeColor','None')
            end
        else
            [reconstructed, xiRefinedGrid, etaRefinedGrid] = mimeticFEM.ReconstructZeroForm2D(discreteZeroForm,gridType,nReconstruction);
            maxValue = max(reconstructed(:));
            minValue = min(reconstructed(:));
            if ~isempty(normalizedMaxValue)
                reconstructed = reconstructed * (normalizedMaxValue/max(abs([minValue maxValue])));
            end

            for element = 1:nElements
                %[reconstructed xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm(:,element),g{element},nReconstruction,gridType);
                [xGrid, yGrid] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);
                surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructed(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            end
        end
    end
    
    shading interp
end
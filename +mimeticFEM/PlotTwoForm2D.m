function PlotTwoForm2D(discreteTwoForm,mesh,nReconstruction,gridType,figureNumber,varargin)
%PlotTwoForm2D Plots a discrete TwoForm.
%
%   USAGE
%   -----
%       PlotTwoForm2D(discreteTwoForm,mesh,
%                      nReconstruction,gridType,figureNumber)
%
%       PlotTwoForm2D(discreteZeroForm,mesh,
%                      nReconstruction,gridType,figureNumber,
%                      'IsSubPlot',subPlotNumber,
%                      'NormalizePlot',normalizedMaxValue,
%                      'HigherDimension',higherDimension)
%
%   INPUTS
%   ------
%
%       discreteTwoForm :: The 2-form discretized in each element.
%                          discreteTwoForm(i,j) is the coefficient
%                          associated to the ith basis function of the jth
%                          element. Therefore columns contain the
%                          coefficients in each element. Since we are in
%                          2D, for spectral elements of order p, we have
%                          p basis in each direction, therefore there
%                          are p*p basis functions (coefficients).
%                          (type: float64, size: [p p n])
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
%   Copyright 2011-2015 Artur Palha

%   Revisions:  2015-06-09 (apalha) Updated documentation.
%               2011-11-25 (apalha) First implementation.

    % the discrete two form contains the degrees of freedom within an
    % element in the columns, therefore there are as many elements as
    % columns in discreteTwoForm
    n = size(discreteTwoForm,2);
    
    figure(figureNumber)
    
    % do the optional tests
    if any(strcmp('IsSubPlot',varargin))
        optIndex = find(strcmp('IsSubPlot',varargin));
        subplot(varargin{optIndex+1}(1),varargin{optIndex+1}(2),varargin{optIndex+1}(3))
    else
        clf(figureNumber)
    end
    
    if any(strcmp('NormalizePlot',varargin))
        optIndex = find(strcmp('NormalizePlot',varargin));
        normalizedMaxValue = varargin{optIndex+1};
    else
        normalizedMaxValue = [];
    end
    
    if any(strcmp('HigherDimension',varargin))
        optIndex = find(strcmp('HigherDimension',varargin));
        higherDimension = varargin{optIndex+1};
    else
        higherDimension = false;
    end
    
    if any(strcmp('PlotType',varargin))
        optIndex = find(strcmp('PlotType',varargin));
        PlotType = varargin{optIndex+1};
    else
        PlotType = 'surf';
    end
    
    if any(strcmp('ContourLines',varargin))
        optIndex = find(strcmp('ContourLines',varargin));
        ContourLines = varargin{optIndex+1};
    else
        ContourLines = linspace(0,1,11);
    end
    
    if mesh.dim(2) == 3
        axis([mesh.xBounds mesh.yBounds]);
        hold on

        if n == 1
            for element = 1:n
                [reconstructed,xiRefinedGrid,etaRefinedGrid] =  mimeticFEM2.ReconstructTwoForm2D(discreteTwoForm(:,element),gridType,mesh,nReconstruction);
                [xGrid,yGrid,zGrid] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);
                surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(zGrid,nReconstruction,nReconstruction),reshape(reconstructed,nReconstruction,nReconstruction),'EdgeColor','None')
            end
        else
            [reconstructed,xiRefinedGrid,etaRefinedGrid] =  mimeticFEM2.ReconstructTwoForm2D(discreteTwoForm,gridType,mesh,nReconstruction);
            maxValue = max(reconstructed(:));
            minValue = min(reconstructed(:));
            if ~isempty(normalizedMaxValue)
                reconstructed = reconstructed * (normalizedMaxValue/max(abs([minValue maxValue])));
            end

            for element = 1:n
                %[reconstructed xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm(:,element),g{element},nReconstruction,gridType);
                [xGrid,yGrid,zGrid] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);
                surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(zGrid,nReconstruction,nReconstruction),reshape(reconstructed(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
            end
        end
        view(37.5,70)
    else
        axis([mesh.xBounds mesh.yBounds]);
        hold on

        if n == 1
            for element = 1:n
                [reconstructed,xiRefinedGrid,etaRefinedGrid] = mimeticFEM2.ReconstructTwoForm2D(discreteTwoForm(:,element),gridType,mesh,nReconstruction);
                [xGrid,yGrid] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);
                if strcmp(PlotType,'surf')
                    surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructed,nReconstruction,nReconstruction),'EdgeColor','None')
                elseif strcmp(PlotType,'contour')
                    contour(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructed,nReconstruction,nReconstruction),ContourLines)
                end
                    
            end
        else
            [reconstructed,xiRefinedGrid,etaRefinedGrid] = mimeticFEM2.ReconstructTwoForm2D(discreteTwoForm,gridType,mesh,nReconstruction);
            maxValue = max(reconstructed(:));
            minValue = min(reconstructed(:));
            if ~isempty(normalizedMaxValue)
                reconstructed = reconstructed * (normalizedMaxValue/max(abs([minValue maxValue])));
            end

            for element = 1:n
                %[reconstructed xiRefinedGrid etaRefinedGrid] = ReconstructTwoForm2D(discreteTwoForm(:,element),g{element},nReconstruction,gridType);
                [xGrid,yGrid] = mesh.mapping(element,xiRefinedGrid,etaRefinedGrid);
                if strcmp(PlotType,'surf')
                    surf(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructed(:,element),nReconstruction,nReconstruction),'EdgeColor','None')
                elseif strcmp(PlotType,'contour')
                    contour(reshape(xGrid,nReconstruction,nReconstruction),reshape(yGrid,nReconstruction,nReconstruction),reshape(reconstructed(:,element),nReconstruction,nReconstruction),ContourLines)
                end
            end
        end
    end
    if strcmp(PlotType,'surf')
        shading interp
    end
end
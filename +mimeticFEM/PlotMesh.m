function PlotMesh(mesh,p,gridType,nPlotPoints,figureNumber,displayElementNumber)
%PlotMesh plots the mesh.
%
%   USAGE
%   -----
%       PlotMesh(mesh,figureNumber)
%
%   INPUTS
%   ------
%       mesh :: The mesh to plot.
%               (type: mimeticFEM2.Mesh, size: single value)
%       p    :: The polynomial degree of the mesh.
%               (type: int32, size: single value)
%       gridType :: The type of grid: 'Lobatto' or 'EGauss'.
%                   (type: string, size: single value)
%       nPlotPoints :: The number of points used to plot each edge.
%                      (type: int32, size: single value)
%       figureNumber :: A figure number for an already existing figure. If []
%                    is given then a new figure is generated.
%                    (type: int32, size: single value)
%       displayElementNumber :: Flag to display the element's number in the
%                               plot.
%                               (type: logic, size: single value)
%
%   OUTPUTS
%   -------  
%
%   Copyright 2009 Artur Palha

%   Revisions:  2009-11-25 (apalha) First implementation.
%               2014-12-05 (apalha) Replaced myAxis where the x and y bounds
%                                   were given together, by xBounds and
%                                   yBounds. This was done to make the
%                                   inputs more homogeneous between
%                                   functions.
    
    if ~isempty(figureNumber) % use the figure handle given by the user
        figure(figureNumber)
    else % or create a new figure
        figure()
    end
    
    hold on
    
    % check if gridType is a valid one
    if ~mimeticFEM.TestPolyType(gridType)
        disp(sprintf(':: %s :: is not a valid type of grid', gridType));
        return
    end
    
    % generate the fine points to plot
    s = linspace(-1,1,nPlotPoints)';
    
    % plot the higher order lines
    if strcmp(gridType,'Lobatto')
        rNodes = mimeticFEM.LobattoQuad(p);
    else 
        rNodes = mimeticFEM.EGaussQuad(p);
    end
    
    if mesh.dim(2) == 2
        % Bottom
        [x, y] = mesh.mapping([],s,-ones(size(s)));
        plot(reshape(x,[],mesh.numElements),reshape(y,[],mesh.numElements),'k-','LineWidth',1.5)
        % Top
        [x, y] = mesh.mapping([],s,ones(size(s)));
        plot(reshape(x,[],mesh.numElements),reshape(y,[],mesh.numElements),'k-','LineWidth',1.5)
        % Left
        [x, y] = mesh.mapping([],-ones(size(s)),s);
        plot(reshape(x,[],mesh.numElements),reshape(y,[],mesh.numElements),'k-','LineWidth',1.5)
        % Right
        [x, y] = mesh.mapping([], ones(size(s)),s);
        plot(reshape(x,[],mesh.numElements),reshape(y,[],mesh.numElements),'k-','LineWidth',1.5)
    
    elseif mesh.dim(2) == 3
        % Bottom
        [x, y, z] = mesh.mapping([],s,-ones(size(s)));
        plot3(reshape(x,[],mesh.numElements),...
             reshape(y,[],mesh.numElements),...
             reshape(z,[],mesh.numElements), 'k-','LineWidth',1.5)
        % Top
        [x, y, z] = mesh.mapping([],s,ones(size(s)));
        plot3(reshape(x,[],mesh.numElements),...
             reshape(y,[],mesh.numElements),...
             reshape(z,[],mesh.numElements), 'k-','LineWidth',1.5)
        % Left
        [x, y, z] = mesh.mapping([],-ones(size(s)),s);
        plot3(reshape(x,[],mesh.numElements),...
             reshape(y,[],mesh.numElements),...
             reshape(z,[],mesh.numElements), 'k-','LineWidth',1.5)
        % Right
        [x, y, z] = mesh.mapping([], ones(size(s)),s);
        plot3(reshape(x,[],mesh.numElements),...
             reshape(y,[],mesh.numElements),...
             reshape(z,[],mesh.numElements), 'k-','LineWidth',1.5)
    end
    
    color = [0.3,0.3,0.3]; % grey
    
    if mesh.dim(2) == 2
        for pLine=2:p
            % horizontal lines
            [x, y] = mesh.mapping([],s,rNodes(pLine)*ones(size(s)));
            plot(reshape(x,[],mesh.numElements),reshape(y,[],mesh.numElements),'-','Color',color)
            % vertical lines
            [x, y] = mesh.mapping([],rNodes(pLine)*ones(size(s)),s);
            plot(reshape(x,[],mesh.numElements),reshape(y,[],mesh.numElements),'-','Color',color)
        end
    elseif mesh.dim(2) == 3
        for pLine=2:p
            % horizontal lines
            [x, y, z] = mesh.mapping([],s,rNodes(pLine)*ones(size(s)));
            plot3(reshape(x,[],mesh.numElements),...
                  reshape(y,[],mesh.numElements),...
                  reshape(z,[],mesh.numElements),...
                  '-','Color',color)
            % vertical lines
            [x, y, z] = mesh.mapping([],rNodes(pLine)*ones(size(s)),s);
            plot3(reshape(x,[],mesh.numElements),...
                  reshape(y,[],mesh.numElements),...
                  reshape(z,[],mesh.numElements),...
                  '-','Color',color)
        end
    end
    
    % plot the grid
    if displayElementNumber
        if mesh.dim(2) == 2
            for element=1:mesh.numElements
                % place the elment number in the middle
                [xMid, yMid] = mesh.mapping(element,0,0);
                text(xMid,yMid,...
                   sprintf('%d',element),'HorizontalAlignment','center','VerticalAlignment','middle','Fontsize',14,'FontWeight','Bold','FontName','Times')
            end
        elseif mesh.dim(2) == 3
            for element=1:mesh.numElements
                % place the elment number in the middle
                [xMid, yMid, zMid] = mesh.mapping(element,0,0);
                text(xMid,yMid,zMid,...
                   sprintf('%d',element),'HorizontalAlignment','center','VerticalAlignment','middle','Fontsize',14,'FontWeight','Bold','FontName','Times')
            end
        end
    end
    
    %plot([-1 1 1 -1 -1], [-1 -1 1 1 -1],'-k','LineWidth',3) % the outer boundary
    if ~isempty(mesh.xBounds)
        xlim(mesh.xBounds)
    end
    if ~isempty(mesh.yBounds)
        ylim(mesh.yBounds)
    end
end
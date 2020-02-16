classdef IrregularSquare %< mimeticFEM2.Mesh
    %IrregularSquare Structured irregular mesh of a square
    %
    % IrregularSquare Properties:
    %   type - (type: str, size: single value)
    %          The type of mesh, in this case 'CrazyMesh'.
    %   dim - (type: int, size: [2, 1])
    %         The spatial dimensions of the mesh. The spatial dimension of the
    %         canonical space and of the physical space, in this case [2, 3].
    %   numElements - (type: int, size: single value)
    %                 The number of elements of the mesh.
    % SurfaceOfSphere Methods:
    %   mapping - Computes the mapping between (xi, eta) and (x, y).
    %   dXdXi - Computes the \frac{\partial X}{\partial \xi} derivative of
    %           the mapping.
    %   dXdEta - Computes the \frac{\partial X}{\partial \eta} derivative
    %            of the mapping.
    %   dYdXi - Computes the \frac{\partial Y}{\partial \xi} derivative of
    %           the mapping.
    %   dYdEta - Computes the \frac{\partial Y}{\partial \eta} derivative
    %            of the mapping.
    %   g - Computes the determinant of the metric det(g_{ij}) induced by
    %       the mapping.
    %   g11 - Computes the g^{11} component of the metric induced by the
    %         mapping.
    %   g22 - Computes the g^{22} component of the metric induced by the
    %         mapping.
    %   g12 - Computes the g^{12} component of the metric induced by the
    %         mapping.
    %
    %
    % Copyright 2018 Artur Palha

    %   Revisions:  2018-07-07 First implementation. (apalha)
    %               

    
    
    properties
        % type - (str, single value)
        %        The type of mesh, in this case 'SurfaceOfSphere'.
        type
        % dim - (int, single value) 
        %       The spatial dimension of the mesh, in this case [2, 3].
        dim
        % numElements - (int, single value)
        %               The number of element of the mesh.
        numElements
        % nodesOfElement - (float64, [4, 2, numElements])
        %                  The 2D coordinates of the nodes of the mesh. These
        %                  coordinates are local coordinates to the
        %                  subdomain to which the element belongs to.
        nodesOfElements
        % subdomainOfElements - (float64, [numElements, 1])
        %                       The index of the subdomain to which the
        %                       element belongs to.
        subdomainOfElements
        n
        xBounds
        yBounds
    end
    
    methods
        function obj = IrregularSquare(n, randomizeNodes)
            obj.numElements = n*n;
            obj.n = [n, n];
            obj.dim = [2, 2];
            obj.nodesOfElements = obj.computeNodesOfElements(randomizeNodes);
            obj.xBounds = [];
            obj.yBounds = [];
        end
        
        function nodesOfElements = computeNodesOfElements(obj, randomizeNodes)
            
            [nodesCoordinatesX, nodesCoordinatesY] = meshgrid(linspace(-1, 1, obj.n(1)+1), linspace(-1, 1, obj.n(2)+1));

            if randomizeNodes > 0
%                 dx_random = 0.0;%randomizeNodes;
%                 dy_random = randomizeNodes;
%                 nodesCoordinatesX(2:end-1, 2:end-1) = nodesCoordinatesX(2:end-1, 2:end-1) + dx_random*(2*rand(size(nodesCoordinatesX(2:end-1, 2:end-1))) - 1.0) * 2.0/(obj.n(1));
%                 nodesCoordinatesY(2:end-1, 2:end-1) = nodesCoordinatesY(2:end-1, 2:end-1) + dy_random*(2*rand(size(nodesCoordinatesY(2:end-1, 2:end-1))) - 1.0) * 2.0/(obj.n(2));
                iDeformed = 2:2:obj.n(1);
                jDeformed = 2:2:obj.n(2);
                nodesCoordinatesY(iDeformed, jDeformed) = nodesCoordinatesY(iDeformed, jDeformed) + randomizeNodes * 2.0/(obj.n(2));
                jDeformed = 1:2:obj.n(1);
                nodesCoordinatesY(iDeformed, jDeformed) = nodesCoordinatesY(iDeformed, jDeformed) - randomizeNodes * 2.0/(obj.n(2));
            end

            V1_x = reshape(nodesCoordinatesX(1:end-1, 1:end-1), 1, 1, []);
            V2_x = reshape(nodesCoordinatesX(1:end-1, 2:end), 1, 1, []);
            V3_x = reshape(nodesCoordinatesX(2:end, 2:end), 1, 1, []);
            V4_x = reshape(nodesCoordinatesX(2:end, 1:end-1), 1, 1, []);

            V1_y = reshape(nodesCoordinatesY(1:end-1, 1:end-1), 1, 1, []);
            V2_y = reshape(nodesCoordinatesY(1:end-1, 2:end), 1, 1, []);
            V3_y = reshape(nodesCoordinatesY(2:end, 2:end), 1, 1, []);
            V4_y = reshape(nodesCoordinatesY(2:end, 1:end-1), 1, 1, []);

            nodesOfElements = [V1_x, V1_y;...
                               V2_x, V2_y;...
                               V3_x, V3_y;...
                               V4_x, V4_y];
        end
    end    
end


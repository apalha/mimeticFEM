classdef CrazyMesh %< mimeticFEM2.Mesh
    %CrazyMesh Structured test wavy mesh (crazy mesh)
    %   Implements the wavy mesh (crazy mesh) as defined in (S4.1.1.4, [1])
    %
    % CrazyMesh Properties:
    %   type - (type: str, size: single value)
    %          The type of mesh, in this case 'CrazyMesh'.
    %   dim - (type: int, size: single value)
    %         The spatial dimension of the mesh, in this case 2.
    %   nodes - (type: double, size: [numNodes,2])
    %           An array whose columns are the x and y coordinates of the
    %           nodes of the mesh. 
    %   numNodes - (type: int, size: single value)
    %              The number of nodes of the mesh.
    %   numElements - (type: int, size: single value)
    %                 The number of element of the mesh.
    %   cc - (type: double, size: single value)
    %        Determines the curvature of the mesh (crazyness
    %        coefficient)
    %   xBounds - (type: double, size: [1,2])
    %             The bounds of the computational domain in the x-direction.
    %   yBounds - (type: double, size: [1,2])
    %             The bounds of the computational domain in the y-direction.
    %
    % CrazyMesh Methods:
    %   mapping - Computes the mapping between (xi,eta) and (x,y).
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
    %   elementBoundaryNodes - Computes n evenly spaced nodes over each of 
    %                          the four edges of each element
    %
    % 
    % [1] Palha, A., High order mimetic discretizations, PhD thesis.
    %
    % Copyright 2009 Artur Palha

    %   Revisions:  2017-21-02 First implementation. (apalha)
    %               

    
    
    properties
        % type - (str, single value)
        %        The type of mesh, in this case 'CrazyMesh'.
        type
        % dim - (type: int, size: [2, 1])
        %       The spatial dimensions of the mesh. The spatial dimension of the
        %       canonical space and of the physical space, in this case [2, 3].
        dim
        % nodes - (double, [numNodes,2])
        %         An array whose columns are the x and y coordinates of the
        %         nodes of the mesh. 
        nodes
        % numNodes - (int, single value)
        %            The number of nodes of the mesh.
        numNodes
        % numElements - (int, single value)
        %               The number of element of the mesh.
        numElements
        % cc - (double, single value)
        %      Determines the curvature of the mesh (crazyness
        %      coefficient). cc \in [0,0.2].
        %      0  : straight mesh.
        %      >0 : increasingly curved mesh.
        cc
        % xBounds - (type: double, size: [1,2])
        %           The bounds of the computational domain in the x-direction.
        xBounds
        % yBounds - (type: double, size: [1,2])
        %           The bounds of the computational domain in the y-direction.
        yBounds
        % n - (type: int, size: [1,2])
        %     The number of nodes in the x-direction (n(1)) and in the
        %     y-direction (n(2)).
        n
    end
    
    methods
        function obj = CrazyMesh(cc,n,xBounds,yBounds)
            obj.cc = cc;
            obj.numElements = prod(n);
            obj.xBounds = xBounds;
            obj.yBounds = yBounds;
            obj.n = n;
            obj.dim = [2, 2];
        end
    end
    
end


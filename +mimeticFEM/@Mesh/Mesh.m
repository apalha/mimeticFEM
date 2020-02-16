classdef (Abstract) Mesh
    %Mesh Abstract definition of Mesh.
    %   All meshes are defined as subclasses o this class.
    
    properties
        type
        dim
        nodes
        numNodes
        numElements
    end
    
    methods
        [x,y] = mapping(obj,element,xi,eta) % computes the mapping between (xi,eta) and (x,y)
        dXdXiEvaluated = dXdXi(obj,element,xi,eta) % computes the \frac{\partial X}{\partial \xi} derivative of the mapping
        dXdEtaEvaluated = dXdEta(obj,element,xi,eta) % computes the \frac{\partial X}{\partial \eta} derivative of the mapping
        dYdXiEvaluated = dYdXi(obj,element,xi,eta) % computes the \frac{\partial Y}{\partial \xi} derivative of the mapping
        dYdEtaEvaluated = dYdEta(obj,element,xi,eta) % computes the \frac{\partial Y}{\partial \eta} derivative of the mapping
        gEvaluated = g(obj,element,xi,eta) % computes the determinant of the metric det(g_{ij}) induced by the mapping
        g11Evluated = g11(obj,element,xi,eta) % computes the g^{11} component of the metric induced by the mapping
        g22Evaluated = g22(obj,element,xi,eta) % computes the g^{22} component of the metric induced by the mapping
        g12Evaluated = g12(obj,element,xi,eta) % computes the g^{12} component of the metric induced by the mapping
        nodesEvaluated = elementBoundaryNodes(obj,element,n) % computes n evenly spaced nodes over each of the four edges of each element
    end
    
end


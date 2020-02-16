function dXdEtaEvaluated = dXdEta(obj,elements,xi,eta)
% dXdEta computes the dXdEta term of the metric for
% at the local positions (xi,eta) for the element. For the IrregularSquare grid.
%
%   USAGE
%   -----
%       dXdEtaEvaluated = obj.dXdEta(elements,xi,eta)
%
%   INPUTS
%   ------
%       elements :: the elements numbers
%                   If elements is a scalar then x, y, and z have dimensions
%                   (N,M), if element is a vector of dimension Q then
%                   x, y, and z have dimensions (N,M,Q), where (N,M) are the
%                   dimensions of xi and eta. If element = [] then the
%                   mapping is computed for all elements
%                   (type: int32, size: vector)
%       xi :: the local (xi) coordinates of the nodes where to
%             evaluate the mapping.
%             (type: float64, size: array [N,M])
%       eta :: the local (eta) coordinates of the nodes where to
%              evaluate the mapping.
%              (type: float64, size: array [N,M])
%
%
%   OUTPUTS
%   -------
%       dXdEtaEvaluated  :: the evaluation of the eta derivative of the mapping
%                          between the computational and physical spaces at 
%                          the points (xi,eta) of the element.
%                          (type: float64, size: array [N,M] (Q=1) or [N,M,Q] (Q>1))
%
%   Copyright 2009-2018 Artur Palha

%   Revisions:  2018-07-07 First implementation. (apalha)
    
    % if element is empty [] then all the mapping is computed for all
    % elements
    if isempty(elements)
        elements = 1:obj.numElements;
    end
    
    % if elements is an array it must be reshaped so that the x and y
    % coordinates are an array with three indices (i,j,k). i,j correspond
    % to the indices of the xi and eta arrays and the k index to the
    % elements.
    if ~isscalar(elements)
        elements = reshape(elements, 1, 1, []);
    end
        
    %% Compute the mapping from the canonical domain to the non-affine quad on the surfaces of the cube
    
    % compute the (s,t) coordinates of the element
    s = (xi+1.0)*0.5;
    t = (eta+1.0)*0.5;
    
    % evaluate the boundary mappings
    gamma1xs = bsxfun(@plus, bsxfun(@times, s, (obj.nodesOfElements(2, 1, elements) - obj.nodesOfElements(1, 1, elements))), obj.nodesOfElements(1, 1, elements));
    gamma1ys = bsxfun(@plus, bsxfun(@times, s, (obj.nodesOfElements(2, 2, elements) - obj.nodesOfElements(1, 2, elements))), obj.nodesOfElements(1, 2, elements));
    gamma2xt = bsxfun(@plus, bsxfun(@times, t, (obj.nodesOfElements(3, 1, elements) - obj.nodesOfElements(2, 1, elements))), obj.nodesOfElements(2, 1, elements));
    gamma2yt = bsxfun(@plus, bsxfun(@times, t, (obj.nodesOfElements(3, 2, elements) - obj.nodesOfElements(2, 2, elements))), obj.nodesOfElements(2, 2, elements));
    gamma3xs = bsxfun(@plus, bsxfun(@times, s, (obj.nodesOfElements(3, 1, elements) - obj.nodesOfElements(4, 1, elements))), obj.nodesOfElements(4, 1, elements));
    gamma3ys = bsxfun(@plus, bsxfun(@times, s, (obj.nodesOfElements(3, 2, elements) - obj.nodesOfElements(4, 2, elements))), obj.nodesOfElements(4, 2, elements));
    gamma4xt = bsxfun(@plus, bsxfun(@times, t, (obj.nodesOfElements(4, 1, elements) - obj.nodesOfElements(1, 1, elements))), obj.nodesOfElements(1, 1, elements));
    gamma4yt = bsxfun(@plus, bsxfun(@times, t, (obj.nodesOfElements(4, 2, elements) - obj.nodesOfElements(1, 2, elements))), obj.nodesOfElements(1, 2, elements));

    
    gamma1x0 = obj.nodesOfElements(1, 1, elements);
    gamma1y0 = obj.nodesOfElements(1, 2, elements);
    gamma1x1 = obj.nodesOfElements(2, 1, elements);
    gamma1y1 = obj.nodesOfElements(2, 2, elements);
    gamma3x0 = obj.nodesOfElements(4, 1, elements);
    gamma3y0 = obj.nodesOfElements(4, 2, elements);
    gamma3x1 = obj.nodesOfElements(3, 1, elements);
    gamma3y1 = obj.nodesOfElements(3, 2, elements);


    %% Compute dX1dEta: the derivative of the mapping from the canonical domain to the surface of the cube
    dX1dEtaEvaluated = bsxfun(@times, (1.0 - s), gamma3x0 - gamma1x0) + bsxfun(@times, s, gamma3x1 - gamma1x1) - gamma1xs + gamma3xs -...
        bsxfun(@times, (1.0 - s), gamma3x0 - gamma1x0) - bsxfun(@times, s, gamma3x1 - gamma1x1);
    
    %% Compute dY1dEta: the derivative of the mapping from the canonical domain to the surface of the cube
    dY1dEtaEvaluated = bsxfun(@times, (1.0 - s), gamma3y0 - gamma1y0) + bsxfun(@times, s, gamma3y1 - gamma1y1) - gamma1ys + gamma3ys -...
        bsxfun(@times, (1.0 - s), gamma3y0 - gamma1y0) - bsxfun(@times, s, gamma3y1 - gamma1y1);

    
    %% Compute dXdEta
    dXdEtaEvaluated = dX1dEtaEvaluated;
    
end
function [x, y] = mapping(obj, elements, xi, eta)
% Mapping computes mapping between (xi,eta) and (x,y) for the
% Irregular Square.
%
% Computes the mapping between (xi, eta) and (x, y),
% the computational domain and the physical domain, respectively, for a
% given element number (element).
%
%
%   USAGE
%   -----
%       [x, y] = obj.mapping(element, xi, eta)
%
%       Computes the mapping between the nodes (xi,eta) of the
%       computational domain and the physical domain.
%
%       element can be for example: A, [A1,A2,A3] or [] (all elements)
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
%       x  :: the evaluation of the x component of the mapping function at
%             at the nodes (xi,eta) in the element (element).
%             (type: float64, size: array [N,M] (Q=1) or [N,M,Q] (Q>1))
%       y  :: the evaluation of the y component of the mapping function at
%             at the nodes (xi,eta) in the element (element).
%
%   Copyright 2018 Artur Palha

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
        
    % compute the (s,t) coordinates of the element
    s = (xi+1.0)*0.5;
    t = (eta+1.0)*0.5;
    
    % evaluate the mappings
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
    
    % compute the mapping between the canonical domain and the element
    % inside the subdomain in 2D
    x = bsxfun(@times, (1.0-s), gamma4xt) + bsxfun(@times, s, gamma2xt) + bsxfun(@times, (1.0-t), gamma1xs) + bsxfun(@times, t, gamma3xs) -...
        bsxfun(@times, (1.0-s), (bsxfun(@times, (1.0-t), gamma1x0) + bsxfun(@times, t, gamma3x0))) - bsxfun(@times, s, (bsxfun(@times, (1.0-t), gamma1x1) + bsxfun(@times, t, gamma3x1)));
    
    y = bsxfun(@times, (1.0-s), gamma4yt) + bsxfun(@times, s, gamma2yt) + bsxfun(@times, (1.0-t), gamma1ys) + bsxfun(@times, t,gamma3ys) -...
        bsxfun(@times, (1.0-s), (bsxfun(@times, (1.0-t), gamma1y0) + bsxfun(@times, t, gamma3y0))) - bsxfun(@times, s, (bsxfun(@times, (1.0-t), gamma1y1) + bsxfun(@times, t, gamma3y1))); 
      
end
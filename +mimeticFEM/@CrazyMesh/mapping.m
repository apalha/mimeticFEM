function [x,y] = mapping(obj,element,xi,eta)
% Mapping computes crazy mapping between (xi,eta) and (x,y)
%
% Computes the mapping between (xi,eta) and (x,y),
% the computational domain and the physical domain, respectively, for a
% given element number (element) and the number of elements in x and y
% directions (n).
%
%   USAGE
%   -----
%       [x, y] = obj.mapping(element,xi,eta)
%
%       Computes the mapping between the nodes (xi,eta) of the
%       computational domain and the physical domain.
%
%       element can be for example: A, [A1,A2,A3] or [] (all elements)
%
%   INPUTS
%   ------
%       element :: the element number
%                  If element is a scalar then x and y have dimensions
%                  (N,M), if element is a vector of dimension Q then
%                  x and y have dimensions (N,M,Q), where (N,M) are the
%                  dimensions of xi and eta. If element = [] then the
%                  mapping is computed for all elements
%                  (type: int32, size: vector)
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
%             (type: float64, size: array [N,M] (Q=1) or [N,M,Q] (Q>1))
%
%   Copyright 2009 Artur Palha

%   Revisions:  2009-11-25 First implementation. (apalha)
%               2017-01-17 Arbitrary number of elements can be computed at
%                          the same time, speeding up the computation.
    
    % the logical spacing of the elements
    deltaX = 2.0/obj.n(1);
    deltaY = 2.0/obj.n(2);
    
    % if element is empty [] then all the mapping is computed for all
    % elements
    if isempty(element)
        element = 1:obj.numElements;
    end
    
    % if element is an array it must be reshaped so that the x and y
    % coordinates are an array with three indices (i,j,k). i,j correspond
    % to the indices of the xi and eta arrays and the k index to the
    % elements.
    if ~isscalar(element)
        element = reshape(element,1,1,[]);
    end
    
    % x and y indices of the element
    iX = ceil(element/obj.n(2))-1;
    iY = mod(element-1,obj.n(2));
    
    xLeft = -1+deltaX*iX;
    yLeft = -1+deltaY*iY;
    
    x = (bsxfun(@plus,(xi+1)*0.5*deltaX,xLeft) +...
        obj.cc*sin(pi*bsxfun(@plus,(eta+1)*0.5*deltaY,yLeft)).*...
        sin(pi*bsxfun(@plus,(xi+1)*0.5*deltaX,xLeft)) + 1.0)*...
        (obj.xBounds(2)-obj.xBounds(1))*0.5 + obj.xBounds(1);
    y = (bsxfun(@plus,(eta+1)*0.5*deltaY,yLeft) +...
        obj.cc*sin(pi*bsxfun(@plus,(eta+1)*0.5*deltaY,yLeft)).*...
        sin(pi*bsxfun(@plus,(xi+1)*0.5*deltaX,xLeft)) + 1.0)*...
        (obj.yBounds(2)-obj.yBounds(1))*0.5 + obj.yBounds(1);
end
function dXdXiEvaluated = dXdXi(obj,element,xi,eta)
% dXdXi computes the dXdXi term of the metric for
% at the local positions (xi,eta) for the element. For the crazy grid.
%
%   USAGE
%   -----
%       dXdXiEvaluated = obj.dXdXi(element,xi,eta)
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
%       dXdXiEvaluated  :: the evaluation of the xi derivative of the mapping
%                          between the computational and physical spaces at 
%                          the points (xi,eta) of the element.
%                          (type: float64, size: array [N,M] (Q=1) or [N,M,Q] (Q>1))
%
%   Copyright 2009-2017 Artur Palha

%   Revisions:  2009-11-25 First implementation. (apalha)
%               2015-06-09 Updated documentation. (apalha) 
%               2017-01-21 Arbitrary number of elements can be computed at
%                          the same time, speeding up the computation. (apalha)
    
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
    
    dXdXiEvaluated = 0.5*deltaX*(obj.xBounds(2)-obj.xBounds(1))*0.5 +...
                        pi*deltaX*0.5*obj.cc*...
                        sin(bsxfun(@plus,pi*(eta+1)*0.5*deltaY,pi*yLeft)).*...
                        cos(bsxfun(@plus,pi*(xi+1)*0.5*deltaX,pi*xLeft))*...
                        (obj.xBounds(2)-obj.xBounds(1))*0.5;
end
function g22Evaluated = g22(obj,element,xi,eta)
% g22 computes the g22 term of the metric at the local positions
%   (xi,eta) for the element. For the SurfaceOfSphere.
%
%   USAGE
%   -----
%       g22Evaluated = obj.g22(element,xi,eta)
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
%       g22Evaluated  :: the evaluation of the g22 term of the
%                        metric at the points (xi,eta) of the
%                        element.
%                        (type: float64, size: array [N,M] (Q=1) or [N,M,Q] (Q>1))
%
%   Copyright 2009-2018 Artur Palha

%   Revisions:  2018-07-04 First implementation. (apalha)

    
    g22Evaluated = (obj.dXdXi(element,xi,eta).*...
                    obj.dXdXi(element,xi,eta) +...
                    obj.dYdXi(element,xi,eta).*...
                    obj.dYdXi(element,xi,eta) +...
                    obj.dZdXi(element,xi,eta).*...
                    obj.dZdXi(element,xi,eta))./...
                    (obj.g(element,xi,eta).^2); 

end
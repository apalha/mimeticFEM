function gEvaluated = g(obj, element, xi, eta)
% g computes the g term of the metric at the local positions
%   (xi,eta) for the element. For the SurfaceOfSphere.
%
%   USAGE
%   -----
%       gEvaluated = obj.g(element,xi,eta)
%
%   INPUTS
%   ------
%       element :: the element number
%                  (type: int32, size: single value)
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
%       gEvaluated  :: the evaluation of the g term of the
%                        metric at the points (xi,eta) of the
%                        element.
%                        (type: float64, size: array [N,M])
%
%   Copyright 2009-2018 Artur Palha

%   Revisions:  2018-07-04 (apalha) First implementation.
    
    element = obj.element;
    
 %   g0 = sqrt((obj.dXdXi(element,0,0).*obj.dYdEta(element,0,0)...
 %             -obj.dXdEta(element,0,0).*obj.dYdXi(element,0,0)).^2);
    g0 = 1.0;
    
    gEvaluated = sqrt((obj.dXdXi(element,xi,eta).*obj.dYdEta(element,xi,eta)...
                      -obj.dXdEta(element,xi,eta).*obj.dYdXi(element,xi,eta)).^2)/g0;
end

function gEvaluated = g(obj,element,xi,eta)
% g computes the g term of the metric at the local positions
%   (xi,eta) for the element. For the crazy grid.
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
%   Copyright 2009 Artur Palha

%   Revisions:  2009-11-25 (apalha) First implementation.
    
    gEvaluated = abs((obj.dXdXi(element,xi,eta).*...
                      obj.dYdEta(element,xi,eta) -...
                      obj.dXdEta(element,xi,eta).*...
                      obj.dYdXi(element,xi,eta)));

end

function g12Evaluated = g12(element,xi,eta,n)
% g12 computes the g12 term of the metric at the local positions
%   (xi,eta) for the element. For the pataki grid.
%
%   USAGE
%   -----
%       g12Evaluated = g12(element,xi,eta,n,epsilon,kappa,delta)
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
%       n :: the number of elements in x and y directions
%            (type: int32, size: array [1,2])
%
%
%   OUTPUTS
%   -------
%       g12Evaluated  :: the evaluation of the g12 term of the
%                        metric at the points (xi,eta) of the
%                        element.
%                        (type: float64, size: array [N,M])
%
%   Copyright 2015 Artur Palha

%   Revisions:  2015-07-30 (apalha) Updated documentation.
    
    g12Evaluated = (-mimeticFEM2.streak.dPhiYdEta(element,xi,eta,n).*...
                     mimeticFEM2.streak.dPhiYdXi(element,xi,eta,n)...
                    -mimeticFEM2.streak.dPhiXdEta(element,xi,eta,n).*...
                     mimeticFEM2.streak.dPhiXdXi(element,xi,eta,n))./...
                     (mimeticFEM2.streak.g(element,xi,eta,n).^2);                      
end
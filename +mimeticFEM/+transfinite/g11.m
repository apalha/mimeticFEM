function g11Evaluated = g11(element,xi,eta,n,gamma,dgamma)
% g11 computes the g11 term of the metric at the local positions
%   (xi,eta) for the element. For a transfinite grid.
%
%   USAGE
%   -----
%       g11Evaluated = g11(element,xi,eta,n,gamma,dgamma)
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
%       gamma :: the matlab functions for the four boundaries that define the
%                domain, as given in [1], (36). This is a cell containing
%                four functions for gamma_{1}, gamma_{2}, gamma_{3} and
%                gamma_{4}. These are functions of s with s \in [0,1]. If
%                we define four vertices in the boundary of our domain and
%                we number them in counter-clockwise fashion as v_{1}, v_{2},
%                v_{3}, v_{4}, then: 
%                   - gamma_{1} defines a path from v_{1} to v_{2}
%                   - gamma_{2} defines a path from v_{2} to v_{3}
%                   - gamma_{3} defines a path from v_{4} to v_{3}
%                   - gamma_{4} defines a path from v_{1} to v_{4}
%                (type: cell of matlab functins, size: 4)
%       dgamma :: the matlab functions for the derivatives of the
%                 boundaries
%                 that define the domain, assumed it is a transfinite 
%                 mapping as given in [1], (36). This is a cell
%                 containing four functions for (dgamma_{1}xds, dgamma_{1}yds),
%                 (dgamma_{2}yds, dgamma_{2}yds), (dgamma_{3}xds, dgamma_{3}yds),
%                 (dgamma_{4}xds and dgamma_{4}yds). These are functions of s
%                 with s \in [0,1]. If
%                 we define four vertices in the boundary of our domain and
%                 we number them in counter-clockwise fashion as v_{1}, v_{2},
%                 v_{3}, v_{4}, then: 
%                    - gamma_{1} defines a path from v_{1} to v_{2}
%                    - gamma_{2} defines a path from v_{2} to v_{3}
%                    - gamma_{3} defines a path from v_{4} to v_{3}
%                    - gamma_{4} defines a path from v_{1} to v_{4}
%
%                 the cell contains then:
%                    - dgamma_{1} is dgamma_{1}xds,dgamma_{1}yds
%                    - dgamma_{2} is dgamma_{2}xds,dgamma_{2}yds
%                    - dgamma_{3} is dgamma_{3}xds,dgamma_{3}yds
%                    - dgamma_{4} is dgamma_{4}xds,dgamma_{4}yds
%                 (type: cell of matlab functins, size: 4)
%
%
%   OUTPUTS
%   -------
%       g11Evaluated  :: the evaluation of the g11 term of the
%                        metric at the points (xi,eta) of the
%                        element.
%                        (type: float64, size: array [N,M])
%
%   Copyright 2015 Artur Palha

%   Revisions:  2015-08-05 (apalha) Updated documentation.

    g11Evaluated = (mimeticFEM2.transfinite.dPhiXdEta(element,xi,eta,n,gamma,dgamma).*...
                    mimeticFEM2.transfinite.dPhiXdEta(element,xi,eta,n,gamma,dgamma) +...
                    mimeticFEM2.transfinite.dPhiYdEta(element,xi,eta,n,gamma,dgamma).*...
                    mimeticFEM2.transfinite.dPhiYdEta(element,xi,eta,n,gamma,dgamma))./...
                    (mimeticFEM2.transfinite.g(element,xi,eta,n,gamma,dgamma).^2);
end
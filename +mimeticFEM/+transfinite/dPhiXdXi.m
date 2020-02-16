function dPhiXdXiEvaluated = dPhiXdXi(element,xi,eta,n,gamma,dgamma)
% dPhiXdXi computes the dPhiXdXi term of the metric for
% at the local positions (xi,eta) for the element. This is done for a
% transfinite mapping.
%
% The transfinite mapping is between the unit square [0,1]x[0,1],
% subdivided into n(1) x n(2) elements of equal sizes, and the domain
% bounded by the four line segments \Gamma_{1}, \Gamma_{2}, \Gamma_{3} and
% \Gamma_{4}.
%
%   The global mapping is generated according to the known transfinite mapping [1]: 
%       (s,t) \in [0,1]x[0,1] --> (x,y)
%
%   (x,y) = T(s,t) = (1-s)\Gamma_{4}(t) + s\Gamma_{2}(t) + (1-t)\Gamma_{1}(s) +
%                    + t \Gamma_{3}(s) - (1-s)((1-t)\Gamma_{1}(0) + t\Gamma_{3}(0)) -
%                    - s((1-t)\Gamma_{1}(1) + t\Gamma_{3}(1))
%
%   The mesh corresponds to a subdivision of (s,t) \in [0,1]x[0,1]
%   into n(1) x n(2) elements of equal sizes in x- and y-direction.
%
%   For the element k with 1 <= k <= n(1)*n(2) (the 2D numbering (i,j) goes such
%   that k = (i-1)n(2) + j) the mapping is:
%
%       (\xi,\eta) \in [-1,1]x[-1,1] --> (r,z)
%
%   First we introduce the quantities:
%
%   \Delta_{s} = 1.0/n(1);
%   \Delta_{t} = 1.0/n(2);
%
%   s_{L} = \Delta_{s}*i;
%   t_{L} = \Delta_{t}*j;
%
%   s(\xi,k) = ((\xi+1)*0.5*\Delta_{s}+s_{L}(k))
%   t(\eta,k) = ((\xi+1)*0.5*\Delta_{t}+t_{L}(k))
%
%   Therefore the local mapping from each element k into the physical
%   domain is given by:
%
%   (x,y) = T(s(\xi,k),t(\eta,k))
%
%   USAGE
%   -----
%       dPhiXdXiEvaluated = dPhiXdXi(element,xi,eta,n,gamma,dgamma)
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
%       dPhiXdXiEvaluated  :: the evaluation of the eta derivative of the mapping
%                              between the computational and physical spaces at 
%                              the points (xi,eta) of the element.
%                              (type: float64, size: array [N,M])
%
%
%   REFERENCES
%   ----------
%       [1] Kopriva, D., Implementing spectral methods for partial differential
%           equations: algorithms for scientists and engineers, Springer, 2009.
%
%   Copyright 2015 Artur Palha

%   Revisions:  2015-08-05 (apalha) First implementation.

    % compute the element spacing in the (s,t) space
    deltaS = 1.0/n(1);
    deltaT = 1.0/n(2);
    
    % x and y indices of the element
    iS = ceil(element/n(2))-1;
    iT = mod(element-1,n(2));
    
    % compute the lower left vertex coordinates of the element
    sLeft = deltaS*iS;
    tLeft = deltaT*iT;
    
    % compute the (s,t) coordinates of the element
    s = (xi+1.0)*deltaS*0.5 + sLeft;
    t = (eta+1.0)*deltaT*0.5 + tLeft;
    
    % evaluate the mappings
    %[gamma1xs,gamma1ys] = gamma{1}(s);
    [gamma2xt,~] = gamma{2}(t);
    %[gamma3xs,gamma3ys] = gamma{3}(s);
    [gamma4xt,~] = gamma{4}(t);
    
    [gamma1x0,~] = gamma{1}(0.0);
    [gamma1x1,~] = gamma{1}(1.0);
    [gamma3x0,~] = gamma{3}(0.0);
    [gamma3x1,~] = gamma{3}(1.0);
    
    % evaluate the derivatives of the mappings
    [dgamma1xds,~] = dgamma{1}(s);
    %[dgamma2xdt,dgamma2ydt] = dgamma{2}(t);
    [dgamma3xds,~] = dgamma{3}(s);
    %[dgamma4xdt,dgamma4ydt] = dgamma{4}(t);


    dPhiXdXiEvaluated = -gamma4xt + gamma2xt + (1.0-t).*dgamma1xds + t.*dgamma3xds +...
                         ((1.0-t).*gamma1x0 + t.*gamma3x0) - ((1.0-t).*gamma1x1 + t.*gamma3x1);
    
    dPhiXdXiEvaluated = dPhiXdXiEvaluated*deltaT*0.5;
end







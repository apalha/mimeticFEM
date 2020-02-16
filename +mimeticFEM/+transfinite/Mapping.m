function [x,y] = Mapping(element,xi,eta,n,gamma)
% Mapping computes the transfinite mapping between the unit square [0,1]x[0,1]
% subdivided into n(1) x n(2) elements of equal sizes and the domain
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
%   
%   USAGE
%   -----
%       [x, y] = Mapping(element,xi,eta,n,gamma)
%
%       Computes the mapping between the nodes (xi,eta) of the
%       computational domain and the physical domain.
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
%
%   OUTPUTS
%   -------
%       x  :: the evaluation of the r component of the mapping function at
%             at the nodes (xi,eta) in the element (element).
%             (type: float64, size: array [N,M])
%       y  :: the evaluation of the z component of the mapping function at
%             at the nodes (xi,eta) in the element (element).
%             (type: float64, size: array [N,M])
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
    [gamma1xs,gamma1ys] = gamma{1}(s);
    [gamma2xt,gamma2yt] = gamma{2}(t);
    [gamma3xs,gamma3ys] = gamma{3}(s);
    [gamma4xt,gamma4yt] = gamma{4}(t);
    
    [gamma1x0,gamma1y0] = gamma{1}(0.0);
    [gamma1x1,gamma1y1] = gamma{1}(1.0);
    [gamma3x0,gamma3y0] = gamma{3}(0.0);
    [gamma3x1,gamma3y1] = gamma{3}(1.0);
    
    % compute the mapping as given in [1]
    x = (1.0-s).*gamma4xt + s.*gamma2xt + (1.0-t).*gamma1xs + t.*gamma3xs -...
        (1.0-s).*((1.0-t).*gamma1x0 + t.*gamma3x0) - s.*((1.0-t).*gamma1x1 + t.*gamma3x1);
    
    y = (1.0-s).*gamma4yt + s.*gamma2yt + (1.0-t).*gamma1ys + t.*gamma3ys -...
        (1.0-s).*((1.0-t).*gamma1y0 + t.*gamma3y0) - s.*((1.0-t).*gamma1y1 + t.*gamma3y1);    
end







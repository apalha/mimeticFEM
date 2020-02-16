function [x,y] = gamma2(s,region)
% gamma2 computes the parametric curve of the right part of the streak 
%   domain as defined in [1].
%
%   
%   USAGE
%   -----
%       [x, y] = gamma2(s,region)
%
%       Computes the parameterization between the points s of the
%       computational domain and the physical domain for the three
%       different regions.
%
%   INPUTS
%   ------
%       s :: the local coordinates of the parameter where to
%             evaluate the curve. 
%       region :: defines which of the three regions of the streak is to be
%                 considered.
%                 1: top
%                 2: middle
%                 3: bottom
%                (type: int, size: single value)
%
%
%   OUTPUTS
%   -------
%       x  :: the evaluation of the r component of the mapping function at
%             at the points s.
%             (type: double)
%       y  :: the evaluation of the z component of the mapping function at
%             at the points s.
%             (type: double)
%
%   REFERENCES
%   ----------
%       [1] Put the reference from where you took the mapping.
%
%   Copyright 2015 Artur Palha

%   Revisions:  2015-08-05 (apalha) First implementation.
    
    % left points where the regions end or start
    % points are numbered from lower to top
    xL1 = 0.0;
    xL2 = 0.0;
    xL3 = 0.0;
    xL4 = 0.0;
    yL1 = 0.0;
    yL2 = 0.6954451150103;
    yL3 = 0.7958260743101;
    yL4 = 1.0;
    
    % right points where the regions end or start
    % points are numbered from lower to top
    xR1 = 1.0;
    xR2 = 1.0;
    xR3 = 1.0;
    xR4 = 1.0;
    yR1 = 0.0;
    yR2 = 0.2324555320337;
    yR3 = 0.3937253933194;
    yR4 = 1.0;
    
    if region == 1
        x = s*(xR4-xR3) + xR3;
        y = s*(yR4-yR3) + yR3;
    elseif region == 2
        x = s*(xR3-xR2) + xR2;
        y = s*(yR3-yR2) + yR2;
    elseif region == 3
        x = s*(xR2-xR1) + xR1;
        y = s*(yR2-yR1) + yR1;
    end
    
end

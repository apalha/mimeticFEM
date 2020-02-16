function [x,y] = gamma1(s,region)
% gamma1 computes the parametric curve of the bottom part of the streak 
%   domain as defined in [1].
%
%   
%   USAGE
%   -----
%       [x, y] = gamma1(s,region)
%
%       Computes the parameterization between the points s of the
%       computational domain and the physical domain for the three
%       different regions.
%
%   INPUTS
%   ------
%       s :: the local coordinates of the parameter where to
%            evaluate the curve. s \in [0,1].
%            (type: double, size: any array)
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
    
    
    if region == 1
        alpha1 = 1.6542264134055116;
        alpha2 = 0.72273424781341566;
        R = 1.2;
        x0 = 0.1;
        y0 = -0.4;
        
        x = R*cos(s*(alpha2-alpha1) + alpha1) + x0;
        y = R*sin(s*(alpha2-alpha1) + alpha1) + y0;
    elseif region == 2
        alpha1 = 1.6618311048323118;
        alpha2 = 0.61255473833933904;
        R = 1.1;
        x0 = 0.1;
        y0 = -0.4;
        
        x = R*cos(s*(alpha2-alpha1) + alpha1) + x0;
        y = R*sin(s*(alpha2-alpha1) + alpha1) + y0;
    elseif region == 3
        x = s;
        y = zeros(size(s));
    end
    
end

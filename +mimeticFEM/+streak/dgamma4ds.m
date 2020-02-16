function [dgamma4xds_evaluated,dgamma4yds_evaluated] = dgamma4ds(s,region)
% dgamma4ds computes the derivative of the left boundary
%    with respect to s, for the Pataki grid given in [1].
%
%
%   
%   USAGE
%   -----
%       [dgamma4xds_evaluated,dgamma4xds_evaluated] = dgamma4ds(s,region)
%
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
%       dgamma4xds_evaluated  :: the evaluation of the x component of the derivative
%                                of the parametric curve of the left
%                                boundary.
%                                (type: double)
%       dgamma4yds_evaluated  :: the evaluation of the y component of the derivative
%                                of the parametric curve of the left
%                                boundary.
%                                (type: double)
%
%   REFERENCES
%   ----------
%       [1] 
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
%         x = s*(xL4-xL3) + xL3;
%         y = s*(yL4-yL3) + yL3;
        dgamma4xds_evaluated = (xL4-xL3);
        dgamma4yds_evaluated = (yL4-yL3);
    elseif region == 2
%         x = s*(xL3-xL2) + xL2;
%         y = s*(yL3-yL2) + yL2;
        dgamma4xds_evaluated = (xL3-xL2);
        dgamma4yds_evaluated = (yL3-yL2);
    elseif region == 3
%         x = s*(xL2-xL1) + xL1;
%         y = s*(yL2-yL1) + yL1;
        dgamma4xds_evaluated = (xL2-xL1);
        dgamma4yds_evaluated = (yL2-yL1);
    end

end
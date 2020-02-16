function [dgamma1xds_evaluated,dgamma1yds_evaluated] = dgamma1ds(s,region)
% dgamma1ds computes the derivative of the bottom boundary
%    with respect to s, for the Pataki grid given in [1].
%
%
%   
%   USAGE
%   -----
%       [dgamma1xds_evaluated,dgamma1yds_evaluated] = dgamma1ds(s,region)
%
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
%   OUTPUTS
%   -------
%       dgamma1xds_evaluated  :: the evaluation of the x component of the derivative
%                                of the parametric curve of the bottom
%                                boundary.
%                                (type: double)
%       dgamma1yds_evaluated  :: the evaluation of the y component of the derivative
%                                of the parametric curve of the bottom
%                                boundary.
%                                (type: double)
%
%   REFERENCES
%   ----------
%       [1] 
%
%   Copyright 2015 Artur Palha

%   Revisions:  2015-08-05 (apalha) First implementation.
    
    if region == 1
        alpha1 = 1.6542264134055116;
        alpha2 = 0.72273424781341566;
        R = 1.2;
        x0 = 0.1;
        y0 = -0.4;

%         x = R*cos(s*(alpha2-alpha1) + alpha1) + x0;
%         y = R*sin(s*(alpha2-alpha1) + alpha1) + y0;  
        dgamma1xds_evaluated =-(alpha2-alpha1)*R*sin(s*(alpha2-alpha1) + alpha1);
        dgamma1yds_evaluated = (alpha2-alpha1)*R*cos(s*(alpha2-alpha1) + alpha1);
    elseif region == 2
        alpha1 = 1.6618311048323118;
        alpha2 = 0.61255473833933904;
        R = 1.1;
        x0 = 0.1;
        y0 = -0.4;
        
%         x = R*cos(s*(alpha2-alpha1) + alpha1) + x0;
%         y = R*sin(s*(alpha2-alpha1) + alpha1) + y0;
        dgamma1xds_evaluated =-(alpha2-alpha1)*R*sin(s*(alpha2-alpha1) + alpha1);
        dgamma1yds_evaluated = (alpha2-alpha1)*R*cos(s*(alpha2-alpha1) + alpha1);
    elseif region == 3
%         x = s;
%         y = zeros(size(s));
        
        dgamma1xds_evaluated = 1;
        dgamma1yds_evaluated = 0;
    end

end
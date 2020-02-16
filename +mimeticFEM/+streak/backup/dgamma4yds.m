function dgamma4yds_evaluated = dgamma4yds(s,epsilon,kappa,delta)
% dgamma4yds computes the derivative of the y component of the left boundary
%    with respect to s, for the Pataki grid given in [1].
%
%
%   
%   USAGE
%   -----
%       dgamma4yds_evaluated = dgamma4yds(s,epsilon,kappa,delta)
%
%
%   INPUTS
%   ------
%       element :: the element number
%                  (type: int32, size: single value)
%       s :: the local coordinates of the parameter where to
%             evaluate the parameterization.
%             (type: float64, size: array [N,M])
%       epsilon :: as given in [1], related to the size of the horizontal
%                  axis of the plasma shape.
%                  (type: float64, size: single value)
%       kappa :: as given in [1], related to the size of the vertical
%                axis of the plasma shape (kappa*epsilon).
%                (type: float64, size: single value)  
%       delta :: as given in [1], related to the misalignment between the
%                top of the plasma shape and the  axis.
%                axis of the plasma shape.
%                (type: float64, size: single value)
%
%
%
%   OUTPUTS
%   -------
%       dgamma4yds_evaluated  :: the evaluation of the y component of the derivative
%                                of the parametric curve of the left
%                                boundary.
%                                (type: float64, size: array [N,M])
%
%   REFERENCES
%   ----------
%       [1] Pataki, A., Cerfon, A. J., Freidberg, J. P., Greengard, L., & O?Neil, M. (2013).
%           A fast, high-order solver for the Grad?Shafranov equation. 
%           Journal of Computational Physics, 243, 28?45. doi:10.1016/j.jcp.2013.02.045
%
%   Copyright 2015 Artur Palha

%   Revisions:  2015-08-05 (apalha) First implementation.
    
    alpha = asin(delta);
    s = -(2.0*pi/4.0)*s + 5.0*pi/4.0;
    
    dgamma4yds_evaluated = epsilon*kappa*cos(s)*2.0*pi/4.0;


end
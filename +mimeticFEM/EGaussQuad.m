function [x,w] = EGaussQuad(p)
%EGaussQuad Returns the p+1 extened Gauss nodes and weights of extended Gauss
% quadrature.
%
%   Extended Gauss nodes of order p are simply the (p-1) Gauss nodes of
%   order (p-2) plus the two end nodes: -1 and 1. The weights are the (p-1)
%   Gauss quadrature weights plus 0 at the ends.
%
%   USAGE
%   -----
%       [x w] = EGaussQuad(p) 
%
%       Gives the p+1 nodes and weights of the Extended Gauss quadrature of
%       order p. 
%   
%   INPUTS
%   ------
%       p :: quadrature order
%            (type: int32, size: single value)
%
%   OUTPUTS
%   -------
%       x :: nodes of Extended Gauss quadrature of order p. x \in [-1,1].
%            (type: float64, size: [1, p+1])
%       w :: weights of Extended Gauss quadrature of order p.
%            (type: float64, size: [1, p+1])
%
%   Copyright 2009 Artur Palha

%   Revisions:  2009-11-26 (apalha) First implementation.
%               2014-12-03 (apalha) Converted EGaussNodes to EGaussQuad,
%                                   for that added the quadrature weights
%                                   to this function.

    % allocate space for the nodes and weights
    x = zeros(p+1,1);
    w = zeros(p+1,1);
    
    if p > 1,
        [x(2:end-1), w(2:end-1)] = mimeticFEM.GaussQuad(p-2);
    end
    
    % add the end nodes
    x(1) = -1.0;
    x(end) = 1.0;
    
    % no need to add the quadrature weights associated to the end nodes -1
    % and 1, since they are 0 and the weights vector, w, is initialized with 0's.
end
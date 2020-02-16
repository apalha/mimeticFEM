function [x w] = LobattoQuad(p)
%LobattoQuad Returns the p+1 Lobatto points and weights of Gauss-Lobatto
% quadrature.
%
%
% For the computation of the nodes it uses a Newton method
% up to machine precision. As initial guess it uses the Chebychev roots.
% for the roots.
%
%   USAGE
%   -----
%       [x w] = LobattoQuad(p) 
%
%       Gives the p+1 nodes and weights of the Gauss-Lobatto quadrature of
%       order p. 
%   
%   INPUTS
%   ------
%       p :: quadrature order
%            (type: int32, size: single value)
%
%   OUTPUTS
%   -------
%       x :: nodes of Gauss-Lobatto quadrature of order p. x \in [-1,1].
%            (type: float64, size: [1, p+1])
%       w :: weights of Gauss-Lobatto quadrature of order p.
%            (type: float64, size: [1, p+1])
%
%   Copyright 2009 Artur Palha

%   Revisions:  2009-11-25 (apalha) First implementation.

    n = p+1;
    x = cos(pi*(0:n-1)/p)';
    P = zeros(n,n);
    xold = 2.0;
    while max(abs(x-xold)) > eps,
        xold = x;
        P(:,1) = 1.0;
        P(:,2) = x;
        for k=3:n,
            P(:,k) = ((2*(k-1)-1).*x.*P(:,k-1) - (k-2)*P(:,k-2))/(k-1);
        end
        x = xold - (x.*P(:,n) - P(:,n-1))./(n*P(:,n));
    end
    x = fliplr(x')';
    w = 2.0./(p*n*(P(:,n)).^2);
    
end
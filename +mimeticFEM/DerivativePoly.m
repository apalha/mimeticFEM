function result = DerivativePoly(x, p, polyType)
%DerivativePoly returns the derivative of the polynomial interpolant basis
%   functions at the points x.
%
%   The derivative of the Lagrange interpolating basis functions
%   (l_{n}^{p}(x)) are given by:
%   
%   \frac{dl_{n}(x)}{dx} = \sum_{i=1}_{i\neq n}^{p+1}\prod_{j\neq n}_{j\neq i}
%                          \frac{1}{x_{n}-x_{i}}\frac{x-x_{j}}{x_{n}-x_{j}}
%
%   This computation is very slow computed like this. Hence, a clever way
%   is followed.
%   
%       1- compute the derivatives of the basis functions at the nodal
%          points (using DERIVATIVEPOLYNODES)
%       2- compute the basis functions at the x points)
%       3- update the values of the derivatives at the x points
%
%   In this way, instead of having an algorithm of O(p^{2}) one has an
%   algorithm of O(p).
%
%   It returns a 2-dimensional matrix with the values of the derivative of
%   the polynomials of order p of type polyType.
%
%   If x is a vector of length N it returns a 2d matrix whose rows are the 
%   values of the evaluated polynomial, P(x), in x:
%                 -                                                  -
%       result = | dP_{1}(x(1))   dP_{1}(x(2))   ...   dP_{1}(x(N))   |
%                | dP_{2}(x(1))   dP_{2}(x(2))   ...   dP_{2}(x(N))   |
%                |                      ...                           |
%                | dP_{p+1}(x(1)) dP_{p+1}(x(2)) ...   dP_{p+1}(x(N)) |
%                 -                                                  -  
%
%   USAGE
%   -----
%       result = DerivativePoly(x, p, polyType)
%
%           Computes the derivatives of the Lagrange interpolants associated
%           to the polyType at the nodes x.
%
%   INPUTS
%   ------
%       x :: Locations where to evaluate the derivatives.
%            x \in [-1,1].
%            (type: float64, size: [N,1], [1,N])
%       p :: The order of the polynomials.
%            (type: int32, size: single value)
%       polyType :: the type of polynomial to compute the derivative.
%                   Valid values: 'Lobatto', 'Gauss', 'EGauss'.
%                   (type: string, size: single string)
%
%   OUTPUTS
%   -------
%       result :: The p+1 derivatives polynomials evaluated at the x points.
%                 (type: float64, size: [p+1, N])
%
%   Copyright 2009 Artur Palha

%   Revisions:  2009-11-25 (apalha) First implementation.
%               2014-12-03 (apalha) Removed pre-allocation of result.
%                                   Removed the final reshaping since x is
%                                   now always a vector. If x is an array,
%                                   one should use x(:) as input.

    % check if polyType is a valid one
    if ~mimeticFEM.TestPolyType(polyType)
        fprintf(':: %s :: is not a valid type of polynomial', polyType);
        return
    end
    
    % % shape of x
    % sizeOfx = size(x);
    % lengthOfx = numel(x);
    
    % transform x in a vector
    x = x(:);
    
    % % preallocate memory space for result
    % result = zeros(p+1, lengthOfx);
    
    % compute the derivatives at the nodal points
    nodalDerivatives = mimeticFEM.DerivativePolyNodes(p, polyType, [true false]);
    
    % update the derivatives to the x points
    
    % compute the basis polynomials in the x points
    polyEval = (eval(sprintf('mimeticFEM.%sPoly(%s, %s)', polyType, 'x', 'p')))';
    
    % compute the derivatives in the x points
    result = (polyEval * nodalDerivatives)';
    
    % if min(sizeOfx) ~= 1
    %    result = reshape(result, [p+1 sizeOfx]);
    % end
end
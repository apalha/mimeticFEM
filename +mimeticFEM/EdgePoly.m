function result = EdgePoly(x, p, polyType)
%EdgePoly Returns the edge basis polynomials associated to the polyType, 
%   evaluated at the x points.
%
%   The edge basis functions are given by, see [1]:
%
%       edge_{i}(x) = -\sum_{j=1}^{i} dh_{j}(x), i=1,...,p
%
%   dh_{j}(x) are computed using DerivativePoly.
%
%   It returns a 2-dimensional matrix with the values of the edge basis
%   polynomials evaluated at x.
%
%   If x is a vector of length N it returns a 2d matrix whose rows are the 
%   values of the evaluated polynomial, P(x), in x:
%                 -                                               -
%       result = | P_{1}(x(1))   P_{1}(x(2))   ...   P_{1}(x(N))   |
%                | P_{2}(x(1))   P_{2}(x(2))   ...   P_{2}(x(N))   |
%                |                   ...                           |
%                | P_{p+1}(x(1)) P_{p+1}(x(2)) ...   P_{p+1}(x(N)) |
%                 -                                               -  
%
%   If x=[] then it computes the edge basis polynomial basis at polyType
%   nodes.
%
%   USAGE
%   -----
%       result = EdgePoly(x, p, polyType)
%
%           Computes the edge basis polynomials of order (p-1) at the x nodes 
%       
%       result = EdgeFunction([], p, polyType)
%
%           Computes the edges basis polynomials of order (p-1) at the x
%           nodes of the polyType. For example, if polytype = 'Gauss', the
%           edge basis polynomials will be evaluated at the Gauss nodes.
%
%   INPUTS
%   ------
%       x :: Locations where to evaluate the edge basis polynomials.
%            x \in [-1,1].
%            (type: float64, size: [N,1], [1,N] or [])
%       p :: The order of the edge basis polynomials.
%            (type: int32, size: single value)
%       polyType :: the type of polynomial used to construct the edge basis
%                   polynomials.
%                   Valid values: 'Lobatto', 'Gauss', 'EGauss'.
%                   (type: string, size: single string)
%
%   OUTPUTS
%   -------
%       result :: The p polynomials evaluated at the x points.
%                 (type: float64, size: [p, N] or [p,p+1] if x=[])
%   
%
%   [1] Gerritsma, M.: Edge functions for spectral element methods. 
%       Submitted to the proceedings of ICOSAHOM 2009
%
%   Copyright 2009 Artur Palha

%   Revisions:  2009-12-04 (apalha) First implementation.
%               2014-12-03 (apalha) Converted to the same format as
%                                   LobattoPoly, GaussPoly and EGaussPoly.
    
    % check if polyType is a valid one
    if ~mimeticFEM.TestPolyType(polyType)
        disp(sprintf(':: %s :: is not a valid type of polynomial', polyType));
        return
    end
    
    if isempty(x)
        % if x is empty then compute the derivatives at the nodes of
        % polyType
        derivatives = mimeticFEM.DerivativePolyNodes(p,polyType);
    else
        % if x is not empty compute the derivatives at x
        derivatives = mimeticFEM.DerivativePoly(x, p, polyType);
    end
    
    % compute the edge basis functions
    result = cumsum(derivatives);   
    result = -result(1:(end-1),:);
end
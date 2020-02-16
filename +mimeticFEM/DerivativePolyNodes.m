function result = DerivativePolyNodes(p, polyType, options)
%DerivativePoly returns the derivative of the polynomial interpolant basis
%   functions at the nodal points.
%
%   The derivative of the Lagrange interpolating basis functions
%   (l_{n}^{p}(x)) are given by:
%   
%   \frac{dl_{n}(x)}{dx} = \sum_{i=1}_{i\neq n}^{p+1}\prod_{j\neq n}_{j\neq i}
%                          \frac{1}{x_{n}-x_{i}}\frac{x-x_{j}}{x_{n}-x_{j}}
%
%
%   For computation at the nodes a more efficient and accurate formula can
%   be used, see [1]:
%       
%             | \frac{c_{k}}{c_{j}}\frac{1}{x_{k}-x_{j}},          k \neq j
%             |
%   d_{kj} = <
%             | \sum_{l=1,l\neq k}^{p+1}\frac{1}{x_{k}-x_{l}},     k = j
%             |
%   
%   with
%
%       c_{k} = \prod_{l=1,l\neq k}^{p+1} (x_{k}-x_{l})   
%
%   It returns a 2-dimensional matrix with the values of the derivative of
%   the polynomials of order p of type polyType.
%                 -                                                    -
%       result = | dP_{1}(x(1))   dP_{1}(x(2))   ...   dP_{1}(x(p+1))   |
%                | dP_{2}(x(1))   dP_{2}(x(2))   ...   dP_{2}(x(p+1))   |
%                |                      ...                             |
%                | dP_{p+1}(x(1)) dP_{p+1}(x(2)) ...   dP_{p+1}(x(p+1)) |
%                 -                                                    -  
%
%   USAGE
%   -----
%       result = DerivativePolyNodes(p, polyType, [options])
%
%           Computes the derivatives of the Lagrange interpolants associated
%           to the polyType at the polyType nodes.
%
%   INPUTS
%   ------
%       p :: The order of the polynomials.
%            (type: int32, size: single value)
%       polyType :: The type of polynomial to compute the derivative.
%                   Valid values: 'Lobatto', 'Gauss', 'EGauss'.
%                   (type: string, size: single string)
%       options :: Options that define which formula to use for the
%                  computation of the derivative and if the results should
%                  be transposed or not.
%                  options(1): false --> use formula (6) from [1]
%                              true  --> use formula (9) from [1]
%                  
%                  options(2): true  --> transpose results
%                              false --> do not transpose results
%                  (type: bool, size: [1,2])
%
%   OUTPUTS
%   -------
%       result :: The p+1 derivatives polynomials evaluated at the nodal points.
%                 (type: float64, size: [p+1, p+1])
%
%
%   [1] Costa, B., Don, W. S.: On the computation of high order
%       pseudospectral derivatives, Applied Numerical Mathematics, vol.33
%       (1-4), pp. 151-159
%
%
%   Copyright 2009 Artur Palha

%   Revisions:  2009-11-25 (apalha) First implementation.
%               2014-12-03 (apalha) Removed pre-allocation of result.
%                                   Replaced repmats by bsxfun for smaller
%                                   memory footprint.
    
    if nargin == 2
        options = [true true];
    end
    
    % check if polyType is a valid one
    if ~mimeticFEM.TestPolyType(polyType)
        fprintf(':: %s :: is not a valid type of polynomial', polyType);
        return
    end
    
    % % preallocate memory space for result
    % result = zeros(p+1, p+1);
    
    % compute the nodes of the type of polynomial
    if strcmp(polyType,'EGauss')
        if p >=2
            roots = [-1; mimeticFEM.GaussQuad(p-2); 1];
        else
            roots = [-1 1]';
        end
    else
        roots = eval(sprintf('mimeticFEM.%sQuad(%d)', polyType, p));
    end
    % compute all (ri - rj) = xi_xj(i,j)
    % xi_xj = repmat(roots, [1 p+1]) - repmat(roots', [p+1 1]);
    xi_xj = bsxfun(@minus,roots,roots');
    
    % transform all the diagonals to 1
    % xi_xj(linspace(1,(p+1)^2,(p+1))) = 1.0;
    xi_xj(1:p+2:(p+1)^2) = 1.0;
        
    % compute (ci's)
    c = prod(xi_xj, 2);
    
    % compute ck/cj = ck_cj(k,j) matrix
    % ci_cj = repmat(c, [1 p+1]) ./ repmat(c', [p+1 1]);
    ci_cj = bsxfun(@rdivide,c,c');
    
    % compute the off diagonal results
    result = ci_cj ./ xi_xj;
    
    % compute the diagonal components
    
    if ~options(1)
        % using formula (6) of [1]
        xi_xj = 1.0 ./ xi_xj;
    
        % transform the diagonal elements to 0.0
        xi_xj(linspace(1,(p+1)^2,(p+1))) = 0.0;
    
        % update the diagonal result with the values of formula (6)
        % result(linspace(1,(p+1)^2,(p+1))) = sum(xi_xj, 2);
        result(1:p+2:(p+1)^2) = sum(xi_xj, 2);
    
    else
        % using formula (9) of [1]
    
        % put all diagonal values equal to 0.0
        % result(linspace(1,(p+1)^2,(p+1))) = 0.0;
        result(1:p+2:(p+1)^2) = 0.0;
    
        % compute the diagonal values
        % result(linspace(1,(p+1)^2,(p+1))) = -sum(result,2);
        result(1:p+2:(p+1)^2) = -sum(result,2);
    end
    
    if options(2)
        % transpose the result
        result = result';
    end
end
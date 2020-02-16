function result = EGaussPoly(x,p)
%EGaussPoly Returns the p+1 Extended Gauss Lagrange interpolant basis
%  functions evaluated at x.
%
%   USAGE
%   -----
%       result = EGaussPoly(x,p)
%
%       Gives the p+1 Extended Gauss basis interpolants evaluated at x.
%
%   INPUTS
%   ------
%       x :: Locations where to evaluate the polynomial basis functions.
%            x \in [-1,1].
%            (type: float64, size: [N,1] or [1,N])
%       p :: The order of the polynomial basis functions.
%            (type: int32, size: single value)
%
%   OUTPUTS
%   -------
%       result :: The (p+1) polynomials evaluated at the x points.
%                 (type: float64, size: [p+1, N])
%
%  Copyright 2009-11-25 Artur Palha


%   Revisions:  2009-11-25 (apalha) First implementation.
    
    sizeOfx = size(x);
    if sizeOfx(1) == 1,
        sizeOfx = sizeOfx(2);
    elseif sizeOfx(2) == 1,
        sizeOfx = sizeOfx(1);
        x = x';
    end
    
    if length(x) == 0
        result = speye(p+1);
        return
    end
    
    % allocate memory space for the result
%    result = zeros([p+1 sizeOfx]);
    result = ones([p+1 sizeOfx]);
    % compute EGauss roots
    roots = mimeticFEM.EGaussQuad(p);
    %roots = ChebyshevRoots(p)';
    
    % auxiliary matrix to select roots for each polynomial
    rootSelection = true(p+1,p+1);
    rootSelection(linspace(1,(p+1)^2,p+1)) = false;
    
    % Compute each polynomial n using the formula:
    % \frac{\prod_{i=1\\i\neq n}^{p+1}(x-r_{i})}{\prod_{i=1\\i\neq
    % j}^{p+1}(r_{n}-r_{i})}
    %
    % For the top product one uses the built in function poly, based upon
    % the roots. For the bottom part, one just computed the product.
%     if p ~= 1,
%         for n=1:p+1,
%             repmatRoots = repmat(roots(rootSelection(n,:)),[1 sizeOfx]);
%             result(n,:) = prod((repmat(x,[p,1])-repmatRoots)./ (roots(n)-repmatRoots));
%         end
%     else
%         for n=1:p+1,
%             repmatRoots = repmat(roots(rootSelection(n,:)),[1 sizeOfx]);
%             result(n,:) = (repmat(x,[p,1])-repmatRoots)./ (roots(n)-repmatRoots);
%         end
%     end
    
    % fast implementation using Carlo Castoldi's code:
    % http://www.mathworks.com/matlabcentral/fileexchange/899-lagrange-polynomial-interpolation
    for i=1:p+1
      for j=1:p+1
         if (i~=j)
            result(i,:)=result(i,:).*(x-roots(j))/(roots(i)-roots(j));
         end
      end
    end
end
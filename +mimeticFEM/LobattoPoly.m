function result = LobattoPoly(x,p)
%LobattoPoly Returns the p+1 Lobatto Lagrange interpolant basis functions,
% evaluated at x.
%
% It returns a 2-dimensional matrix with the values of the Lobatto basis
% interpolants evaluated at x.
%
% If x is a vector of length N it returns a 2d matrix whose rows are the 
% values of the evaluated polynomial, P(x), in x:
%             -                                               -
%   result = | P_{1}(x(1))   P_{1}(x(2))   ...   P_{1}(x(N))   |
%            | P_{2}(x(1))   P_{2}(x(2))   ...   P_{2}(x(N))   |
%            |                   ...                           |
%            | P_{p+1}(x(1)) P_{p+1}(x(2)) ...   P_{p+1}(x(N)) |
%             -                                               -  
%
% If x=[] then it computes the Lobatto polynomial basis at the Lobatto nodes,
% that is, the result is a sparse identity matrix.
%
%   USAGE
%   -----
%       result = LobattoPoly(x,p)
%
%       Gives the p+1 Lobatto basis interpolants evaluated at x.
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
%   
%
%   Copyright 2009 Artur Palha

%   Revisions:  2009-11-25 (apalha) First implementation.
%               2014-12-03 (apalha) Used the transpose in the for loop and
%                                   then returned the transpose. This
%                                   speeds ups the copmutation by a factor
%                                   of almost 2.
    
    sizeOfx = size(x);
    if sizeOfx(1) == 1,
        sizeOfx = sizeOfx(2);
    elseif sizeOfx(2) == 1,
        sizeOfx = sizeOfx(1);
        x = x';
    end
    
    if isempty(x)
        result = speye(p+1);
        return
    end
    
    % allocate memory space for the result
%    result = zeros([p+1 sizeOfx]);
    %result = ones([p+1 sizeOfx]);
    result = ones([sizeOfx p+1]);
    % compute lobatto roots
    roots = mimeticFEM.LobattoQuad(p);
    %roots = ChebyshevRoots(p)';
    
    % auxiliary matrix to select roots for each polynomial
    %rootSelection = true(p+1,p+1);
    %rootSelection(linspace(1,(p+1)^2,p+1)) = false;
    
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
            result(:,i)=result(:,i).*(x'-roots(j))/(roots(i)-roots(j));
         end
      end
    end
    result = result';
end
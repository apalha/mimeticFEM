function [epsilon_1] = HDivBasisPrimal(xi, eta, p, tensorGrid)
%HDivBasisPrimal Returns 2 x p x (p+1) 2D vector basis of HCurl primal space
%
%
% Computes the 2D HDiv primal space basis functions by tensor product of 1D
% edge polynomial basis functions and Lobatto polynomial basis functions
% at the points defined by the
% coordinates xi and eta. If tensorGrid is set to true then xi and eta are
% the xi and eta coordinates of a tensor grid. If tensorGrid is false then
% the xi and eta are the coordinate of the points where to evaluate the
% basis functions.
%
%   USAGE
%   -----
%       [epsilon_1] = HDivBasisPrimal(xi, eta, p, false)
%
%       Gives the 2 x p x (p+1) primal basis functions of HCurl space, evaluated
%       at the points (xi(k), eta(k)).
%
%       [epsilon_1] = HDivBasisPrimal(xi, eta, p, true)
%
%       Gives the 2 x p x (p+1) primal basis functions of HCurl space, evaluated
%       at the points (xi(i), eta(j)).
%
%   
%   INPUTS
%   ------
%       xi         :: xi coordinates of the points where to evaluate the basis
%                     functions.
%                     (type: float64, size: Nxi)
%       eta        :: eta coordinates of the points where to evaluate the basis
%                     functions.
%                     (type: float64, size: Neta)
%       p          :: polynomial degree of basis
%                     (type: int32, size: single value)
%       tensorGrid :: defines if xi and eta values are a tensor grid (true) or
%                     single points (false). If false Nxi must be equal to Neta
%                     (type: bool, size: single value)
%
%   OUTPUTS
%   -------
%       epsilon_1 :: evaluation of the primal HDiv basis functions at the points
%                    - [xi(k), eta(k)] if tensorGrid = true
%                    - [xi(i), eta(j)] if tensorGrid = false
%                       epsilon(i, j) corresponds to the basis 
%                          i = (i_xi-1) x p + i_eta if i <= p x (p+1) --> dxi basis
%                          i = (i_xi-1) x (p+1) + i_eta if i > p x (p+1) --> deta basis
%                       and node
%                          j
%                    (type: float64, size: tensorGrid = true  => [2xpx(p+1), Nxi x Neta]
%                                          tensorGrid = false => [2xpx(p+1), Nxi = Neta])
%
%   Copyright 2020 Artur Palha

%   Revisions:  2020-02-15 (apalha) First implementation.

    if tensorGrid  % use the nodes as tensor grid
        % Compute the HDiv basis functions based on tensor product of edge basis
        % and Lobatto basis over the p+1 Gauss-Lobatto-Legendre nodes in the interval [-1, 1].
        % epsilon_1 will contain the evaluation of all 2 x p x (p+1) basis functions of degree p
        % as a matrix. Each row is associated to a basis, i.e. epsilon_2(j,:) are
        % the values of the basis function j evaluated at the points xiPlot(:).
        nNodes = length(xi(:))^2;
        nBasis = 2*p*(p+1);
        epsilon_1 = zeros(nBasis, nNodes);
        epsilon_1(1:p*(p+1), :) = kron(mimeticFEM.LobattoPoly(xi, p), mimeticFEM.EdgePoly(eta, p, 'Lobatto'));  % dxi basis
        epsilon_1((p*(p+1)+1):end, :) = kron(mimeticFEM.EdgePoly(xi, p, 'Lobatto'), mimeticFEM.LobattoPoly(eta, p));  % deta basis
    
    else % use the nodes as single nodes
        nNodes = length(xi(:));
        nBasis = 2*p*(p+1);
        epsilon_1 = zeros(nBasis, nNodes);
        
        % dxi basis
        epsilon_1_xi = mimeticFEM.LobattoPoly(eta(:), p);
        epsilon_1_eta = mimeticFEM.EdgePoly(xi(:), p, 'Lobatto');
        
        for kXi = 1:(p+1)
            for kEta = 1:p
                kBasis = (kXi-1)*p + kEta;
                epsilon_1(kBasis, :) = epsilon_1_xi(kXi, :) .* epsilon_1_eta(kEta, :);
            end
        end
        
        % deta basis
        epsilon_1_xi = mimeticFEM.EdgePoly(eta(:), p, 'Lobatto');
        epsilon_1_eta = mimeticFEM.LobattoPoly(xi(:), p);
        
        for kXi = 1:p
            for kEta = 1:(p+1)
                kBasis = (kXi-1)*(p+1) + kEta + p * (p+1);
                epsilon_1(kBasis, :) = epsilon_1_xi(kXi, :) .* epsilon_1_eta(kEta, :);
            end
        end
    end
    
end
function [epsilon_0] = H1BasisPrimal(xi, eta, p, tensorGrid)
%H1BasisPrimal Returns (p+1)x(p+1) 2D basis of H1 primal space
%
%
% Computes the 2D H1 primal space basis functions by tensor product of 1D
% Lobatto polynomial basis functions at the points defines by the
% coordinates xi and eta. If tensorGrid is set to true then xi and eta are
% the xi and eta coordinates of a tensor grid. If tensorGrid is false then
% the xi and eta are the coordinate of the points where to evaluate the
% basis functions.
%
%   USAGE
%   -----
%       [epsilon_0] = H1BasisPrimal(xi, eta, p, false)
%
%       Gives the (p+1)x(p+1) primal basis functions of H1 space, evaluated
%       at the points (xi(k), eta(k)).
%
%       [epsilon_0] = H1BasisPrimal(xi, eta, p, true)
%
%       Gives the (p+1)x(p+1) primal basis functions of H1 space, evaluated
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
%       epsilon_0 :: evaluation of the primal H1 basis functions at the points
%                    - [xi(k), eta(k)] if tensorGrid = true
%                    - [xi(i), eta(j)] if tensorGrid = false
%                       epsilon(i, j) corresponds to the basis 
%                          i = (i_xi-1)x(p+1) + i_eta
%                       and node
%                          j 
%                    (type: float64, size: tensorGrid = true  => [(p+1)x(p+1), Nxi x Neta]
%                                          tensorGrid = false => [(p+1)x(p+1), Nxi = Neta])
%
%   Copyright 2020 Artur Palha

%   Revisions:  2020-02-15 (apalha) First implementation.

    if tensorGrid  % use the nodes as tensor grid
        % Compute the H1 basis functions based on tensor product of Lagrange interpolants over
        % the p+1 Gauss-Lobatto-Legendre nodes in the interval [-1, 1].
        % epsilon_0 will contain the evaluation of all p+1 basis functions of degree p
        % as a matrix. Each row is associated to a basis, i.e. epsilon_0(j,:) are
        % the values of the basis function j evaluated at the points xiPlot(:).
        epsilon_0 = kron(mimeticFEM.LobattoPoly(xi, p), mimeticFEM.LobattoPoly(eta, p));
    
    else % use the nodes as single nodes
        nNodes = length(xi(:));
        nBasis = (p+1)*(p+1);
        epsilon_0 = zeros(nBasis, nNodes);
        
        epsilon_0_xi = mimeticFEM.LobattoPoly(xi(:), p);
        epsilon_0_eta = mimeticFEM.LobattoPoly(eta(:), p);
        
        for kXi = 1:(p+1)
            for kEta = 1:(p+1)
                kBasis = (kXi-1)*(p+1) + kEta;
                epsilon_0(kBasis, :) = epsilon_0_xi(kXi, :) .* epsilon_0_eta(kEta, :);
            end
        end
    end
    
end
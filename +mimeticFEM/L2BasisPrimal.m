function [epsilon_2] = L2BasisPrimal(xi, eta, p, tensorGrid)
%L2BasisPrimal Returns p x p 2D basis of L2 primal space
%
%
% Computes the 2D L2 primal space basis functions by tensor product of 1D
% edge polynomial basis functions at the points defined by the
% coordinates xi and eta. If tensorGrid is set to true then xi and eta are
% the xi and eta coordinates of a tensor grid. If tensorGrid is false then
% the xi and eta are the coordinate of the points where to evaluate the
% basis functions.
%
%   USAGE
%   -----
%       [epsilon_2] = L2BasisPrimal(xi, eta, p, false)
%
%       Gives the p x p primal basis functions of L2 space, evaluated
%       at the points (xi(k), eta(k)).
%
%       [epsilon_2] = L2BasisPrimal(xi, eta, p, true)
%
%       Gives the p x p primal basis functions of L2 space, evaluated
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
%       epsilon_2 :: evaluation of the primal L2 basis functions at the points
%                    - [xi(k), eta(k)] if tensorGrid = true
%                    - [xi(i), eta(j)] if tensorGrid = false
%                       epsilon(i, j) corresponds to the basis 
%                          i = (i_xi-1) x p + i_eta
%                       and node
%                          j
%                    (type: float64, size: tensorGrid = true  => [p x p, Nxi x Neta]
%                                          tensorGrid = false => [p x p, Nxi = Neta])
%
%   Copyright 2020 Artur Palha

%   Revisions:  2020-02-15 (apalha) First implementation.

    if tensorGrid  % use the nodes as tensor grid
        % Compute the L2 basis functions based on tensor product of edge basis over
        % the p+1 Gauss-Lobatto-Legendre nodes in the interval [-1, 1].
        % epsilon_2 will contain the evaluation of all pxp basis functions of degree p
        % as a matrix. Each row is associated to a basis, i.e. epsilon_0(j,:) are
        % the values of the basis function j evaluated at the points xiPlot(:).
        epsilon_2 = kron(mimeticFEM.EdgePoly(xi, p, 'Lobatto'), mimeticFEM.EdgePoly(eta, p, 'Lobatto'));
    
    else % use the nodes as single nodes
        nNodes = length(xi(:));
        nBasis = p*p;
        epsilon_2 = zeros(nBasis, nNodes);
        
        epsilon_2_xi = mimeticFEM.EdgePoly(xi(:), p, 'Lobatto');
        epsilon_2_eta = mimeticFEM.EdgePoly(eta(:), p, 'Lobatto');
        
        for kXi = 1:(p+1)
            for kEta = 1:(p+1)
                kBasis = (kXi-1)*(p+1) + kEta;
                epsilon_2(kBasis, :) = epsilon_2_xi(kXi, :) .* epsilon_2_eta(kEta, :);
            end
        end
    end
    
end
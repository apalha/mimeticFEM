%% mFEM_04_inner_products_2D
% 
% mimeticFEM tutorial demonstrating the computation of inner products
% between basis functions
%
%   Copyright 2020 Artur Palha

%   Revisions:  2020-02-15 (apalha) First implementation.

%% Input parameters
p = 3;  % the polynomial degree of the basis functions
        % note that this is the polynomial degree of the nodal basis
        % functions, the edge basis functions will be one degree lower
        % this occurs in order to have a proper de Rham sequence
pInt = 4;  % degree of quadrature rule

                                             
%% Inner product between H1 basis functions

%% Using Gauss Quadrature

% Compute Gauss nodes and quadrature weights
[gaussNodes, gaussWeights] = mimeticFEM.GaussQuad(pInt);
gaussWeights2D = kron(gaussWeights, gaussWeights'); % convert the 1D weights into 2D weights
gaussWeightsMatrix = spdiags(reshape(kron(gaussWeights, gaussWeights'), [], 1),...
                         0, (pInt+1)^2, (pInt+1)^2);  % convert the 2D weights matrix into a sparse diagonal matrix
                                                      % to simplifly the computation of the inner product

% Evaluate the H1 basis at the gaussNodes
% No optimization is done here, we could separate the basis into the xi and
% eta components and perform two 1D integrals, which would be faster
epsilon_0 = mimeticFEM.H1BasisPrimal(gaussNodes, gaussNodes, p, true);

% Compute the inner product
M_0_gauss = epsilon_0 * gaussWeightsMatrix * epsilon_0';

%% Using Lobatto Quadrature

% Compute Lobatto nodes and quadrature weights
[lobattoNodes, lobattoWeights] = mimeticFEM.LobattoQuad(pInt);
lobattoWeights2D = kron(lobattoWeights, lobattoWeights'); % convert the 1D weights into 2D weights
lobattoWeightsMatrix = spdiags(reshape(kron(lobattoWeights, lobattoWeights'), [], 1),...
                         0, (pInt+1)^2, (pInt+1)^2);  % convert the 2D weights matrix into a sparse diagonal matrix
                                                      % to simplifly the computation of the inner product

% Evaluate the H1 basis at the lobattoNodes
% No optimization is done here, we could separate the basis into the xi and
% eta components and perform two 1D integrals, which would be faster
epsilon_0 = mimeticFEM.H1BasisPrimal(lobattoNodes, lobattoNodes, p, true);

% Compute the inner product
M_0_lobatto = epsilon_0 * lobattoWeightsMatrix * epsilon_0';

%% Difference between Gauss and Lobatto Quadrature
% This different is zero when the quadratures become exact integration,
% this occurs when pInt > p for Gauss quadrature and pInt > p+1 
% for Lobatto quadrature 
figure()
surf(M_0_lobatto - M_0_gauss)
shading interp




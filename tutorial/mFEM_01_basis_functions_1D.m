%% mFEM_01_basis_functions_1D
% 
% mimeticFEM tutorial demonstrating the use of available basis functions in
% 1D
%
% Evaluates the 1D basis functions and plots them for the canonical element
% in the interval [-1, 1]. For plotting on the physical domain,
% proper reconstruction is needed and the Jacobians needs to be taken into
% consideration. This will be in another tutorial.
%
%   Copyright 2020 Artur Palha

%   Revisions:  2020-02-14 (apalha) First implementation.


%% Input parameters
p = 3;  % the polynomial degree of the basis functions
        % note that this is the polynomial degree of the nodal basis
        % functions, the edge basis functions will be one degree lower
        % this occurs in order to have a proper de Rham sequence
        
nPlotPoints = 1001;  % the number of points used to plot the basis functions

%% Pre-computations
xiPlot = linspace(-1, 1, nPlotPoints);  % the points where to evaluate and
                                        % plot the basis functions

%% Nodal basis functions

% Compute the nodal basis functions based on Lagrange interpolants over
% the p+1 Gauss-Lobatto-Legendre nodes in the interval [-1, 1].
% epsilon_0 will contain the evaluation of all p+1 basis functions of degree p
% as a matrix. Each row is associated to a basis, i.e. epsilon_0(j,:) are
% the values of the basis function j evaluated at the points xiPlot(:).
epsilon_0 = mimeticFEM.LobattoPoly(xiPlot, p);

% Plot the basis functions all in one plot
figure()
plot(xiPlot, epsilon_0, '-')
title('Nodal basis functions')
xlabel('\xi')
ylabel('\epsilon^{0}(\xi)')
ylim([-0.25, 1.25])

% Add the Gauss-Lobatto-Legendre nodes just for illustration
xiGLL = mimeticFEM.LobattoQuad(p);

hold on
plot(xiGLL, 1, 'xk')


%% Edge basis functions

% Compute the edge basis functions based on Lagrange interpolants over
% the p+1 Gauss-Lobatto-Legendre nodes in the interval [-1, 1].
% epsilon_1 will contain the evaluation of all p basis functions of degree
% p-1 as a matrix. Each row is associated to a basis, i.e. epsilon_1(j,:) are
% the values of the basis function j evaluated at the points xiPlot(:).
polyType = 'Lobatto';  % this is the type of polynomial (sub-grid) that is used
                       % to generate the edge basis functions. The idea is
                       % that in 1D these edge basis functions and the
                       % associated nodal basis functions form a discrete de
                       % Rham sequence. For this reason if a Lobatto
                       % short for Gauss-Lobatto-Legendre) polynomial is
                       % used for the nodal basis then the edge basis must
                       % be constructed as 'Lobatto' type. The polynomial
                       % degree p is the same because internally the
                       % degrees as matched to that if a sequence of degree
                       % p is chosen we simply need to use p on both nodal
                       % and edge basis, it is simply a convention of the
                       % code.
epsilon_1 = mimeticFEM.EdgePoly(xiPlot, p, polyType);

% Plot the basis functions all in one plot
figure()
plot(xiPlot, epsilon_1, '-')
title('Edge basis functions')
xlabel('\xi')
ylabel('\epsilon^{1}(\xi)')

% Add the Gauss-Lobatto-Legendre nodes just for illustration
xiGLL = mimeticFEM.LobattoQuad(p);

hold on
% plot(xiGLL, 1, 'xk')

% Add the interval associated to each basis function
% Over the interval [xi_{k}, xi_{k+1}] (with xi_{k} the kth-GLL node) the
% basis functions have the kronecker delta property:
% int_{xi_{k}}^{xi_{k+1}} \epsilon_{m}^{1}(\xi) d\xi = \delta _{k,m}
yPlotBounds = ylim;
for kPatch = 1:p
    faceAlpha = mod(kPatch, 2) * 0.1 + mod(kPatch + 1, 2) * 0.15;
    xiCoordinates = rectpulse(xiGLL(kPatch:(kPatch+1)), 2);
    yCoordinates = [yPlotBounds, fliplr(yPlotBounds)];
   patch(xiCoordinates, yCoordinates, 'black', 'EdgeColor', 'none', 'FaceAlpha', faceAlpha) ;
end
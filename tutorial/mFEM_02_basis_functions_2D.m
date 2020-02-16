%% mFEM_02_basis_functions_2D
% 
% mimeticFEM tutorial demonstrating the use of available basis functions in
% 2D
%
% Evaluates the 2D basis functions and plots them for the canonical element
% in the interval [-1, 1]x[-1,1]. For plotting on the physical domain,
% proper reconstruction is needed and the Jacobians needs to be taken into
% consideration. This will be in another tutorial.
%
% We plot here:
%   - basis functions associated to scalar fields in H^1
%   - basis functions associated to scalar fields in H(curl)
%   - basis functions associated to scalar fields in H(div)
%   - basis functions assoiated to scalar fields in L^2
%
%   Copyright 2020 Artur Palha

%   Revisions:  2020-02-15 (apalha) First implementation.

%% Input parameters
p = 3;  % the polynomial degree of the basis functions
        % note that this is the polynomial degree of the nodal basis
        % functions, the edge basis functions will be one degree lower
        % this occurs in order to have a proper de Rham sequence

% Number of plot points for scalar basis functions
nPlotPointsScalar_xi = 201;  % the number of points used to plot the basis functions in xi-direction
nPlotPointsScalar_eta = 201;  % the number of points used to plot the basis functions in eta-direction

% Number of plot points for vector basis functions
nPlotPointsVector_xi = 21;  % the number of points used to plot the basis functions in xi-direction
nPlotPointsVector_eta = 21;  % the number of points used to plot the basis functions in eta-direction

plot_offset = 2.25;  % the offset in xi and eta directions between each plot


%% Pre-computations
xiPlotScalar = linspace(-1, 1, nPlotPointsScalar_xi);  % the points where to evaluate and
                                                       % plot the basis functions in the xi
                                                       % direction
etaPlotScalar = linspace(-1, 1, nPlotPointsScalar_eta);  % the points where to evaluate and
                                                         % plot the basis functions in the
                                                         % eta direction
                                        
xiPlotVector = linspace(-1, 1, nPlotPointsVector_xi);  % the points where to evaluate and
                                                       % plot the basis functions in the xi
                                                       % direction
etaPlotVector = linspace(-1, 1, nPlotPointsVector_eta);  % the points where to evaluate and
                                                         % plot the basis functions in the
                                                         % eta direction

                                             
%% H1 basis functions

% Compute the H1 basis functions based on tensor product of Lagrange interpolants over
% the p+1 Gauss-Lobatto-Legendre nodes in the interval [-1, 1].
% epsilon_0 will contain the evaluation of all (p+1) x (p+1) basis functions of degree p
% as a matrix. Each row is associated to a basis, i.e. epsilon_0(j,:) are
% the values of the basis function j evaluated at the points xiPlot(:).
epsilon_0 = mimeticFEM.H1BasisPrimal(xiPlotScalar, etaPlotScalar, p, true);

% Plot the basis functions all in one plot
figure()
for kXi = 1:p+1
    for kEta = 1:p+1
        kBasis = (kXi - 1) * (p + 1) + kEta;
        xiPlotOffset = plot_offset * (kXi - 1);
        etaPlotOffset = plot_offset * (kEta - 1);
        surf(xiPlotScalar + xiPlotOffset, etaPlotScalar + etaPlotOffset, reshape(epsilon_0(kBasis, :), nPlotPointsScalar_xi, nPlotPointsScalar_eta))
        shading interp
        hold on
    end
end

title('H^{1} basis functions')
xlabel('\xi')
ylabel('\eta')
zlabel('\epsilon^{0}(\xi, \eta)')
view([-35 80])
set(gca,'xtick', (0:p+1) * plot_offset)
set(gca,'xticklabel', 1:p+1)
set(gca,'ytick',(0:p+1) * plot_offset)
set(gca,'yticklabel', 1:p+1)


%% HCurl basis functions

% Compute the Hcurl basis functions based on tensor product of Lagrange interpolants over
% the p+1 Gauss-Lobatto-Legendre nodes and edge basis in the interval [-1, 1].
% epsilon_1 will contain the evaluation of all 2 x p x (p+1) basis functions of degree p
% as a matrix. Each row is associated to a basis, i.e. epsilon_1(j,:) are
% the values of the basis function j evaluated at the points xiPlot(:).
% Since this is a vector basis, we can say that the first px(p+1) basis are
% in the xi direction and the remaining px(p+1) basis are in the eta
% direction, this becomes slightly different for non-orthogonal elements.
[epsilon_1] = mimeticFEM.HCurlBasisPrimal(xiPlotVector, etaPlotVector, p, true);

% Plot the dxi basis functions all in one plot
figure()
[xiPlotGrid, etaPlotGrid] = meshgrid(xiPlotVector, etaPlotVector);
for kXi = 1:p
    for kEta = 1:p+1
        % Plot the basis
        kBasis = (kXi - 1) * (p + 1) + kEta;
        xiPlotOffset = plot_offset * (kXi - 1);
        etaPlotOffset = plot_offset * (kEta - 1);
        quiver(xiPlotGrid(:) + xiPlotOffset, etaPlotGrid(:) + etaPlotOffset, epsilon_1(kBasis, :)', zeros(size(epsilon_1(kBasis, :)')), 'color', 'k')
        hold on
        
        % Plot a box around
        plot([-1 + xiPlotOffset, 1 + xiPlotOffset, 1 + xiPlotOffset, -1 + xiPlotOffset, -1 + xiPlotOffset],...
             [-1 + etaPlotOffset, -1 + etaPlotOffset, 1 + etaPlotOffset, 1 + etaPlotOffset, -1 + etaPlotOffset], '-k')
    end
end

title('HCurl basis functions (xi-direction)')
xlabel('\xi')
ylabel('\eta')
axis equal
set(gca,'xtick', (0:p+1) * plot_offset)
set(gca,'xticklabel', 1:p+1)
set(gca,'ytick',(0:p+1) * plot_offset)
set(gca,'yticklabel', 1:p+1)

% Plot the deta basis functions all in one plot
figure()
[xiPlotGrid, etaPlotGrid] = meshgrid(xiPlotVector, etaPlotVector);
for kXi = 1:p
    for kEta = 1:p+1
        % Plot the basis
        kBasis = (kXi - 1) * (p + 1) + kEta + p * (p+1);
        xiPlotOffset = plot_offset * (kXi - 1);
        etaPlotOffset = plot_offset * (kEta - 1);
        quiver(xiPlotGrid(:) + xiPlotOffset, etaPlotGrid(:) + etaPlotOffset, zeros(size(epsilon_1(kBasis, :)')), epsilon_1(kBasis, :)', 'color', 'k')
        hold on
        
        % Plot a box around
        plot([-1 + xiPlotOffset, 1 + xiPlotOffset, 1 + xiPlotOffset, -1 + xiPlotOffset, -1 + xiPlotOffset],...
             [-1 + etaPlotOffset, -1 + etaPlotOffset, 1 + etaPlotOffset, 1 + etaPlotOffset, -1 + etaPlotOffset], '-k')
    end
end

title('HCurl basis functions (eta-direction)')
xlabel('\xi')
ylabel('\eta')
axis equal
set(gca,'xtick', (0:p+1) * plot_offset)
set(gca,'xticklabel', 1:p+1)
set(gca,'ytick',(0:p+1) * plot_offset)
set(gca,'yticklabel', 1:p+1)


%% HDiv basis functions

% Compute the HDiv basis functions based on tensor product of Lagrange interpolants over
% the p+1 Gauss-Lobatto-Legendre nodes and edge basis in the interval [-1, 1].
% epsilon_1 will contain the evaluation of all 2 x p x (p+1) basis functions of degree p
% as a matrix. Each row is associated to a basis, i.e. epsilon_1(j,:) are
% the values of the basis function j evaluated at the points xiPlot(:).
% Since this is a vector basis, we can say that the first px(p+1) basis are
% in the xi direction and the remaining px(p+1) basis are in the eta
% direction, this becomes slightly different for non-orthogonal elements.
[epsilon_1] = mimeticFEM.HDivBasisPrimal(xiPlotVector, etaPlotVector, p, true);

% Plot the dxi basis functions all in one plot
figure()
[xiPlotGrid, etaPlotGrid] = meshgrid(xiPlotVector, etaPlotVector);
for kXi = 1:(p+1)
    for kEta = 1:p
        % Plot the basis
        kBasis = (kXi - 1) * p + kEta;
        xiPlotOffset = plot_offset * (kXi - 1);
        etaPlotOffset = plot_offset * (kEta - 1);
        quiver(xiPlotGrid(:) + xiPlotOffset, etaPlotGrid(:) + etaPlotOffset, epsilon_1(kBasis, :)', zeros(size(epsilon_1(kBasis, :)')), 'color', 'k')
        hold on
        
        % Plot a box around
        plot([-1 + xiPlotOffset, 1 + xiPlotOffset, 1 + xiPlotOffset, -1 + xiPlotOffset, -1 + xiPlotOffset],...
             [-1 + etaPlotOffset, -1 + etaPlotOffset, 1 + etaPlotOffset, 1 + etaPlotOffset, -1 + etaPlotOffset], '-k')
    end
end

title('HDiv basis functions (xi-direction)')
xlabel('\xi')
ylabel('\eta')
axis equal
set(gca,'xtick', (0:p+1) * plot_offset)
set(gca,'xticklabel', 1:p+1)
set(gca,'ytick',(0:p+1) * plot_offset)
set(gca,'yticklabel', 1:p+1)

% Plot the deta basis functions all in one plot
figure()
[xiPlotGrid, etaPlotGrid] = meshgrid(xiPlotVector, etaPlotVector);
for kXi = 1:p
    for kEta = 1:(p+1)
        % Plot the basis
        kBasis = (kXi - 1) * (p + 1) + kEta + p*(p+1);
        xiPlotOffset = plot_offset * (kXi - 1);
        etaPlotOffset = plot_offset * (kEta - 1);
        quiver(xiPlotGrid(:) + xiPlotOffset, etaPlotGrid(:) + etaPlotOffset, zeros(size(epsilon_1(kBasis, :)')), epsilon_1(kBasis, :)', 'color', 'k')
        hold on
        
        % Plot a box around
        plot([-1 + xiPlotOffset, 1 + xiPlotOffset, 1 + xiPlotOffset, -1 + xiPlotOffset, -1 + xiPlotOffset],...
             [-1 + etaPlotOffset, -1 + etaPlotOffset, 1 + etaPlotOffset, 1 + etaPlotOffset, -1 + etaPlotOffset], '-k')
    end
end

title('HDiv basis functions (eta-direction)')
xlabel('\xi')
ylabel('\eta')
axis equal
set(gca,'xtick', (0:p+1) * plot_offset)
set(gca,'xticklabel', 1:p+1)
set(gca,'ytick',(0:p+1) * plot_offset)
set(gca,'yticklabel', 1:p+1)

%% L2 basis functions

% Compute the L2 basis functions based on tensor product of Lagrange interpolants over
% the p+1 Gauss-Lobatto-Legendre nodes in the interval [-1, 1].
% epsilon_2 will contain the evaluation of all pxp basis functions of
% degree p-1 as a matrix. Each row is associated to a basis, i.e. epsilon_2(j,:) are
% the values of the basis function j evaluated at the points xiPlot(:).
epsilon_2 = mimeticFEM.L2BasisPrimal(xiPlotScalar, etaPlotScalar, p, true);

% Plot the basis functions all in one plot
figure()
for kXi = 1:p
    for kEta = 1:p
        kBasis = (kXi - 1) * p + kEta;
        xiPlotOffset = plot_offset * (kXi - 1);
        etaPlotOffset = plot_offset * (kEta - 1);
        surf(xiPlotScalar + xiPlotOffset, etaPlotScalar + etaPlotOffset, reshape(epsilon_2(kBasis, :), nPlotPointsScalar_xi, nPlotPointsScalar_eta))
        shading interp
        hold on
    end
end

title('L^{2} basis functions')
xlabel('\xi')
ylabel('\eta')
zlabel('\epsilon^{0}(\xi, \eta)')
view([-35 80])
set(gca,'xtick', (0:p) * plot_offset)
set(gca,'xticklabel', 1:p)
set(gca,'ytick',(0:p) * plot_offset)
set(gca,'yticklabel', 1:p)
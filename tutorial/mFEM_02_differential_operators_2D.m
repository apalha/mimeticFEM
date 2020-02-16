%% mFEM_03_differential_operators_2D
% 
% mimeticFEM tutorial demonstrating the use of available differential
% operators in 2D
%
% Evaluates the exterior derivative of all 2D basis functions 
% and plots them for the canonical element
% in the interval [-1, 1]x[-1,1]. For plotting on the physical domain,
% proper reconstruction is needed and the Jacobians needs to be taken into
% consideration. This will be in another tutorial.
%
% We plot here:
%   - exterior derivative (gradient) of basis functions associated to scalar fields in H^1
%   - exterior derivative (curl) of basis functions associated to scalar fields in H(curl)
%   - exterior derivative (divergence) of basis functions associated to scalar fields in H(div)
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

                                             
%% Exterior derivative of H1 basis functions (for the section) H1 -\nabla-> Hcurl


% Compute the exterior derivative of  H1 basis functions
% To do this we use the property:
%   (E^{1,0})^{T} epsilon_{1} = epsilon_{0}
% Where E^{1,0} is the incidence matrix from nodes to edges, i.e. the
% discrete exterior derivative for 0-forms
E10 = mimeticFEM.dZero(p);  % the incidence matrix E^{10}
epsilon_1 = mimeticFEM.HCurlBasisPrimal(xiPlotVector, etaPlotVector, p, true);  

d_epsilon_0_dxi = E10(1:p*(p+1), :)' * epsilon_1(1:p*(p+1), :);
d_epsilon_0_deta = E10((p*(p+1)+1):end, :)' * epsilon_1((p*(p+1)+1):end, :);

figure()
[xiPlotGrid, etaPlotGrid] = meshgrid(xiPlotVector, etaPlotVector);
for kXi = 1:(p+1)
    for kEta = 1:p+1
        % Plot the basis
        kBasis = (kXi - 1) * (p + 1) + kEta;
        xiPlotOffset = plot_offset * (kXi - 1);
        etaPlotOffset = plot_offset * (kEta - 1);
        quiver(xiPlotGrid(:) + xiPlotOffset, etaPlotGrid(:) + etaPlotOffset, d_epsilon_0_dxi(kBasis, :)', d_epsilon_0_deta(kBasis, :)', 'color', 'k')
        hold on
        
        % Plot a box around
        plot([-1 + xiPlotOffset, 1 + xiPlotOffset, 1 + xiPlotOffset, -1 + xiPlotOffset, -1 + xiPlotOffset],...
             [-1 + etaPlotOffset, -1 + etaPlotOffset, 1 + etaPlotOffset, 1 + etaPlotOffset, -1 + etaPlotOffset], '-k')
    end
end

title('\nabla H1 basis functions')
xlabel('\xi')
ylabel('\eta')
axis equal
set(gca,'xtick', (0:p+1) * plot_offset)
set(gca,'xticklabel', 1:p+1)
set(gca,'ytick',(0:p+1) * plot_offset)
set(gca,'yticklabel', 1:p+1)


%% Exterior derivative of H1 basis functions (for the section) H1 -d-> Hdiv


% Compute the exterior derivative of  H1 basis functions
% To do this we use the property:
%   (E^{1,0})^{T} epsilon_{1} = epsilon_{0}
% With some swap of axis
% Where E^{1,0} is the incidence matrix from nodes to edges, i.e. the
% discrete exterior derivative for 0-forms
E10 = mimeticFEM.dZero(p);  % the incidence matrix E^{10}
epsilon_1 = mimeticFEM.HDivBasisPrimal(xiPlotVector, etaPlotVector, p, true);  

d_epsilon_0_dxi = E10((p*(p+1)+1):end, :)' * epsilon_1(1:p*(p+1), :);
d_epsilon_0_deta = -E10(1:p*(p+1), :)' * epsilon_1((p*(p+1)+1):end, :);

figure()
[xiPlotGrid, etaPlotGrid] = meshgrid(xiPlotVector, etaPlotVector);
for kXi = 1:(p+1)
    for kEta = 1:p+1
        % Plot the basis
        kBasis = (kXi - 1) * (p + 1) + kEta;
        xiPlotOffset = plot_offset * (kXi - 1);
        etaPlotOffset = plot_offset * (kEta - 1);
        quiver(xiPlotGrid(:) + xiPlotOffset, etaPlotGrid(:) + etaPlotOffset, d_epsilon_0_dxi(kBasis, :)', d_epsilon_0_deta(kBasis, :)', 'color', 'k')
        hold on
        
        % Plot a box around
        plot([-1 + xiPlotOffset, 1 + xiPlotOffset, 1 + xiPlotOffset, -1 + xiPlotOffset, -1 + xiPlotOffset],...
             [-1 + etaPlotOffset, -1 + etaPlotOffset, 1 + etaPlotOffset, 1 + etaPlotOffset, -1 + etaPlotOffset], '-k')
    end
end

title('\nabla\times H1 basis functions')
xlabel('\xi')
ylabel('\eta')
axis equal
set(gca,'xtick', (0:p+1) * plot_offset)
set(gca,'xticklabel', 1:p+1)
set(gca,'ytick',(0:p+1) * plot_offset)
set(gca,'yticklabel', 1:p+1)

%% div of HDiv basis functions

% Compute the divergence of HDiv basis functions
% To do this we use the property:
%   (E^{2,1})^{T} epsilon_{2} = epsilon_{1}
% Where E^{1,0} is the incidence matrix from nodes to edges, i.e. the
% discrete exterior derivative for 1-forms
E21 = mimeticFEM.dOne(p);  % the incidence matrix E^{21}
epsilon_2 = mimeticFEM.L2BasisPrimal(xiPlotScalar, etaPlotScalar, p, true);  

d_epsilon_1_dxi = E21(:, 1:(p*(p+1)))' * epsilon_2;
d_epsilon_1_deta = E21(:, (p*(p+1)+1):end)' * epsilon_2;

% Plot the divergence of the dxi basis functions all in one plot
figure()
for kXi = 1:p
    for kEta = 1:p
        kBasis = (kXi - 1) * p + kEta;
        xiPlotOffset = plot_offset * (kXi - 1);
        etaPlotOffset = plot_offset * (kEta - 1);
        surf(xiPlotScalar + xiPlotOffset, etaPlotScalar + etaPlotOffset, reshape(d_epsilon_1_dxi(kBasis, :), nPlotPointsScalar_xi, nPlotPointsScalar_eta))
        shading interp
        hold on
    end
end

title('Div Hdiv dxi basis functions')
xlabel('\xi')
ylabel('\eta')
zlabel('\nabla\cdot\epsilon^{1}(\xi, \eta)')
view([-35 80])
set(gca,'xtick', (0:p) * plot_offset)
set(gca,'xticklabel', 1:p)
set(gca,'ytick',(0:p) * plot_offset)
set(gca,'yticklabel', 1:p)

% Plot the divergence of the deta basis functions all in one plot
figure()
for kXi = 1:p
    for kEta = 1:p
        kBasis = (kXi - 1) * p + kEta;
        xiPlotOffset = plot_offset * (kXi - 1);
        etaPlotOffset = plot_offset * (kEta - 1);
        surf(xiPlotScalar + xiPlotOffset, etaPlotScalar + etaPlotOffset, reshape(d_epsilon_1_deta(kBasis, :), nPlotPointsScalar_xi, nPlotPointsScalar_eta))
        shading interp
        hold on
    end
end

title('Div Hdiv deta basis functions')
xlabel('\xi')
ylabel('\eta')
zlabel('\nabla\cdot\epsilon^{1}(\xi, \eta)')
view([-35 80])
set(gca,'xtick', (0:p) * plot_offset)
set(gca,'xticklabel', 1:p)
set(gca,'ytick',(0:p) * plot_offset)
set(gca,'yticklabel', 1:p)


%% curl of HCurl basis functions

% Compute the divergence of HDiv basis functions
% To do this we use the property:
%   (E^{2,1})^{T} epsilon_{2} = epsilon_{1}
% Where E^{1,0} is the incidence matrix from nodes to edges, i.e. the
% discrete exterior derivative for 1-forms
E21 = mimeticFEM.dOne(p);  % the incidence matrix E^{21}
epsilon_2 = mimeticFEM.L2BasisPrimal(xiPlotScalar, etaPlotScalar, p, true);  

d_epsilon_1_dxi = E21(:, 1:(p*(p+1)))' * epsilon_2;
d_epsilon_1_deta = E21(:, (p*(p+1)+1):end)' * epsilon_2;

% Plot the divergence of the dxi basis functions all in one plot
figure()
for kXi = 1:p
    for kEta = 1:p
        kBasis = (kXi - 1) * p + kEta;
        xiPlotOffset = plot_offset * (kXi - 1);
        etaPlotOffset = plot_offset * (kEta - 1);
        surf(xiPlotScalar + xiPlotOffset, etaPlotScalar + etaPlotOffset, reshape(d_epsilon_1_dxi(kBasis, :), nPlotPointsScalar_xi, nPlotPointsScalar_eta))
        shading interp
        hold on
    end
end

title('Div Hdiv dxi basis functions')
xlabel('\xi')
ylabel('\eta')
zlabel('\nabla\cdot\epsilon^{1}(\xi, \eta)')
view([-35 80])
set(gca,'xtick', (0:p) * plot_offset)
set(gca,'xticklabel', 1:p)
set(gca,'ytick',(0:p) * plot_offset)
set(gca,'yticklabel', 1:p)

% Plot the divergence of the deta basis functions all in one plot
figure()
for kXi = 1:p
    for kEta = 1:p
        kBasis = (kXi - 1) * p + kEta;
        xiPlotOffset = plot_offset * (kXi - 1);
        etaPlotOffset = plot_offset * (kEta - 1);
        surf(xiPlotScalar + xiPlotOffset, etaPlotScalar + etaPlotOffset, reshape(d_epsilon_1_deta(kBasis, :), nPlotPointsScalar_xi, nPlotPointsScalar_eta))
        shading interp
        hold on
    end
end

title('Div Hdiv deta basis functions')
xlabel('\xi')
ylabel('\eta')
zlabel('\nabla\cdot\epsilon^{1}(\xi, \eta)')
view([-35 80])
set(gca,'xtick', (0:p) * plot_offset)
set(gca,'xticklabel', 1:p)
set(gca,'ytick',(0:p) * plot_offset)
set(gca,'yticklabel', 1:p)
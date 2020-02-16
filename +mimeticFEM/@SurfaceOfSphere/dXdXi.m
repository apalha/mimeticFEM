function dXdXiEvaluated = dXdXi(obj,elements,xi,eta)
% dXdXi computes the dXdXi term of the metric for
% at the local positions (xi,eta) for the element. For the SurfaceOfSphere grid.
%
%   USAGE
%   -----
%       dXdXiEvaluated = obj.dXdXi(elements,xi,eta)
%
%   INPUTS
%   ------
%       elements :: the elements numbers
%                   If elements is a scalar then x, y, and z have dimensions
%                   (N,M), if element is a vector of dimension Q then
%                   x, y, and z have dimensions (N,M,Q), where (N,M) are the
%                   dimensions of xi and eta. If element = [] then the
%                   mapping is computed for all elements
%                   (type: int32, size: vector)
%       xi :: the local (xi) coordinates of the nodes where to
%             evaluate the mapping.
%             (type: float64, size: array [N,M])
%       eta :: the local (eta) coordinates of the nodes where to
%              evaluate the mapping.
%              (type: float64, size: array [N,M])
%
%
%   OUTPUTS
%   -------
%       dXdXiEvaluated  :: the evaluation of the xi derivative of the mapping
%                          between the computational and physical spaces at 
%                          the points (xi,eta) of the element.
%                          (type: float64, size: array [N,M] (Q=1) or [N,M,Q] (Q>1))
%
%   Copyright 2009-2018 Artur Palha

%   Revisions:  2018-06-27 First implementation. (apalha)
    
    % if element is empty [] then all the mapping is computed for all
    % elements
    if isempty(elements)
        elements = 1:obj.numElements;
    end
    
    % if elements is an array it must be reshaped so that the x and y
    % coordinates are an array with three indices (i,j,k). i,j correspond
    % to the indices of the xi and eta arrays and the k index to the
    % elements.
    if ~isscalar(elements)
        elements = reshape(elements, 1, 1, []);
    end
    
    % First check to which subdomain each element belongs to
    elementsSubDomain = ceil(elements/(obj.n(1)*obj.n(2)));
    
    %% Compute the mapping from the canonical domain to the non-affine quad on the surfaces of the cube
    
    % compute the (s,t) coordinates of the element
    s = (xi+1.0)*0.5;
    t = (eta+1.0)*0.5;
    
    % evaluate the boundary mappings
    gamma1xs = bsxfun(@plus, bsxfun(@times, s, (obj.nodesOfElements(2, 1, elements) - obj.nodesOfElements(1, 1, elements))), obj.nodesOfElements(1, 1, elements));
    gamma1ys = bsxfun(@plus, bsxfun(@times, s, (obj.nodesOfElements(2, 2, elements) - obj.nodesOfElements(1, 2, elements))), obj.nodesOfElements(1, 2, elements));
    gamma2xt = bsxfun(@plus, bsxfun(@times, t, (obj.nodesOfElements(3, 1, elements) - obj.nodesOfElements(2, 1, elements))), obj.nodesOfElements(2, 1, elements));
    gamma2yt = bsxfun(@plus, bsxfun(@times, t, (obj.nodesOfElements(3, 2, elements) - obj.nodesOfElements(2, 2, elements))), obj.nodesOfElements(2, 2, elements));
    gamma3xs = bsxfun(@plus, bsxfun(@times, s, (obj.nodesOfElements(3, 1, elements) - obj.nodesOfElements(4, 1, elements))), obj.nodesOfElements(4, 1, elements));
    gamma3ys = bsxfun(@plus, bsxfun(@times, s, (obj.nodesOfElements(3, 2, elements) - obj.nodesOfElements(4, 2, elements))), obj.nodesOfElements(4, 2, elements));
    gamma4xt = bsxfun(@plus, bsxfun(@times, t, (obj.nodesOfElements(4, 1, elements) - obj.nodesOfElements(1, 1, elements))), obj.nodesOfElements(1, 1, elements));
    gamma4yt = bsxfun(@plus, bsxfun(@times, t, (obj.nodesOfElements(4, 2, elements) - obj.nodesOfElements(1, 2, elements))), obj.nodesOfElements(1, 2, elements));

    
    gamma1x0 = obj.nodesOfElements(1, 1, elements);
    gamma1y0 = obj.nodesOfElements(1, 2, elements);
    gamma1x1 = obj.nodesOfElements(2, 1, elements);
    gamma1y1 = obj.nodesOfElements(2, 2, elements);
    gamma3x0 = obj.nodesOfElements(4, 1, elements);
    gamma3y0 = obj.nodesOfElements(4, 2, elements);
    gamma3x1 = obj.nodesOfElements(3, 1, elements);
    gamma3y1 = obj.nodesOfElements(3, 2, elements);
    
    % evaluate the derivative of the boundary mappings
    % maybe need to multiply by 1/2 all these terms because the mapping is
    % first [-1, 1] --> [0, 1]
    dgamma1xs_ds = obj.nodesOfElements(2, 1, elements) - obj.nodesOfElements(1, 1, elements);
    dgamma1ys_ds = obj.nodesOfElements(2, 2, elements) - obj.nodesOfElements(1, 2, elements);
    dgamma2xt_dt = obj.nodesOfElements(3, 1, elements) - obj.nodesOfElements(2, 1, elements);
    dgamma2yt_dt = obj.nodesOfElements(3, 2, elements) - obj.nodesOfElements(2, 2, elements);
    dgamma3xs_ds = obj.nodesOfElements(3, 1, elements) - obj.nodesOfElements(4, 1, elements);
    dgamma3ys_ds = obj.nodesOfElements(3, 2, elements) - obj.nodesOfElements(4, 2, elements);
    dgamma4xt_dt = obj.nodesOfElements(4, 1, elements) - obj.nodesOfElements(1, 1, elements);
    dgamma4yt_dt = obj.nodesOfElements(4, 2, elements) - obj.nodesOfElements(1, 2, elements);

    
    % compute the mapping between the canonical domain and the element
    % inside the subdomain in 2D
    xi_subdomains = bsxfun(@times, (1.0-s), gamma4xt) + bsxfun(@times, s, gamma2xt) + bsxfun(@times, (1.0-t), gamma1xs) + bsxfun(@times, t, gamma3xs) -...
        bsxfun(@times, (1.0-s), (bsxfun(@times, (1.0-t), gamma1x0) + bsxfun(@times, t, gamma3x0))) - bsxfun(@times, s, (bsxfun(@times, (1.0-t), gamma1x1) + bsxfun(@times, t, gamma3x1)));
    
    eta_subdomains = bsxfun(@times, (1.0-s), gamma4yt) + bsxfun(@times, s, gamma2yt) + bsxfun(@times, (1.0-t), gamma1ys) + bsxfun(@times, t,gamma3ys) -...
        bsxfun(@times, (1.0-s), (bsxfun(@times, (1.0-t), gamma1y0) + bsxfun(@times, t, gamma3y0))) - bsxfun(@times, s, (bsxfun(@times, (1.0-t), gamma1y1) + bsxfun(@times, t, gamma3y1))); 

    %% Compute dX2dXi: the derivative of the mapping from the subdomain to the surface of the sphere
    
    dX2dXiEvaluated = zeros(size(xi_subdomains));
    
    elementsToComputeIndex = 1:length(elements);
    
    % subdomain 1
    elementsInSubdomain_1 = elementsToComputeIndex(elementsSubDomain == 1);
    dX2dXiEvaluated(:, :, elementsInSubdomain_1) = sqrt(0.5 - (eta_subdomains(:, :, elementsInSubdomain_1).^2)/6);
    
    % subdomain 2
    elementsInSubdomain_2 = elementsToComputeIndex(elementsSubDomain == 2);
    dX2dXiEvaluated(:, :, elementsInSubdomain_2) = 0.5 * (-xi_subdomains(:, :, elementsInSubdomain_2) +...
        (2/3)*xi_subdomains(:, :, elementsInSubdomain_2).*(eta_subdomains(:, :, elementsInSubdomain_2).^2)) ./...
        sqrt(1.0 - 0.5*(xi_subdomains(:, :, elementsInSubdomain_2).^2) - 0.5*(eta_subdomains(:, :, elementsInSubdomain_2).^2) +...
        (((xi_subdomains(:, :, elementsInSubdomain_2).*eta_subdomains(:, :, elementsInSubdomain_2)).^2)/3));
    
    % subdomain 3
    elementsInSubdomain_3 = elementsToComputeIndex(elementsSubDomain == 3);
    dX2dXiEvaluated(:, :, elementsInSubdomain_3) = -(1/6)*eta_subdomains(:, :, elementsInSubdomain_3).*xi_subdomains(:, :, elementsInSubdomain_3)./...
        sqrt(0.5 - (xi_subdomains(:, :, elementsInSubdomain_3).^2)/6);
    
    % subdomain 4
    elementsInSubdomain_4 = elementsToComputeIndex(elementsSubDomain == 4);
    dX2dXiEvaluated(:, :, elementsInSubdomain_4) = -(1/6)*eta_subdomains(:, :, elementsInSubdomain_4).*xi_subdomains(:, :, elementsInSubdomain_4)./...
        sqrt(0.5 - (xi_subdomains(:, :, elementsInSubdomain_4).^2)/6);
    
    % subdomain 5
    elementsInSubdomain_5 = elementsToComputeIndex(elementsSubDomain == 5);
    dX2dXiEvaluated(:, :, elementsInSubdomain_5) = sqrt(0.5 - (eta_subdomains(:, :, elementsInSubdomain_5).^2)/6);
    
    % subdomain 6
    elementsInSubdomain_6 = elementsToComputeIndex(elementsSubDomain == 6);
    dX2dXiEvaluated(:, :, elementsInSubdomain_6) = -0.5 * (-xi_subdomains(:, :, elementsInSubdomain_6) +...
        (2/3)*xi_subdomains(:, :, elementsInSubdomain_6).*(eta_subdomains(:, :, elementsInSubdomain_6).^2)) ./...
        sqrt(1.0 - 0.5*(xi_subdomains(:, :, elementsInSubdomain_6).^2) - 0.5*(eta_subdomains(:, :, elementsInSubdomain_6).^2) +...
        (((xi_subdomains(:, :, elementsInSubdomain_6).*eta_subdomains(:, :, elementsInSubdomain_6)).^2)/3));
    
    
    %% Compute dX2dEta: the derivative of the mapping from the subdomain to the surface of the sphere
    
    dX2dEtaEvaluated = zeros(size(xi_subdomains));
    
    elementsToComputeIndex = 1:length(elements);
    
    % subdomain 1
    elementsInSubdomain_1 = elementsToComputeIndex(elementsSubDomain == 1);
    dX2dEtaEvaluated(:, :, elementsInSubdomain_1) = -(1/6)*eta_subdomains(:, :, elementsInSubdomain_1).*xi_subdomains(:, :, elementsInSubdomain_1)./...
        sqrt(0.5 - (eta_subdomains(:, :, elementsInSubdomain_1).^2)/6);
    
    % subdomain 2
    elementsInSubdomain_2 = elementsToComputeIndex(elementsSubDomain == 2);
    dX2dEtaEvaluated(:, :, elementsInSubdomain_2) = 0.5 * (-eta_subdomains(:, :, elementsInSubdomain_2) +...
        (2/3)*eta_subdomains(:, :, elementsInSubdomain_2).*(xi_subdomains(:, :, elementsInSubdomain_2).^2)) ./...
        sqrt(1.0 - 0.5*(xi_subdomains(:, :, elementsInSubdomain_2).^2) - 0.5*(eta_subdomains(:, :, elementsInSubdomain_2).^2) +...
        (((xi_subdomains(:, :, elementsInSubdomain_2).*eta_subdomains(:, :, elementsInSubdomain_2)).^2)/3));
    
    % subdomain 3
    elementsInSubdomain_3 = elementsToComputeIndex(elementsSubDomain == 3);
    dX2dEtaEvaluated(:, :, elementsInSubdomain_3) = sqrt(0.5 - (xi_subdomains(:, :, elementsInSubdomain_3).^2)/6);
    
    % subdomain 4
    elementsInSubdomain_4 = elementsToComputeIndex(elementsSubDomain == 4);
    dX2dEtaEvaluated(:, :, elementsInSubdomain_4) = sqrt(0.5 - (xi_subdomains(:, :, elementsInSubdomain_4).^2)/6);
    
    % subdomain 5
    elementsInSubdomain_5 = elementsToComputeIndex(elementsSubDomain == 5);
    dX2dEtaEvaluated(:, :, elementsInSubdomain_5) = -(1/6)*eta_subdomains(:, :, elementsInSubdomain_5).*xi_subdomains(:, :, elementsInSubdomain_5)./...
        sqrt(0.5 - (eta_subdomains(:, :, elementsInSubdomain_5).^2)/6);
    
    % subdomain 6
    elementsInSubdomain_6 = elementsToComputeIndex(elementsSubDomain == 6);
    dX2dEtaEvaluated(:, :, elementsInSubdomain_6) = -0.5 * (-eta_subdomains(:, :, elementsInSubdomain_6) +...
        (2/3)*eta_subdomains(:, :, elementsInSubdomain_6).*(xi_subdomains(:, :, elementsInSubdomain_6).^2)) ./...
        sqrt(1.0 - 0.5*(xi_subdomains(:, :, elementsInSubdomain_6).^2) - 0.5*(eta_subdomains(:, :, elementsInSubdomain_6).^2) +...
        (((xi_subdomains(:, :, elementsInSubdomain_6).*eta_subdomains(:, :, elementsInSubdomain_6)).^2)/3));

    
    %% Compute dX1dXi: the derivative of the mapping from the canonical domain to the surface of the cube
    dX1dXiEvaluated = -gamma4xt + gamma2xt + bsxfun(@times, (1.0 - t), (gamma1x1 - gamma1x0)) + bsxfun(@times, t, (gamma3x1 - gamma3x0)) +...
        (bsxfun(@times, (1.0 - t), gamma1x0) + bsxfun(@times, t, gamma3x0)) - (bsxfun(@times, (1.0 - t), gamma1x1) + bsxfun(@times, t, gamma3x1));
    
    %% Compute dY1dXi: the derivative of the mapping from the canonical domain to the surface of the cube
    dY1dXiEvaluated = -gamma4yt + gamma2yt + bsxfun(@times, (1.0 - t), (gamma1y1 - gamma1y0)) + bsxfun(@times, t, (gamma3y1 - gamma3y0)) +...
        (bsxfun(@times, (1.0 - t), gamma1y0) + bsxfun(@times, t, gamma3y0)) - (bsxfun(@times, (1.0 - t), gamma1y1) + bsxfun(@times, t, gamma3y1));

    
    %% Compute dXdXi
    dXdXiEvaluated = dX2dXiEvaluated .* dX1dXiEvaluated + dX2dEtaEvaluated .* dY1dXiEvaluated;
    
end
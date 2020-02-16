function [x, y, z] = mapping(obj, elements, xi, eta)
% Mapping computes mapping between (xi,eta) and (x,y,z) for the
% SurfaceOfSphere.
%
% Computes the mapping between (xi, eta) and (x, y, z),
% the computational domain and the physical domain, respectively, for a
% given element number (element).
%
% The surface of the sphere is divided into six subdomains which are
% oriented according to the following diagram:
%
%                               G
%                             -->--
%                         -------------
%                        | ?--> eta = y|
%                      | | |           | |
%                    A v | v           | v B
%                      | | xi = z      | |
%                        |             |
%                         -------------
%                        | ?--> eta = y|
%                      | | |           | |
%                    C v | v           | v D
%                 C    | | xi = x      | |    D
%               -->--    |             |    --<--
%           ------------- ------------- -------------
%          |             |             |             |
%        | | eta = z     | eta = z     |      xi = z | |
%      A ? | ?           | ?           |           ? | ? B
%        | | |           | |           |           | | |
%          | ?--> xi = x | ?--> xi = x |eta = x <--? |
%           ------------- ------------- -------------
%               -->--    |             |    --<--
%                 E    | | eta = x     | |    F
%                    E ? | ?           | ? F
%                      | | |           | |
%                        | ?--> xi = y |
%                         -------------
%                             -->--
%                               G
%
% The mapping is done by first mapping each subdomain to a surface of a
% cube. This is a straightforward map. The surface of the cube is then
% mapped to the surface of the sphere by the folowing expression:
%
%    *** Missing expression!!! ***
%
%   USAGE
%   -----
%       [x, y, z] = obj.mapping(element, xi, eta)
%
%       Computes the mapping between the nodes (xi,eta) of the
%       computational domain and the physical domain.
%
%       element can be for example: A, [A1,A2,A3] or [] (all elements)
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
%       x  :: the evaluation of the x component of the mapping function at
%             at the nodes (xi,eta) in the element (element).
%             (type: float64, size: array [N,M] (Q=1) or [N,M,Q] (Q>1))
%       y  :: the evaluation of the y component of the mapping function at
%             at the nodes (xi,eta) in the element (element).
%             (type: float64, size: array [N,M] (Q=1) or [N,M,Q] (Q>1))
%       z  :: the evaluation of the z component of the mapping function at
%             the nodes (xi,eta) in the element (element).
%             (type: float64, size: array [N,M] (Q=1) or [N,M,Q] (Q>1))
%
%   Copyright 2018 Artur Palha

%   Revisions:  2018-05-13 First implementation. (apalha)
    
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
    
    % compute the (s,t) coordinates of the element
    s = (xi+1.0)*0.5;
    t = (eta+1.0)*0.5;
    
    % evaluate the mappings
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
    
    % compute the mapping between the canonical domain and the element
    % inside the subdomain in 2D
    xi_subdomains = bsxfun(@times, (1.0-s), gamma4xt) + bsxfun(@times, s, gamma2xt) + bsxfun(@times, (1.0-t), gamma1xs) + bsxfun(@times, t, gamma3xs) -...
        bsxfun(@times, (1.0-s), (bsxfun(@times, (1.0-t), gamma1x0) + bsxfun(@times, t, gamma3x0))) - bsxfun(@times, s, (bsxfun(@times, (1.0-t), gamma1x1) + bsxfun(@times, t, gamma3x1)));
    
    eta_subdomains = bsxfun(@times, (1.0-s), gamma4yt) + bsxfun(@times, s, gamma2yt) + bsxfun(@times, (1.0-t), gamma1ys) + bsxfun(@times, t,gamma3ys) -...
        bsxfun(@times, (1.0-s), (bsxfun(@times, (1.0-t), gamma1y0) + bsxfun(@times, t, gamma3y0))) - bsxfun(@times, s, (bsxfun(@times, (1.0-t), gamma1y1) + bsxfun(@times, t, gamma3y1))); 

    % map all elements in each subdomain into the surface of the sphere
    x = zeros(size(xi_subdomains));
    y = zeros(size(xi_subdomains));
    z = zeros(size(xi_subdomains));
    
    elementsToComputeIndex = 1:length(elements);
    
    % subdomain 1
    elementsInSubdomain_1 = elementsToComputeIndex(elementsSubDomain == 1);
    [x(:, :, elementsInSubdomain_1), y(:, :, elementsInSubdomain_1), z(:, :, elementsInSubdomain_1)] = obj.Mapping_1(xi_subdomains(:, :, elementsInSubdomain_1), eta_subdomains(:, :, elementsInSubdomain_1));
    
    % subdomain 2
    elementsInSubdomain_2 = elementsToComputeIndex(elementsSubDomain == 2);
    [x(:, :, elementsInSubdomain_2), y(:, :, elementsInSubdomain_2), z(:, :, elementsInSubdomain_2)] = obj.Mapping_2(xi_subdomains(:, :, elementsInSubdomain_2), eta_subdomains(:, :, elementsInSubdomain_2));
    
    % subdomain 3
    elementsInSubdomain_3 = elementsToComputeIndex(elementsSubDomain == 3);
    [x(:, :, elementsInSubdomain_3), y(:, :, elementsInSubdomain_3), z(:, :, elementsInSubdomain_3)] = obj.Mapping_3(xi_subdomains(:, :, elementsInSubdomain_3), eta_subdomains(:, :, elementsInSubdomain_3));
    
    % subdomain 4
    elementsInSubdomain_4 = elementsToComputeIndex(elementsSubDomain == 4);
    [x(:, :, elementsInSubdomain_4), y(:, :, elementsInSubdomain_4), z(:, :, elementsInSubdomain_4)] = obj.Mapping_4(xi_subdomains(:, :, elementsInSubdomain_4), eta_subdomains(:, :, elementsInSubdomain_4));
    
    % subdomain 5
    elementsInSubdomain_5 = elementsToComputeIndex(elementsSubDomain == 5);
    [x(:, :, elementsInSubdomain_5), y(:, :, elementsInSubdomain_5), z(:, :, elementsInSubdomain_5)] = obj.Mapping_5(xi_subdomains(:, :, elementsInSubdomain_5), eta_subdomains(:, :, elementsInSubdomain_5));
    
    % subdomain 6
    elementsInSubdomain_6 = elementsToComputeIndex(elementsSubDomain == 6);
    [x(:, :, elementsInSubdomain_6), y(:, :, elementsInSubdomain_6), z(:, :, elementsInSubdomain_6)] = obj.Mapping_6(xi_subdomains(:, :, elementsInSubdomain_6), eta_subdomains(:, :, elementsInSubdomain_6));
    
end
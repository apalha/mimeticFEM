function [x,y] = Mapping(element,xi,eta,n)
% Mapping computes streak mapping between (xi,eta) and (x,y) using a
% transfinite mapping.
%
% Computes the mapping between (xi,eta) and (x,y),
% the computational domain and the physical domain, respectively, for a
% given element number (element) and the number of elements in x and y
% directions (n).
%
%
%   The mesh corresponds to a subdivision of each region into n(1)xn(2)
%   elements.
%
%   %% TO CHANGE
%   For the element k with 1 <= k <= n(1)*n(2) (the 2D numbering (i,j) goes such
%   that k = (i-1)n(2) + j) the mapping is:
%
%       (\xi,\eta) \in [-1,1]x[-1,1] --> (r,z)
%
%   
%   USAGE
%   -----
%       [x, y] = Mapping(element,xi,eta,n)
%
%       Computes the mapping between the nodes (xi,eta) of the
%       computational domain and the physical domain.
%
%   INPUTS
%   ------
%       element :: the element number
%                  (type: int32, size: single value)
%       xi :: the local (xi) coordinates of the nodes where to
%             evaluate the mapping.
%             (type: float64, size: array [N,M])
%       eta :: the local (eta) coordinates of the nodes where to
%              evaluate the mapping.
%              (type: float64, size: array [N,M])
%       n :: the number of elements in x and y directions
%            (type: int32, size: array [1,2])
%
%
%
%   OUTPUTS
%   -------
%       x  :: the evaluation of the x component of the mapping function at
%             at the nodes (xi,eta) in the element (element).
%             (type: float64, size: array [N,M])
%       y  :: the evaluation of the y component of the mapping function at
%             at the nodes (xi,eta) in the element (element).
%             (type: float64, size: array [N,M])
%
%   REFERENCES
%   ----------
%       
%
%   Copyright 2015 Artur Palha

%   Revisions:  2015-08-05 (apalha) First implementation.
    
    if isempty(element)
        element = 1:(3*prod(n));
    end
    
    % preallocate memory space for x and y
    sizeNodes = size(xi);
    x = zeros(sizeNodes(1),sizeNodes(2),length(element));
    y = zeros(sizeNodes(1),sizeNodes(2),length(element));
    
    for k=1:length(element)
        ix = floor(element(k)/(3*n(2)));
        iy = mod(element(k),3*n(2));
        if iy == 0 % region 1 (top region)
            region = 1;
            element_region = ix*n(2);
        elseif (iy <= n(2)) % region 3 (bottom region)
            region = 3;
            element_region = ix*n(2)+iy;
        elseif (iy <= 2*n(2)) % region 2 (middle region)
            region = 2;
            element_region = ix*n(2)+iy-n(2);
        else % region 1 (top region)
            region = 1;
            element_region = ix*n(2)+iy-2*n(2);
        end

        % define the boundaries
        gamma1 = @(s) mimeticFEM2.streak.gamma1(s,region); % bottom
        gamma2 = @(s) mimeticFEM2.streak.gamma2(s,region); % right
        gamma3 = @(s) mimeticFEM2.streak.gamma3(s,region); % top
        gamma4 = @(s) mimeticFEM2.streak.gamma4(s,region); % left

        gamma = {gamma1,gamma2,gamma3,gamma4};

        % compute the mapping using the transfinite mapping
        [x(:,:,k),y(:,:,k)] = mimeticFEM2.transfinite.Mapping(element_region,xi,eta,n,gamma);
    end

end
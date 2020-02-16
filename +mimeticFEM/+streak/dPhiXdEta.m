function dPhiXdEtaEvaluated = dPhiXdEta(element,xi,eta,n)
% dPhiXdEta computes the dPhiXdEta term of the metric for
% at the local positions (xi,eta) for the element. For the transfinite pataki grid.
%
%   USAGE
%   -----
%       dPhiXdEtaEvaluated = dPhiXdEta(element,xi,eta,n,region)
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
%   OUTPUTS
%   -------
%       dPhiXdEtaEvaluated  :: the evaluation of the eta derivative of the mapping
%                              between the computational and physical spaces at 
%                              the points (xi,eta) of the element.
%                              (type: float64, size: array [N,M])
%
%   Copyright 2015 Artur Palha

%   Revisions:  2015-08-05 (apalha) Updated documentation.
    
    if isempty(element)
        element = 1:(3*prod(n));
    end
    
    % preallocate memory space for dPhiXdEtaEvaluated
    sizeNodes = size(xi);
    dPhiXdEtaEvaluated = zeros(sizeNodes(1),sizeNodes(2),length(element));
    
    for k=1:length(element)
        ix = floor(element(k)/(3*n(2)));
        iy = mod(element(k),3*n(2));
        if iy == 0 % region 1 (top region)
            region = 1;
            element_region = ix*n(2);
        elseif iy <= n(2) % region 3 (bottom region)
            region = 3;
            element_region = ix*n(2)+iy;
        elseif iy <= 2*n(2) % region 2 (middle region)
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

        % define the derivatives of the boundaries
        dgamma1 = @(s) mimeticFEM2.streak.dgamma1ds(s,region); % bottom
        dgamma2 = @(s) mimeticFEM2.streak.dgamma2ds(s,region); % right
        dgamma3 = @(s) mimeticFEM2.streak.dgamma3ds(s,region); % top
        dgamma4 = @(s) mimeticFEM2.streak.dgamma4ds(s,region); % left

        dgamma = {dgamma1,dgamma2,dgamma3,dgamma4};

        % compute the mapping using the transfinite mapping
        dPhiXdEtaEvaluated(:,:,k) = mimeticFEM2.transfinite.dPhiXdEta(element_region,xi,eta,n,gamma,dgamma);
    end
end
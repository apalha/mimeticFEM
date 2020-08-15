function D = dHRotBasisPrimal(p,varargin)
%dHRotBasisPrimal computes the discrete exterior derivative of elements of
%   HRotBasisPrimal.
%
%   USAGE
%   -----
%       D = dHRotBasisPrimal(p)
%
%   INPUTS
%   ------
%       p :: The order of the polynomial basis functions.
%            (type: int32, size: single value)
%       n :: The dimension of the system. For now, 1 (1D) and 2 (2D)
%       -     can be used. (optional, default:2)
%            (type: int32, size: single value)
%
%   OUTPUTS
%   -------  
%       D :: The matrix representing the discrete exterior derivative.
%            (type: float64, size: [2p(p+1) (p+1)(p+1)] sparse)
%
%   Copyright 2010 Artur Palha

%   Revisions:  2020-08-15 (apalha) First implementation.
    
    % Compute first the DGrad and then swap the axis
    if isempty(varargin)
        DGrad = mimeticFEM.dHGradBasisPrimal(p);
    else
        DGrad = mimeticFEM.dHGradBasisPrimal(p, varargin);
    end
    
    % Swap the axis (and change sign)
    D = zeros(size(DGrad));
    D(1:p*(p+1), :) = DGrad((p*(p+1)+1):end, :);
    D((p*(p+1)+1):end, :) = -DGrad(1:p*(p+1), :);
end
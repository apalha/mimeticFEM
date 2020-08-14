function D = dZero(p,varargin)
%dZero computes the discrete exterior derivative of discrete zero forms in
%one spectral element of order p.
%
%   USAGE
%   -----
%       D = dZero(p)
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

%   Revisions:  2010-01-21 (apalha) First implementation.
    
    % check if the dimension of the manifold is given in the optional
    % arguments
    if ~isempty(varargin)
        % if the dimension is given, use it
        n = varargin{1};
    else
        % if not use 2D as the default
        n = 2;
    end
    
    if n == 2 % construct the 2D discrete exterior derivative
        % xi part
        iiXiMinus = 1:(p*(p+1));
        jjXiMinus = 1:(p*(p+1));
        iiXiPlus = iiXiMinus;
        jjXiPlus = (1+p+1):((p+1)*(p+1));

        % eta part
        iiEtaMinus = (1+p*(p+1)):(2*p*(p+1));
        jjEtaMinus = zeros([p*(p+1) 1]);

        for l = 1:(p+1)
            for k = 1:p
                jjEtaMinus((l-1)*p + k) = (l-1)*(p+1) + k;
            end
        end

        iiEtaPlus = iiEtaMinus;
        jjEtaPlus = jjEtaMinus + 1;

        % values
        valuesMinus = -ones([p*(p+1) 1]);
        valuesPlus = ones([p*(p+1) 1]);

        D = sparse([iiXiMinus'; iiXiPlus'; iiEtaMinus'; iiEtaPlus'], [jjXiMinus'; jjXiPlus'; jjEtaMinus; jjEtaPlus], [valuesMinus; valuesPlus; valuesMinus; valuesPlus]);
    
    else % construct the 1D discrete exterior derivative 
        D = spdiags(repmat([-1, 1],p,2),[0 1],p,p+1);
    end
end
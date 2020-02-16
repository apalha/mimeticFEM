function [x, y, z] = Mapping_2(obj, xi, eta)
%MAPPING_2 Maps the canonical domain onto the domain 2 of the unit sphere

    x = sqrt(1.0 - 0.5*(xi.^2) - 0.5*(eta.^2) + ((xi.^2).*(eta.^2))/3.0);
    y = xi .* sqrt(0.5 - (eta.^2)/6.0);
    z = eta .* sqrt(0.5 - (xi.^2)/6.0);
end


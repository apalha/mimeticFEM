function [x, y, z] = Mapping_6(obj, xi, eta)
%MAPPING_6 Maps the canonical domain onto the domain 6 of the unit sphere

    x = -sqrt(1.0 - 0.5*(eta.^2) - 0.5*(xi.^2) + ((xi.^2).*(eta.^2))/3.0);
    y = eta .* sqrt(0.5 - (xi.^2)/6.0);
    z = xi .* sqrt(0.5 - (eta.^2)/6.0);
end


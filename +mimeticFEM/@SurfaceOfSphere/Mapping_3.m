function [x, y, z] = Mapping_3(obj, xi, eta)
%MAPPING_3 Maps the canonical domain onto the domain 3 of the unit sphere

    x = eta .* sqrt(0.5 - (xi.^2)/6.0);
    y = sqrt(1.0 - 0.5*(eta.^2) - 0.5*(xi.^2) + ((xi.^2).*(eta.^2))/3.0);
    z = xi .* sqrt(0.5 - (eta.^2)/6.0);
end


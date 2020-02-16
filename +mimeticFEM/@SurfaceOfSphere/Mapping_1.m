function [x, y, z] = Mapping_1(obj, xi, eta)
%MAPPING_1 Maps the canonical domain onto the element 1 of the unit sphere

    x = xi .* sqrt(0.5 - (eta.^2)/6.0);
    y = -sqrt(1.0 - 0.5*(eta.^2) - 0.5*(xi.^2) + ((xi.^2).*(eta.^2))/3.0);
    z = eta .* sqrt(0.5 - (xi.^2)/6.0);
end


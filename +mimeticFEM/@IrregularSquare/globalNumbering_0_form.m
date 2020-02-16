function [zeroFormGN] = globalNumbering_0_form(obj, p)
% globalNumbering_0_form computes the local to global numbering of 0-forms
%  of polynomial degree p on the SurfaceOfSphere mesh.
%
% The numbering is done along eta first and then xi. According to the
% schematics below:
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
%   USAGE
%   -----
%       [zeroFormGN] = obj.globalNumbering_0_form(p)
%
%       Computes the lobal to global numbering for 0-forms on the
%       SurfaceOfSphere mesh.
%
%   INPUTS
%   ------
%       p :: the polynomial degree of the 0-forms.
%            (type: int32, size: sv)
%
%   OUTPUTS
%   -------
%       zeroFormGN  :: the local to global numbering of zero-forms.
%                      (type: int32, size: [(p+1)^2, obj.numElements])
%
%
%   Copyright 2018 Artur Palha

%   Revisions:  2018-06-17 First implementation. (apalha)
    
    %% Numbering of 0-form
    zeroFormGN = zeros(p+1, p+1, obj.numElements);

    %% Subdomain 1

    numberOfNewDofs = (obj.n(1)*p + 1) * (obj.n(2)*p + 1);
    newDofs = reshape(1:numberOfNewDofs, (obj.n(2)*p + 1), (obj.n(1)*p + 1));

    for element=1:obj.n(1)*obj.n(2)
        iElement = ceil(element/obj.n(2));
        jElement = element - (iElement-1)*obj.n(2);
        iStart = (iElement-1)*p + 1;
        iEnd = iStart + p;
        jStart = (jElement-1)*p + 1;
        jEnd = jStart + p;
        zeroFormGN(:, :, element) = newDofs(jStart:jEnd, iStart:iEnd);
    end

    lastDof = max(zeroFormGN(:));

    %% Subdomain 2

    numberOfNewDofs = (obj.n(1)*p + 1) * (obj.n(2)*p + 1) - (obj.n(2)*p + 1); % remove the dofs on the left because they were already numbered in subdomain 1
    newDofs = zeros((obj.n(2)*p + 1), (obj.n(1)*p + 1));
    newDofs(:,2:end) = reshape(1:numberOfNewDofs, (obj.n(2)*p + 1), obj.n(1)*p) + lastDof;

    for elementInSubdomain=1:obj.n(1)*obj.n(2)
        iElement = ceil(elementInSubdomain/obj.n(2));
        jElement = elementInSubdomain - (iElement-1)*obj.n(2);
        iStart = (iElement-1)*p + 1;
        iEnd = iStart + p;
        jStart = (jElement-1)*p + 1;
        jEnd = jStart + p;
        element = elementInSubdomain + obj.n(1)*obj.n(2);
        zeroFormGN(:, :, element) = newDofs(jStart:jEnd, iStart:iEnd);
    end

    % add the dofs already numbered
    for elementInSubdomain=1:obj.n(2)
        iElement = ceil(elementInSubdomain/obj.n(2));
        jElement = elementInSubdomain - (iElement-1)*obj.n(2);
        elementSlave = elementInSubdomain + obj.n(1)*obj.n(2);
        elementMaster = (obj.n(1)-1)*obj.n(2) + jElement;
        zeroFormGN(:, 1, elementSlave) = zeroFormGN(:, p+1, elementMaster);
    end

    lastDof = max(zeroFormGN(:));


    %% Subdomain 3

    numberOfNewDofs = (obj.n(1)*p + 1) * (obj.n(2)*p + 1) - (obj.n(1)*p + 1); % remove the dofs on the left because they were already numbered in subdomain 1
    newDofs = zeros((obj.n(2)*p + 1), (obj.n(1)*p + 1));
    newDofs(1:end-1,:) = reshape(1:numberOfNewDofs, obj.n(2)*p, (obj.n(1)*p+1)) + lastDof;

    for elementInSubdomain=1:obj.n(1)*obj.n(2)
        iElement = ceil(elementInSubdomain/obj.n(2));
        jElement = elementInSubdomain - (iElement-1)*obj.n(2);
        iStart = (iElement-1)*p + 1;
        iEnd = iStart + p;
        jStart = (jElement-1)*p + 1;
        jEnd = jStart + p;
        element = elementInSubdomain + 2*obj.n(1)*obj.n(2);
        zeroFormGN(:, :, element) = newDofs(jStart:jEnd, iStart:iEnd);
    end

    % add the dofs already numbered
    for iElement=1:obj.n(1)
        jElement = obj.n(2);
        elementInSubdomain = (iElement-1)*obj.n(2) + jElement;
        elementSlave = elementInSubdomain + 2*obj.n(1)*obj.n(2);
        elementMaster = obj.n(1)*obj.n(2) + (obj.n(1)-1)*obj.n(2) + iElement;
        zeroFormGN(p+1, :, elementSlave) = zeroFormGN(:, p+1, elementMaster);
    end

    lastDof = max(zeroFormGN(:));

    %% Subdomain 4

    numberOfNewDofs = (obj.n(1)*p + 1) * (obj.n(2)*p + 1) - 2*(obj.n(2)*p + 1) - (obj.n(1)*p - 1); % remove the dofs on the left because they were already numbered in subdomain 1
    newDofs = zeros((obj.n(2)*p + 1), (obj.n(1)*p + 1));
    newDofs(1:end-1,2:end-1) = reshape(1:numberOfNewDofs, obj.n(2)*p, obj.n(1)*p - 1) + lastDof;

    for elementInSubdomain=1:obj.n(1)*obj.n(2)
        iElement = ceil(elementInSubdomain/obj.n(2));
        jElement = elementInSubdomain - (iElement-1)*obj.n(2);
        iStart = (iElement-1)*p + 1;
        iEnd = iStart + p;
        jStart = (jElement-1)*p + 1;
        jEnd = jStart + p;
        element = elementInSubdomain + 3*obj.n(1)*obj.n(2);
        zeroFormGN(:, :, element) = newDofs(jStart:jEnd, iStart:iEnd);
    end

    % add the dofs already numbered
    for elementInSubdomain=1:obj.n(2)  % left dofs
        iElement = ceil(elementInSubdomain/obj.n(2));
        jElement = elementInSubdomain - (iElement-1)*obj.n(2);
        elementSlave = elementInSubdomain + 3*obj.n(1)*obj.n(2);
        elementMaster = (jElement - 1)*obj.n(2) + 1;
        zeroFormGN(:, 1, elementSlave) = zeroFormGN(1, :, elementMaster);
    end

    for elementInSubdomain=((obj.n(1)-1)*obj.n(2) + 1):obj.n(1)*obj.n(2)  % right dofs
        iElement = ceil(elementInSubdomain/obj.n(2));
        jElement = elementInSubdomain - (iElement-1)*obj.n(2);
        elementSlave = elementInSubdomain + 3*obj.n(1)*obj.n(2);
        elementMaster = 2*obj.n(1)*obj.n(2) + jElement;
        zeroFormGN(:, p+1, elementSlave) = fliplr(zeroFormGN(:, 1, elementMaster));
    end

    for iElement=1:obj.n(1)  % top dofs
        jElement = obj.n(2);
        elementInSubdomain = (iElement - 1) * obj.n(2) + jElement;
        elementSlave = elementInSubdomain + 3*obj.n(1)*obj.n(2);
        elementMaster = obj.n(1)*obj.n(2) + (iElement - 1)*obj.n(2) + 1;
        zeroFormGN(p+1, :, elementSlave) = zeroFormGN(1, :, elementMaster);
    end

    lastDof = max(zeroFormGN(:));


    %% Subdomain 5

    numberOfNewDofs = (obj.n(1)*p + 1) * (obj.n(2)*p + 1) - 2*(obj.n(1)*p + 1) - (obj.n(2)*p - 1); 
    newDofs = zeros((obj.n(2)*p + 1), (obj.n(1)*p + 1));
    newDofs(2:end-1,1:end-1) = reshape(1:numberOfNewDofs, obj.n(2)*p - 1, obj.n(1)*p) + lastDof;

    for elementInSubdomain=1:obj.n(1)*obj.n(2)
        iElement = ceil(elementInSubdomain/obj.n(2));
        jElement = elementInSubdomain - (iElement-1)*obj.n(2);
        iStart = (iElement-1)*p + 1;
        iEnd = iStart + p;
        jStart = (jElement-1)*p + 1;
        jEnd = jStart + p;
        element = elementInSubdomain + 4*obj.n(1)*obj.n(2);
        zeroFormGN(:, :, element) = newDofs(jStart:jEnd, iStart:iEnd);
    end

    % add the dofs already numbered
    for iElement = 1:obj.n(1)  % bottom dofs
        jElement = 1;
        elementInSubdomain = (iElement - 1) * obj.n(2) + jElement;
        elementSlave = elementInSubdomain + 4*obj.n(1)*obj.n(2);
        elementMaster = (iElement - 1)*obj.n(2) + obj.n(2);
        zeroFormGN(1, :, elementSlave) = zeroFormGN(p+1, :, elementMaster);
    end

    for iElement = 1:obj.n(1)  % top dofs
        jElement = obj.n(2);
        elementInSubdomain = (iElement - 1) * obj.n(2) + jElement;
        elementSlave = elementInSubdomain + 4*obj.n(1)*obj.n(2);
        elementMaster = 2*obj.n(1)*obj.n(2) + (obj.n(1)-1)*obj.n(2) + iElement;
        zeroFormGN(p+1, :, elementSlave) = zeroFormGN(:, p+1, elementMaster);
    end

    for jElement=1:obj.n(2)  % right dofs
        iElement = obj.n(1);
        elementInSubdomain = (iElement - 1) * obj.n(2) + jElement;
        elementSlave = elementInSubdomain + 4*obj.n(1)*obj.n(2);
        iElementMaster = jElement;
        jElementMaster = obj.n(2);
        elementMaster = obj.n(1)*obj.n(2) + (iElementMaster - 1)*obj.n(2) + jElementMaster;
        zeroFormGN(:, p+1, elementSlave) = zeroFormGN(p+1, :, elementMaster);
    end

    lastDof = max(zeroFormGN(:));


    %% Subdomain 6

    numberOfNewDofs = (obj.n(1)*p + 1) * (obj.n(2)*p + 1) - 2*(obj.n(1)*p + 1) - 2*(obj.n(2)*p - 1); 
    newDofs = zeros((obj.n(2)*p + 1), (obj.n(1)*p + 1));
    newDofs(2:end-1,2:end-1) = reshape(1:numberOfNewDofs, obj.n(2)*p - 1, obj.n(1)*p - 1) + lastDof;

    for elementInSubdomain=1:obj.n(1)*obj.n(2)
        iElement = ceil(elementInSubdomain/obj.n(2));
        jElement = elementInSubdomain - (iElement-1)*obj.n(2);
        iStart = (iElement-1)*p + 1;
        iEnd = iStart + p;
        jStart = (jElement-1)*p + 1;
        jEnd = jStart + p;
        element = elementInSubdomain + 5*obj.n(1)*obj.n(2);
        zeroFormGN(:, :, element) = newDofs(jStart:jEnd, iStart:iEnd);
    end

    % add the dofs already numbered
    for iElement = 1:obj.n(1)  % bottom dofs
        jElement = 1;
        elementInSubdomain = (iElement - 1) * obj.n(2) + jElement;
        elementSlave = elementInSubdomain + 5*obj.n(1)*obj.n(2);
        elementMaster = iElement;
        zeroFormGN(1, :, elementSlave) = zeroFormGN(:, 1, elementMaster);
    end

    for iElement = 1:obj.n(1)  % top dofs
        jElement = obj.n(2);
        elementInSubdomain = (iElement - 1) * obj.n(2) + jElement;
        elementSlave = elementInSubdomain + 5*obj.n(1)*obj.n(2);
        elementMaster = (iElement-1)*obj.n(2) + 1 + 2*obj.n(1)*obj.n(2);
        zeroFormGN(p+1, :, elementSlave) = zeroFormGN(1, :, elementMaster);
    end


    for jElement = 1:obj.n(2)  % left dofs
        iElement = 1;
        elementInSubdomain = (iElement - 1) * obj.n(2) + jElement;
        elementSlave = elementInSubdomain + 5*obj.n(1)*obj.n(2);
        elementMaster = (jElement-1)*obj.n(2) + 1 + 3*obj.n(1)*obj.n(2);
        zeroFormGN(:, 1, elementSlave) = zeroFormGN(1, :, elementMaster);
    end


    for jElement = 1:obj.n(2)  % right dofs
        iElement = obj.n(1);
        elementInSubdomain = (iElement - 1) * obj.n(2) + jElement;
        elementSlave = elementInSubdomain + 5*obj.n(1)*obj.n(2);
        elementMaster = jElement + 4*obj.n(1)*obj.n(2);
        zeroFormGN(:, p+1, elementSlave) = zeroFormGN(:, 1, elementMaster);
    end

end
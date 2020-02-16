function zeroFormDiscrete = DiscretizeZeroForm(f, mesh, p, gridType)
%DiscretizeZeroForm discretizes a zero form.
%
%   USAGE
%   -----
%       zeroFormDiscrete = DiscretizeZeroForm(f, mesh, p, gridType)
%
%   INPUTS
%   ------
%
%       f           :: The 0-form to discretize. A matlab function: f(x,y),
%                      of f(x,y,z), depends on mesh.dim(2).
%                      (type:matlab function, size: single value)
%       mesh :: (type: mimeticFEM2.Mesh, size: single value)
%               A mesh object over which the 2-form is defined. 
%       p           :: The order of the polynomials used for the
%                      discretization.
%                      (type: int64, size: single value)
%       gridType    :: The type of grid used for the discretization, can
%                      be 'Lobatto' or 'EGauss'
%                      (type: string, size: single value)
%
%   OUTPUTS
%   -------  
%       zeroFormDiscrete :: The coefficients associated to the discrete
%                           zero-form, which correspond to the evaluation
%                           of function f at the nodes of the grid.
%                           (type: float64, size: [(p+1)(p+1), n])
%
%   Copyright 2012-2014 Artur Palha

%   Revisions:  2011-01-01 (apalha) First implementation.
%               2012-04-03 (apalha) Added general mapping.
%               2014-12-05 (apalha) Use only n as a single value to be the
%                                   number of elements.
%               2017-03-05 (apalha) Mapping in the inputs was changed to
%                                   mesh object.
%               2018-06-17 (apalha) Adds discretization of surfaces
%                                   embedded in R^3. mesh.dim is now used.

%    
    % check if gridType is a valid one
    if ~mimeticFEM.TestPolyType(gridType)
        disp(sprintf(':: %s :: is not a valid type of grid', gridType));
        return
    end
    
    % compute the nodes of the grid to use, given the gridType
    gridNodes = eval(sprintf('mimeticFEM.%sQuad(%s)', strtrim(gridType), 'p'));
    
    % compute the meshgrid of nodes
    [xi, eta] = meshgrid(gridNodes);
    
    % allocate memory space for discrete zero-form
    zeroFormDiscrete = zeros((p+1)*(p+1),mesh.numElements);
    
    if mesh.dim(2) == 2
        for element=1:mesh.numElements
            [x, y] = mesh.mapping(element,xi,eta);
            zeroFormDiscrete(:,element) = reshape(f(x,y),[(p+1)*(p+1) 1]);
        end
        
    elseif mesh.dim(2) == 3
        for element=1:mesh.numElements
            [x, y, z] = mesh.mapping(element,xi,eta);
            zeroFormDiscrete(:,element) = reshape(f(x, y, z),[(p+1)*(p+1) 1]);
        end
    end
    
end
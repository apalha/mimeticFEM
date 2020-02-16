function twoFormDiscrete = DiscretizeTwoForm(f, mesh, p, pint, gridType)
%DiscretizeTwoForm discretizes a two form.
%
%   USAGE
%   -----
%       twoFormDiscrete = DiscretizeTwoForm(f, mesh, p, pint, gridType)
%
%   INPUTS
%   ------
%
%       f :: The 2-form to discretize. A matlab function: f(x,y).
%            (type:matlab function, size: single value)
%       mesh :: (type: mimeticFEM2.Mesh, size: single value)
%               A mesh object over which the 2-form is defined.
%       p :: The order of the polynomials used for the
%            discretization.
%            (type: int64, size: single value)
%       pint :: The order of the quadrature used for computing the
%               integrals of the 2-form
%               (type: int64, size: single value)
%       gridType :: The type of grid used for the discretization, can
%                   be 'Lobatto' or 'EGauss'
%                      (type: string, size: single value)
%
%   OUTPUTS
%   -------  
%       twoFormDiscrete :: The coefficients associated to the discrete
%                           two-form, which correspond to the integral of
%                           the function f over the sub-element surfaces grid.
%                           (type: float64, size: [p p, mesh.numElements])
%
%   Copyright 2012-2015 Artur Palha

%   Revisions:  2018-06-17 (apalha) Adds discretization of surfaces
%                                   embedded in R^3. mesh.dim is now used.
%               2017-02-22 (apalha) Mapping, g and n (number of elements)
%                                   inputs changed to mesh.
%               2015-06-09 (apalha) Updated documentation.
%               2014-12-05 (apalha) Use only n as a single value to be the
%                                   number of elements.
%               2012-04-03 (apalha) Added general mapping.
%               2011-01-01 (apalha) First implementation.
%               
    
    % check if gridType is a valid one
    if ~mimeticFEM.TestPolyType(gridType)
        disp(sprintf(':: %s :: is not a valid type of grid', gridType));
        return
    end
    
    % compute the nodes of the grid to use, given the gridType
    gridNodes = eval(sprintf('mimeticFEM.%sQuad(%s)', strtrim(gridType), 'p'));
    
    % compute the nodes and the weights of the quadrature to use to
    % approximate the integrals
    % compute quadrature weights and nodes
    [quadNodes,quadWeights] = mimeticFEM.GaussQuad(pint);
    
    % compute the matrix of 2d weights
    quadWeights2d = kron(quadWeights', quadWeights);
    
    % compute the meshgrid of nodes
    [xi,eta] = meshgrid(quadNodes);
    
    % compute the delimiting nodes of each sub-element
    iIndices = reshape(repmat((1:p), p, 1), [], 1);
%     iIndices = rectpulse((1:p)', p);
    jIndices = repmat((1:p)', [p 1]);
    
    nodeSubElementsLowerLeft = [gridNodes(iIndices) gridNodes(jIndices)]; 
    nodeSubElementsUpperRight = [gridNodes(iIndices+1) gridNodes(jIndices+1)];
    
    % compute the scalling factor of the inner integrals, it is just
    % multiplying by the volume ratio, because it is already straight
    subCellSizes = [(nodeSubElementsUpperRight(:,1)-nodeSubElementsLowerLeft(:,1)),...
                        (nodeSubElementsUpperRight(:,2)-nodeSubElementsLowerLeft(:,2))];
    
    xi = 0.5*repmat(xi(:)+1,[1 p*p])*spdiags(subCellSizes(:,1),0,p*p,p*p) + repmat(nodeSubElementsLowerLeft(:,1)',[(pint+1)*(pint+1) 1]);
    eta = 0.5*repmat(eta(:)+1,[1 p*p])*spdiags(subCellSizes(:,2),0,p*p,p*p) + repmat(nodeSubElementsLowerLeft(:,2)',[(pint+1)*(pint+1) 1]);
    
    % compute the integral for each sub-element
    
    % allocate memory space for the result of the integral
    twoFormDiscrete = zeros(p*p,mesh.numElements);
    
    if mesh.dim(2) == 2 
        if mesh.numElements>1
            for element=1:mesh.numElements
                [x,y] = mesh.mapping(element,xi,eta);
                evaluatedg = mesh.g(element,xi,eta);
                twoFormDiscrete(:,element) = quadWeights2d(:)'*(f(x,y).*evaluatedg)*(spdiags(subCellSizes(:,1).*subCellSizes(:,2)*0.25,0,p*p,p*p));
            end
        else
            [x,y] = mesh.mapping(1,xi,eta);
            evaluatedg = mesh.g(1,xi,eta);
            twoFormDiscrete = (quadWeights2d(:)'*(f(x,y).*evaluatedg)*(spdiags(subCellSizes(:,1).*subCellSizes(:,2)*0.25,0,p*p,p*p)))';
        end
    elseif mesh.dim(2) == 3
        if mesh.numElements>1
            for element=1:mesh.numElements
                [x, y, z] = mesh.mapping(element,xi,eta);
                evaluatedg = mesh.g(element,xi,eta);
                twoFormDiscrete(:,element) = quadWeights2d(:)'*(f(x, y, z).*evaluatedg)*(spdiags(subCellSizes(:,1).*subCellSizes(:,2)*0.25,0,p*p,p*p));
            end
        else
            [x, y, z] = mesh.mapping(1,xi,eta);
            evaluatedg = mesh.g(1,xi,eta);
            twoFormDiscrete = (quadWeights2d(:)'*(f(x,y,z).*evaluatedg)*(spdiags(subCellSizes(:,1).*subCellSizes(:,2)*0.25,0,p*p,p*p)))';
        end
    end
    
end
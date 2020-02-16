function oneFormDiscrete = DiscretizeOneForm(f, mesh, p, pint, gridType,n)
%DiscretizeOneForm discretized a 1-form.
%
%   twoFormDiscrete = DiscretizeOneForm(f, mesh, p, pint, gridType,n)
%
%   Where:
%
%       f           :: the 1-form to discretize (matlab function)
%       mesh :: (type: mimeticFEM2.Mesh, size: single value)
%               A mesh object over which the 2-form is defined.
%       p           :: the order of the discretization
%       pint        :: the order used in the evaluation of the integrals
%                      for discretization
%       gridType    :: the type of grid used for the discretization, can
%                      be 'Lobatto' or 'EGauss'
%       n           :: the number of elements
%
%   It returns a vector: twoFormDiscrete.
%
%   Each element of the vector oneFormDiscrete (referred here as fd) is one of 
%   the discrete components of the 1-form f. The discretization is
%   implemented in the following way:
%
%   fd_{i} = \int_{eta_{k}}^{eta_{k+1}}\int_{\xi_{n}}^{\xi^{n+1}} (f o phi)
%   (\xi,\eta) d\xi d\eta
%
%   where: i = (n-1)*p + k,   l,k = 1,...,p
%
%   Note that this integral is done numerically, that is, a Gauss
%   quadrature of order qint is used in both \xi and \eta directions.

%   Copyright 2011 Artur Palha
%   $Revision: 1.0 $  $Date: 2011/11/07 $
    
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
    [quadNodes quadWeights] = mimeticFEM.GaussQuad(pint);
    
    
    % compute the scalling factor of the inner integrals, it is just
    % multiplying by the volume ratio, because it is already straight
    subEdgeSizes = gridNodes(2:end) - gridNodes(1:(end-1));
    
    % compute the integration nodes coordinates in xi
    xiXi = reshape(repmat(repmat(0.5*(quadNodes(:)+1),[1 p])*spdiags(subEdgeSizes(:),0,p,p) + repmat(gridNodes(1:(end-1))',[pint+1 1]),[p+1 1]),[pint+1 p*(p+1)]);
    etaXi = repmat(gridNodes(:)',[pint+1 p]);
    
    % compute the integration nodes coordinates in eta
    xiEta = repmat(rectpulse(gridNodes(:)',p),[pint+1 1]);
    etaEta = repmat(repmat(0.5*(quadNodes(:)+1),[1 p])*spdiags(subEdgeSizes(:),0,p,p) + repmat(gridNodes(1:(end-1))',[(pint+1) 1]),[1 p+1]);
     
    % compute the integral for each sub-element
    
    % allocate memory space for the result of the integral
    oneFormDiscrete = zeros(2*p*(p+1),mesh.numElements);
%    twoFormDiscreteAlt = zeros(p*p,nElements);
    
    if mesh.numElements>1
        for element=1:mesh.numElements
            % the xi component
            
            [x y] = mesh.mapping(element,xiXi,etaXi); % the real coordinates of the quadrature nodes
            
            dPhiXdXiEvaluated = mesh.dXdXi(element,xiXi,etaXi); % the metric terms
            dPhiYdXiEvaluated = mesh.dYdXi(element,xiXi,etaXi); % the metric terms
            dPhiXdEtaEvaluated = mesh.dXdEta(element,xiXi,etaXi); % the metric terms
            dPhiYdEtaEvaluated = mesh.dYdEta(element,xiXi,etaXi); % the metric terms
            
            % compute the form at the integration nodes
            [myFX, myFY] = f(x,y);
            
            oneFormDiscrete(1:(p*(p+1)),element) = quadWeights(:)'*(myFX.*dPhiXdXiEvaluated + myFY.*dPhiYdXiEvaluated)*(spdiags(rectpulse(subEdgeSizes(:),p+1)*0.5,0,p*(p+1),p*(p+1)));
            
            % the eta component
            
            [x, y] = mesh.mapping(element,xiEta,etaEta); % the real coordinates of the quadrature nodes
            
            dPhiXdXiEvaluated = mesh.dXdXi(element,xiEta,etaEta); % the metric terms
            dPhiYdXiEvaluated = mesh.dYdXi(element,xiEta,etaEta); % the metric terms
            dPhiXdEtaEvaluated = mesh.dXdEta(element,xiEta,etaEta); % the metric terms
            dPhiYdEtaEvaluated = mesh.dYdEta(element,xiEta,etaEta); % the metric terms
            
            % compute the form at the integration nodes
            [myFX, myFY] = f(x,y);
            
            oneFormDiscrete(((p*(p+1))+1):end,element) = quadWeights(:)'*(myFX.*dPhiXdEtaEvaluated + myFY.*dPhiYdEtaEvaluated)*(spdiags(repmat(subEdgeSizes(:),[p+1 1])*0.5,0,p*(p+1),p*(p+1)));

%             for n=1:(p*p)
%                 [x y] = myPhi(0.5*(xi+1.0)*subCellSizes(n,1) + nodeSubElementsLowerLeft(n,1),...
%                             0.5*(eta+1.0)*subCellSizes(n,2) + nodeSubElementsLowerLeft(n,2));
%                 evaluatedg = myG(0.5*(xi+1.0)*subCellSizes(n,1) + nodeSubElementsLowerLeft(n,1),...
%                              0.5*(eta+1.0)*subCellSizes(n,2) + nodeSubElementsLowerLeft(n,2));
%                 evaluatedweightsfg = (quadWeights2d.*f(x,y).*evaluatedg);
%                 twoFormDiscrete(n,element) = sum(evaluatedweightsfg(:))*subCellSizes(n,1)*subCellSizes(n,2)*0.25;
%             end
        end
    else
        % the xi component

        [x, y] = mesh.mapping(1,xiXi,etaXi); % the real coordinates of the quadrature nodes

        dPhiXdXiEvaluated = mesh.dXdXi(1,xiXi,etaXi); % the metric terms
        dPhiYdXiEvaluated = mesh.dYdXi(1,xiXi,etaXi); % the metric terms
        dPhiXdEtaEvaluated = mesh.dXdEta(1,xiXi,etaXi); % the metric terms
        dPhiYdEtaEvaluated = mesh.dYdEta(1,xiXi,etaXi); % the metric terms

        % compute the form at the integration nodes
        [myFX, myFY] = f(x,y);

        oneFormDiscrete(1:(p*(p+1)),1) = quadWeights(:)'*(myFX.*dPhiXdXiEvaluated + myFY.*dPhiYdXiEvaluated)*(spdiags(rectpulse(subEdgeSizes(:),p+1)*0.5,0,p*(p+1),p*(p+1)));

        % the eta component

        [x, y] = mesh.mapping(1,xiEta,etaEta); % the real coordinates of the quadrature nodes

        dPhiXdXiEvaluated = mesh.dXdXi(1,xiEta,etaEta); % the metric terms
        dPhiYdXiEvaluated = mesh.dYdXi(1,xiEta,etaEta); % the metric terms
        dPhiXdEtaEvaluated = mesh.dXdEta(1,xiEta,etaEta); % the metric terms
        dPhiYdEtaEvaluated = mesh.dYdEta(1,xiEta,etaEta); % the metric terms

        % compute the form at the integration nodes
        [myFX, myFY] = f(x,y);

        oneFormDiscrete(((p*(p+1))+1):end,1) = quadWeights(:)'*(myFX.*dPhiXdEtaEvaluated + myFY.*dPhiYdEtaEvaluated)*(spdiags(repmat(subEdgeSizes(:),[p+1 1])*0.5,0,p*(p+1),p*(p+1)));

%             for n=1:(p*p)
%                 [x y] = myPhi(0.5*(xi+1.0)*subCellSizes(n,1) + nodeSubElementsLowerLeft(n,1),...
%                             0.5*(eta+1.0)*subCellSizes(n,2) + nodeSubElementsLowerLeft(n,2));
%                 evaluatedg = myG(0.5*(xi+1.0)*subCellSizes(n,1) + nodeSubElementsLowerLeft(n,1),...
%                              0.5*(eta+1.0)*subCellSizes(n,2) + nodeSubElementsLowerLeft(n,2));
%                 evaluatedweightsfg = (quadWeights2d.*f(x,y).*evaluatedg);
%                 twoFormDiscrete(n,element) = sum(evaluatedweightsfg(:))*subCellSizes(n,1)*subCellSizes(n,2)*0.25;
%             end
       
    end
    
end
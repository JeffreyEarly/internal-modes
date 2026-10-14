classdef (Abstract) IMSolver
    % Define the shared protocol for canonical EVP solvers.
    %
    % Concrete solvers own the grid, coordinate mapping, derivative
    % matrices, integration rule, and interpolation of native modes. The
    % base class owns the common generalized-eigenvalue workflow.
    %
    % - Topic: Solve EVPs
    % - Topic: Solve geostrophic zero-APV modes
    % - Topic: Developer topics
    % - Declaration: classdef (Abstract) IMSolver

    methods
        function [z,weights] = nativeQuadratureRule(self,zBounds)
            % Return the solver's native physical quadrature rule.
            %
            % The returned points are increasing in physical coordinate and
            % `weights` are reordered with them. For a WKB-configured
            % `IMSolverSpectral`, the points are Chebyshev--Lobatto points in
            % $$x(z)=\int N(z)\,dz$$ and the weights act directly on values
            % sampled in physical $$z$$.
            %
            % - Topic: Solve EVPs
            % - Declaration: [z,weights] = nativeQuadratureRule(solver,zBounds)
            % - Parameter zBounds: physical integration bounds
            % - Returns z: increasing physical quadrature points
            % - Returns weights: physical-coordinate quadrature weights
            arguments
                self IMSolver
                zBounds (1,2) double
            end

            nativeZ = self.innerProductGrid(zBounds);
            nativeWeights = self.innerProductWeights(nativeZ,zBounds);
            [z,index] = sort(nativeZ(:));
            weights = nativeWeights(index);
        end

        function [z,weights,Dz] = nativeDifferentiationRule(self,zBounds,derivativeOrder)
            % Return one ordered grid for physical quadrature and differentiation.
            %
            % The returned points increase in physical coordinate. `weights`
            % integrate values in that order, and `Dz*values` differentiates
            % columns sampled in the same order. This keeps the grid, weights,
            % and differentiation matrix as one inseparable numerical rule.
            %
            % - Topic: Solve EVPs
            % - Declaration: [z,weights,Dz] = nativeDifferentiationRule(solver,zBounds,derivativeOrder)
            % - Parameter zBounds: physical integration bounds
            % - Parameter derivativeOrder: physical derivative order; default `1`
            % - Returns z: increasing physical grid points
            % - Returns weights: physical-coordinate quadrature weights
            % - Returns Dz: matrix mapping values on `z` to physical derivatives on `z`
            arguments
                self IMSolver
                zBounds (1,2) double
                derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 1
            end

            nativeZ = self.zNative(:);
            quadratureZ = self.innerProductGrid(zBounds);
            if length(quadratureZ) ~= length(nativeZ) || max(abs(sort(quadratureZ(:))-sort(nativeZ))) > 100*eps(max(1,max(abs(nativeZ))))
                error("IMSolver:UnsupportedDifferentiationBounds", ...
                    "nativeDifferentiationRule requires bounds whose quadrature uses the complete native solver grid.");
            end
            nativeWeights = self.innerProductWeights(nativeZ,zBounds);
            nativeDerivative = self.differentiateGridValues(eye(length(nativeZ)),derivativeOrder);
            [z,index] = sort(nativeZ(:));
            weights = nativeWeights(index);
            Dz = nativeDerivative(index,index);
        end

        function basisSet = solveEVP(self, evp, options)
            % Solve an EVP and return a basis set.
            %
            % If the assembled matrices produce no finite real eigenvalues,
            % `solveEVP` throws a matrix-level diagnostic before returning.
            % Before selection, the EVP may reject family-specific
            % numerical representations that are not finite physical modes.
            %
            % - Topic: Solve EVPs
            % - Declaration: basisSet = solveEVP(solver,evp,options)
            % - Parameter evp: canonical EVP descriptor
            % - Parameter options.nModes: number of modes to retain
            % - Returns basisSet: solved basis set
            arguments
                self IMSolver
                evp IMEigenvalueProblem
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 100
            end

            solver = self.configuredForEVP(evp);
            [A, B] = evp.assemble(solver);
            basisSet = solver.solveAssembledEVP(evp,A,B,options.nModes);
        end

        function zRoots = rootsOfNativeMode(self, nativeMode)
            % Return physical roots of one native mode.
            %
            % Concrete solvers may implement this developer hook when their
            % native representation supports accurate root finding. The base
            % implementation reports that mode-root grids are unavailable.
            %
            % - Topic: Developer topics
            % - Declaration: zRoots = rootsOfNativeMode(solver,nativeMode)
            % - Parameter nativeMode: one native mode column
            % - Returns zRoots: physical roots in the solver domain
            % - Developer: true
            arguments
                self IMSolver
                nativeMode (:,1) double {mustBeReal, mustBeFinite}
            end

            zRoots = zeros(0,1); %#ok<NASGU>
            error("IMSolver:UnsupportedModeRoots", ...
                "%s does not implement rootsOfNativeMode for %d native values.", class(self), length(nativeMode));
        end

        function basisSet = solveGeostrophicZeroAPVModes(self, problem)
            % Solve canonical geostrophic zero-APV boundary modes.
            %
            % The operator is factored once for each distinct requested
            % horizontal wavenumber. Every requested unit endpoint response is solved
            % in the same multiple-right-hand-side operation. No
            % generalized-energy coefficient or rotation enters this solve.
            % The regular first-order F/G equations permit N2 to vanish.
            %
            % - Topic: Solve geostrophic zero-APV modes
            % - Declaration: basisSet = solveGeostrophicZeroAPVModes(solver,problem)
            % - Parameter problem: geostrophic zero-APV problem
            % - Returns basisSet: canonical boundary-normalized basis
            arguments
                self IMSolver
                problem IMGeostrophicZeroAPVModes
            end

            solver = self.configuredForGeostrophicZeroAPVModes(problem);
            if isa(solver,"IMSolverFiniteDifference")
                solveSystems = @(matrix,rhs) solver.solveBoundaryValueSystems(matrix,rhs);
                basisSet = IMGeostrophicFiniteDifferenceTools.solveZeroAPV(problem,solver,solveSystems);
                return
            end
            z = solver.zNative;
            n = length(z);
            D0 = solver.physicalDerivativeMatrix(0);
            D1 = solver.physicalDerivativeMatrix(1);
            N2Values = problem.N2(z);
            N2Values = N2Values(:);
            if length(N2Values) ~= n || ~isreal(N2Values) || any(~isfinite(N2Values)) || any(N2Values < 0)
                error("IMGeostrophicZeroAPVModes:InvalidStratification", "N2 must return one finite real nonnegative value for each solver grid point.");
            end

            % The regular first-order pair remains valid in neutral layers.
            % Replace one row per equation with unit endpoint responses.
            % Keeping all rows of a first derivative is rank deficient when
            % its coupling coefficient vanishes throughout the domain.
            baseMatrix = [D1, diag(N2Values/problem.g)*D0; zeros(n), D1];
            surfaceIndex = solver.boundaryIndex("surface");
            bottomIndex = solver.boundaryIndex("bottom");
            surfaceRow = [zeros(1,n), D0(surfaceIndex,:)];
            if problem.surfaceBoundary == "freeSurface"
                surfaceRow(1:n) = -D0(surfaceIndex,:);
            end
            bottomRow = [zeros(1,n), D0(bottomIndex,:)];
            nEndpoints = numel(problem.endpoints);
            rightHandSides = zeros(2*n,nEndpoints);
            for iEndpoint = 1:nEndpoints
                switch problem.endpoints(iEndpoint)
                    case "surface"
                        rightHandSides(surfaceIndex,iEndpoint) = 1;
                    case "bottom"
                        rightHandSides(n+bottomIndex,iEndpoint) = 1;
                end
            end

            nK = numel(problem.k);
            nativeModes = zeros(n,nEndpoints,nK);
            nativeGModes = zeros(n,nEndpoints,nK);
            solvedWavenumbers = zeros(1,nK);
            solvedModes = cell(1,nK);
            nSolved = 0;
            degreeScale = solver.eigenproblemColumnScale(n);
            for iK = 1:nK
                solvedIndex = find(solvedWavenumbers(1:nSolved) == problem.k(iK),1);
                if isempty(solvedIndex)
                    a = problem.g*problem.k(iK)^2/problem.f0^2;
                    nSolved = nSolved + 1;
                    solvedWavenumbers(nSolved) = problem.k(iK);
                    if all(N2Values == 0)
                        % Only constant F and affine G remain. Solve for F
                        % and G at the bottom, avoiding a redundant pencil.
                        sigma = double(problem.surfaceBoundary == "freeSurface");
                        smallMatrix = [-(a*diff(problem.zDomain)+sigma), 1; 0, 1];
                        smallRHS = rightHandSides([surfaceIndex,n+bottomIndex],:);
                        response = solver.solveBoundaryValueSystems(smallMatrix,smallRHS);
                        FValues = ones(n,1)*response(1,:);
                        GValues = ones(n,1)*response(2,:)-a*(z(:)-problem.zDomain(1))*response(1,:);
                        solvedModes{nSolved} = [D0\FValues; D0\GValues];
                    else
                        matrix = baseMatrix;
                        matrix(n+(1:n),1:n) = a*D0;
                        matrix(surfaceIndex,:) = surfaceRow;
                        matrix(n+bottomIndex,:) = bottomRow;
                        % Balance the physical F/G amplitudes before the
                        % polynomial scaling and row equilibration.
                        FScale = sqrt(max(N2Values)/(problem.g*a));
                        columnScale = [FScale*degreeScale, degreeScale];
                        matrix = matrix.*columnScale;
                        rowScale = max(abs(matrix),[],2);
                        rowScale(rowScale == 0) = 1;
                        solvedModes{nSolved} = columnScale.'.*solver.solveBoundaryValueSystems(matrix./rowScale,rightHandSides./rowScale);
                    end
                    solvedIndex = nSolved;
                end
                nativeModes(:,:,iK) = solvedModes{solvedIndex}(1:n,:);
                nativeGModes(:,:,iK) = solvedModes{solvedIndex}(n+(1:n),:);
            end

            metadata = problem.metadata;
            metadata.solutionKind = "geostrophicZeroAPVModes";
            metadata.k = problem.k;
            metadata.endpoints = problem.endpoints;
            metadata.surfaceBoundary = problem.surfaceBoundary;
            metadata.modesPerWavenumber = nEndpoints;
            metadata.factorizations = nSolved;
            basisSet = IMGeostrophicZeroAPVModesBasis(problem=problem, solver=solver, nativeModes=nativeModes, nativeGModes=nativeGModes, metadata=metadata);
        end
    end

    methods (Hidden)
        function weights = innerProductWeights(self, z, zBounds)
            % Return the linear weights used by integrateInnerProduct.
            %
            % This generic implementation obtains the exact weights of the
            % configured solver's integration operator.  It lets callers
            % integrate many columns with one matrix multiplication while
            % preserving each concrete solver's existing quadrature rule.
            arguments
                self IMSolver
                z (:,1) double
                zBounds (1,2) double
            end

            nGrid = length(z);
            weights = zeros(nGrid,1);
            impulse = zeros(nGrid,1);
            for iGrid = 1:nGrid
                impulse(iGrid) = 1;
                weights(iGrid) = self.integrateInnerProduct(z,impulse,zBounds);
                impulse(iGrid) = 0;
            end
        end
    end

    methods (Access = protected)
        function [basisSet,costs] = solveAssembledEVP(self,evp,A,B,nModes,diagnostics)
            % Preserve the shared eigenpair filtering, labels and orientation.
            arguments (Input)
                self IMSolver
                evp IMEigenvalueProblem
                A (:,:) double
                B (:,:) double
                nModes (1,1) double {mustBeInteger,mustBePositive}
                diagnostics (1,1) struct = struct()
            end
            timer = tic;
            % Equilibrate the equations without changing the physical matrices
            % used for classification and diagnostics. Filter null modes in
            % this same equation scaling for scalar and bulk construction.
            [solveA,solveB,map] = evp.prepareSolvePencil(self,A,B);
            columnScale = zeros(1,size(solveA,2));
            offset = 0;
            for fieldSize = map.nativeFieldSizes
                columnScale(offset+(1:fieldSize)) = self.eigenproblemColumnScale(fieldSize);
                offset = offset+fieldSize;
            end
            physicalColumnScale = ones(1,size(solveA,2));
            if isfield(map,"nativeColumnScale"), physicalColumnScale = map.nativeColumnScale; end
            columnScale = columnScale.*physicalColumnScale;
            % Equilibrate after balancing physical field units. Otherwise
            % weak N2 can leave every interior row tiny beside boundary rows.
            rowScale = max(max(abs(solveA.*physicalColumnScale),[],2),max(abs(solveB.*physicalColumnScale),[],2));
            rowScale(rowScale == 0) = 1;
            solveA = solveA./rowScale;
            solveB = solveB./rowScale;
            [V, D] = eig(solveA.*columnScale, solveB.*columnScale);
            % Recover the original native representation before filtering,
            % normalization, or any physical evaluation: c = S*y.
            V = columnScale.'.*V;
            costs.eigensolveSeconds = toc(timer);
            timer = tic;
            eigenvalues = diag(D);
            [V,solveEigenvalues] = evp.prepareGeneralizedEigenpairs(self,solveA,solveB,V,eigenvalues);
            if isempty(solveEigenvalues)
                error("IMSolver:NoValidEigenvalues", "%s", ...
                    self.noValidEigenvalueMessage(evp, A, B, eigenvalues, false(size(eigenvalues))));
            end
            V = map.expandEigenvectors(V);
            if isfield(map,"restoreEigenvalues")
                solveEigenvalues = map.restoreEigenvalues(solveEigenvalues);
            end
            mapped = map.publicEigenvalues(V,solveEigenvalues);
            if isstruct(mapped)
                eigenvalues = mapped.values;
                publicValid = mapped.qualified;
            else
                eigenvalues = mapped;
                publicValid = true(size(eigenvalues));
            end
            selection = evp.selectModes(eigenvalues(:),nModes,self,A,diagnostics=diagnostics);
            if any(~publicValid(selection.sortIndex))
                first = find(~publicValid(selection.sortIndex),1);
                error("IMGeostrophicModes:UnresolvedMode", ...
                    "Requested mode %d failed eigenvalue qualification. Increase nEVP or request fewer modes; higher modes cannot replace it.",selection.modeNumber(first));
            end
            eigenvalues = eigenvalues(selection.sortIndex);
            eigenvalues(selection.modeNumber == 0) = 0;
            V = V(:,selection.sortIndex);
            if map.hasInverseDepth
                selection.modeSelectionDiagnostics.geostrophicInverseDepth = reshape(solveEigenvalues(selection.sortIndex),1,[]);
            end
            basisSet = evp.makeBasisSet(self, V, eigenvalues(:).', ...
                selection.modeNumber, selection.modeSelectionDiagnostics);
            basisSet = basisSet.orientModeSigns();
            costs.finalizationSeconds = toc(timer);
        end

        function values = solveBoundaryValueSystems(~, matrix, rightHandSides)
            % Solve one boundary-value matrix for multiple response columns.
            matrixFactorization = decomposition(matrix);
            values = matrixFactorization \ rightHandSides;
        end

        function scale = eigenproblemColumnScale(~,nColumns)
            % Grid-value discretizations retain their existing unknowns.
            scale = ones(1,nColumns);
        end
    end

    methods (Access = private)
        function message = noValidEigenvalueMessage(~, evp, A, B, eigenvalues, valid)
            normA = norm(A, "fro");
            normB = norm(B, "fro");
            rankA = rank(A);
            rankB = rank(B);
            finite = isfinite(real(eigenvalues)) & isfinite(imag(eigenvalues));
            threshold = 1e-12*max([normA, normB, 1]);
            hint = "";
            if normA <= threshold || normB <= threshold
                hint = " One assembled matrix is nearly zero; check EVP coefficient parameters and physical constants such as f0.";
            end
            message = sprintf('EVP "%s" produced no finite real eigenvalues. total eigenvalues=%d, finite eigenvalues=%d, valid-real eigenvalues=%d, norm(A,"fro")=%.3g, norm(B,"fro")=%.3g, rank(A)=%d, rank(B)=%d.%s', ...
                evp.name, length(eigenvalues), nnz(finite), nnz(valid), normA, normB, rankA, rankB, hint);
        end
    end

    methods (Abstract)
        solver = configuredForEVP(self, evp)

        % Return a solver configured for geostrophic zero-APV modes.
        %
        % Concrete solvers prepare their native grid, coordinate mapping,
        % and derivative matrices for the supplied zero-APV problem.
        %
        % - Topic: Solve geostrophic zero-APV modes
        % - Topic: Developer topics
        % - Declaration: solver = configuredForGeostrophicZeroAPVModes(solver,problem)
        % - Parameter problem: geostrophic zero-APV problem
        % - Returns solver: configured solver
        % - Developer: true
        solver = configuredForGeostrophicZeroAPVModes(self, problem)

        context = context(self)
        values = N2(self, z)
        D = physicalDerivativeMatrix(self, derivativeOrder)
        values = differentiateGridValues(self, values, derivativeOrder)
        index = boundaryIndex(self, location)
        values = evaluateNativeModes(self, nativeModes, z)
        values = evaluatePhysicalDerivative(self, nativeModes, z, derivativeOrder)
        nativeIntegral = integrateGridValuesFromSurface(self, gridValues)
        z = innerProductGrid(self, zBounds)
        value = integrateInnerProduct(self, z, integrand, zBounds)
    end
end

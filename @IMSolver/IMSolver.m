classdef (Abstract) IMSolver
    % Define the shared protocol for canonical EVP solvers.
    %
    % Concrete solvers own the grid, coordinate mapping, derivative
    % matrices, integration rule, and interpolation of native modes. The
    % base class owns the common generalized-eigenvalue workflow.
    %
    % - Topic: Solve EVPs
    % - Topic: Solve surface-geostrophic modes
    % - Topic: Developer topics
    % - Declaration: classdef (Abstract) IMSolver

    methods
        function basisSet = solveEVP(self, evp, options)
            % Solve an EVP and return a basis set.
            %
            % If the assembled matrices produce no finite real eigenvalues,
            % `solveEVP` throws a matrix-level diagnostic before returning.
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
            [V, D] = eig(A, B);
            eigenvalues = diag(D);
            valid = isfinite(real(eigenvalues)) & isfinite(imag(eigenvalues)) ...
                & abs(imag(eigenvalues)) < 1e-8*max(1,abs(real(eigenvalues)));
            if ~any(valid)
                error("IMSolver:NoValidEigenvalues", "%s", ...
                    self.noValidEigenvalueMessage(evp, A, B, eigenvalues, valid));
            end

            V = real(V(:,valid));
            eigenvalues = real(eigenvalues(valid));
            selection = evp.selectModes(eigenvalues(:), options.nModes, solver, A);
            eigenvalues = eigenvalues(selection.sortIndex);
            eigenvalues(selection.modeNumber == 0) = 0;
            V = V(:,selection.sortIndex);
            basisSet = evp.makeBasisSet(solver, V, eigenvalues(:).', ...
                selection.modeNumber, selection.modeSelectionDiagnostics);
            basisSet = basisSet.orientModeSigns();
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

        function basisSet = solveSurfaceGeostrophicModes(self, problem)
            % Solve projected surface-geostrophic boundary modes.
            %
            % `solveSurfaceGeostrophicModes` solves the raw zero-APV
            % endpoint modes stored by `IMSurfaceGeostrophicModes`, forms
            % the boundary-energy projection, and returns an
            % `IMSurfaceGeostrophicModesBasis` with `F`, `G`, and `h`.
            %
            % - Topic: Solve surface-geostrophic modes
            % - Declaration: basisSet = solveSurfaceGeostrophicModes(solver,problem)
            % - Parameter problem: surface-geostrophic boundary-mode problem
            % - Returns basisSet: solved surface-geostrophic basis
            arguments
                self IMSolver
                problem IMSurfaceGeostrophicModes
            end

            solver = self.configuredForSurfaceGeostrophicModes(problem);
            z = solver.zNative;
            n = length(z);
            D0 = solver.physicalDerivativeMatrix(0);
            D1 = solver.physicalDerivativeMatrix(1);
            D2 = solver.physicalDerivativeMatrix(2);
            N2Values = problem.N2(z);
            N2Values = N2Values(:);
            if length(N2Values) ~= n
                error("IMSurfaceGeostrophicModes:InvalidStratification", "N2 must return one value for each solver grid point.");
            end
            if any(~isfinite(N2Values)) || any(N2Values <= 0)
                error("IMSurfaceGeostrophicModes:InvalidStratification", "N2 must be finite and positive on the solver grid.");
            end

            pValues = problem.f0^2 ./ N2Values;
            pzValues = solver.differentiateGridValues(pValues, 1);
            % Assemble raw zero-APV modes in divergence form for conditioning.
            baseMatrix = diag(pValues)*D2 + diag(pzValues)*D1;
            surfaceIndex = solver.boundaryIndex("surface");
            bottomIndex = solver.boundaryIndex("bottom");

            N2Surface = N2Values(surfaceIndex);
            N2Bottom = N2Values(bottomIndex);
            surfaceRow = D1(surfaceIndex,:);
            if problem.surfaceAnomaly == "freeSurface"
                surfaceRow = surfaceRow + (N2Surface/problem.g)*D0(surfaceIndex,:);
            end
            bottomRow = D1(bottomIndex,:);

            includeSurface = isfinite(problem.g0);
            includeBottom = isfinite(problem.gd);
            rawRows = [];
            targetAnomalies = zeros(2,0);
            scaledBoundaryTargets = zeros(2,0);
            if includeSurface
                rawRows(end+1) = 1;
                targetAnomalies(:,end+1) = [problem.g0/N2Surface; 0];
                scaledBoundaryTargets(:,end+1) = [-problem.g0/problem.f0; 0];
            end
            if includeBottom
                rawRows(end+1) = 2;
                targetAnomalies(:,end+1) = [0; problem.gd/N2Bottom];
                scaledBoundaryTargets(:,end+1) = [0; -problem.gd/problem.f0];
            end
            nRawModes = size(targetAnomalies,2);
            nColumns = nRawModes*numel(problem.k);
            nativeModes = zeros(n, nColumns);
            kByMode = zeros(1, nColumns);
            h = zeros(1, nColumns);
            modeNumber = zeros(1, nColumns);
            energyEigenvalues = zeros(1, nColumns);
            mixingCoefficients = zeros(2, nColumns);
            columnIndex = 0;
            for iK = 1:numel(problem.k)
                matrix = baseMatrix - problem.k(iK)^2*D0;
                matrix(surfaceIndex,:) = surfaceRow;
                matrix(bottomIndex,:) = bottomRow;

                rawModes = zeros(n, nRawModes);
                for iRaw = 1:nRawModes
                    rhs = zeros(n, 1);
                    rhs(surfaceIndex) = scaledBoundaryTargets(1,iRaw);
                    rhs(bottomIndex) = scaledBoundaryTargets(2,iRaw);
                    rawModes(:,iRaw) = matrix \ rhs;
                end

                rawModeValues = D0*rawModes;
                surfaceValues = rawModeValues(surfaceIndex,:).';
                bottomValues = rawModeValues(bottomIndex,:).';
                etaSurface = targetAnomalies(1,:).';
                etaBottom = targetAnomalies(2,:).';
                energyMatrix = -problem.f0*(surfaceValues*etaSurface.') + problem.f0*(bottomValues*etaBottom.');
                if includeSurface
                    energyMatrix = energyMatrix + problem.g0*(etaSurface*etaSurface.');
                end
                if includeBottom
                    energyMatrix = energyMatrix + problem.gd*(etaBottom*etaBottom.');
                end
                energyMatrix = 0.5*(energyMatrix + energyMatrix.');

                [C, Gamma] = eig(energyMatrix);
                gamma = real(diag(Gamma));
                localH = 2*problem.k(iK)^2*gamma;
                [localH, sortIndex] = sort(localH, "descend");
                gamma = gamma(sortIndex);
                C = C(:,sortIndex);
                projectedModes = rawModes*C;
                for iMode = 1:nRawModes
                    projectedValues = D0*projectedModes(:,iMode);
                    [~, referenceIndex] = max(abs(projectedValues));
                    if projectedValues(referenceIndex) < 0
                        projectedModes(:,iMode) = -projectedModes(:,iMode);
                        C(:,iMode) = -C(:,iMode);
                    end

                    columnIndex = columnIndex + 1;
                    nativeModes(:,columnIndex) = projectedModes(:,iMode);
                    kByMode(columnIndex) = problem.k(iK);
                    h(columnIndex) = localH(iMode);
                    modeNumber(columnIndex) = iMode;
                    energyEigenvalues(columnIndex) = gamma(iMode);
                    mixingCoefficients(rawRows,columnIndex) = C(:,iMode);
                end
            end

            metadata = problem.metadata;
            metadata.solutionKind = "surfaceGeostrophicModes";
            metadata.k = kByMode;
            metadata.surfaceAnomaly = problem.surfaceAnomaly;
            metadata.modesPerWavenumber = nRawModes;
            basisSet = IMSurfaceGeostrophicModesBasis(problem=problem, solver=solver, nativeModes=nativeModes, k=kByMode, h=h, modeNumber=modeNumber, mixingCoefficients=mixingCoefficients, energyEigenvalues=energyEigenvalues, metadata=metadata);
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

        % Return a solver configured for surface-geostrophic modes.
        %
        % Concrete solvers prepare their native grid, coordinate mapping,
        % and derivative matrices for the supplied SQG problem.
        %
        % - Topic: Solve surface-geostrophic modes
        % - Topic: Developer topics
        % - Declaration: solver = configuredForSurfaceGeostrophicModes(solver,problem)
        % - Parameter problem: surface-geostrophic boundary-mode problem
        % - Returns solver: configured solver
        % - Developer: true
        solver = configuredForSurfaceGeostrophicModes(self, problem)

        context = context(self)
        values = N2(self, z)
        D = physicalDerivativeMatrix(self, derivativeOrder)
        values = differentiateGridValues(self, values, derivativeOrder)
        index = boundaryIndex(self, location)
        values = evaluateNativeModes(self, nativeModes, z)
        values = evaluatePhysicalDerivative(self, nativeModes, z, derivativeOrder)
        z = innerProductGrid(self, zBounds)
        value = integrateInnerProduct(self, z, integrand, zBounds)
    end
end

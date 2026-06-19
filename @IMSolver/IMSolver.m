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
            V = V(:,selection.sortIndex);
            basisSet = evp.makeBasisSet(solver, V, eigenvalues(:).', ...
                selection.modeNumber, selection.modeSelectionDiagnostics);
            basisSet = basisSet.orientModeSigns();
        end

        function basisSet = solveSurfaceGeostrophicModes(self, problem)
            % Solve surface-geostrophic boundary modes.
            %
            % `solveSurfaceGeostrophicModes` solves the boundary-value
            % problem stored by `IMSurfaceGeostrophicModes` and returns an
            % `IMSurfaceGeostrophicModesBasis`.
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

            N2z = solver.differentiateGridValues(N2Values, 1);
            baseMatrix = diag(N2Values)*D2 - diag(N2z)*D1;
            metricMatrix = diag(N2Values.*N2Values)*D0;
            surfaceIndex = solver.boundaryIndex("surface");
            bottomIndex = solver.boundaryIndex("bottom");
            activeIndex = surfaceIndex;
            if problem.boundary == "bottom"
                activeIndex = bottomIndex;
            end

            rhs = zeros(n, 1);
            nativeModes = zeros(n, numel(problem.k));
            for iK = 1:numel(problem.k)
                matrix = baseMatrix - (problem.k(iK)^2/problem.f0^2)*metricMatrix;
                matrix(surfaceIndex,:) = problem.f0*D1(surfaceIndex,:);
                matrix(bottomIndex,:) = problem.f0*D1(bottomIndex,:);
                rhs(:) = 0;
                rhs(activeIndex) = 1;
                nativeModes(:,iK) = matrix \ rhs;
            end

            metadata = problem.metadata;
            metadata.solutionKind = "surfaceGeostrophicModes";
            metadata.boundary = problem.boundary;
            metadata.k = problem.k;
            basisSet = IMSurfaceGeostrophicModesBasis(problem=problem, solver=solver, nativeModes=nativeModes, metadata=metadata);
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

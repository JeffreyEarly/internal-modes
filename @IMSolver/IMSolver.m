classdef (Abstract) IMSolver
    % Define the shared protocol for canonical EVP solvers.
    %
    % Concrete solvers own the grid, coordinate mapping, derivative
    % matrices, integration rule, and interpolation of native modes. The
    % base class owns the common generalized-eigenvalue workflow.
    %
    % - Topic: Solve EVPs
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
            h = evp.hFromEigenvalue(eigenvalues(:).');
            basisSet = evp.makeBasisSet(solver, V, eigenvalues(:).', h, ...
                selection.modeNumber, selection.index);
            basisSet = basisSet.orientModeSigns();
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
        context = context(self)
        values = N2(self, z)
        values = dzLogN2(self, z)
        D = physicalDerivativeMatrix(self, derivativeOrder)
        values = differentiateGridValues(self, values, derivativeOrder)
        index = boundaryIndex(self, location)
        values = evaluateNativeModes(self, nativeModes, z)
        values = evaluatePhysicalDerivative(self, nativeModes, z, derivativeOrder)
        z = innerProductGrid(self, zBounds)
        value = integrateInnerProduct(self, z, integrand, zBounds)
    end
end

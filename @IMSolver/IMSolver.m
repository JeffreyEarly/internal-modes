classdef (Abstract) IMSolver
    % Define the shared protocol for internal-mode solvers.
    %
    % `IMSolver` owns the solver-independent generalized EVP
    % workflow. Concrete subclasses provide the native grid, physical
    % derivative matrices, boundary rows, and native-mode evaluation.
    % Solvers own the numerical medium and discretization. EVPs own the
    % physical constants and combine them with solver context during
    % assembly.
    %
    % - Topic: Solve EVPs
    % - Topic: Developer topics
    % - Declaration: classdef (Abstract) IMSolver

    methods
        function basisSet = solveEVP(self, evp, options)
            % Solve an EVP and return a native-basis solution set.
            %
            % If the assembled matrices produce no finite real eigenvalues,
            % `solveEVP` throws an explanatory diagnostic before returning.
            %
            % - Topic: Solve EVPs
            % - Declaration: basisSet = solveEVP(solver,evp,options)
            % - Parameter evp: physical-coordinate EVP descriptor
            % - Parameter options.nModes: number of modes to retain
            % - Returns basisSet: solved native basis set
            arguments
                self IMSolver
                evp IMEigenvalueProblem
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 100
            end

            [A, B] = evp.assemble(self);
            [V, D] = eig(A, B);
            eigenvalues = diag(D);
            valid = isfinite(real(eigenvalues)) & isfinite(imag(eigenvalues)) & abs(imag(eigenvalues)) < 1e-8*max(1,abs(real(eigenvalues)));
            if ~any(valid)
                error("IMSolver:NoValidEigenvalues", "%s", ...
                    self.noValidEigenvalueMessage(evp, A, B, eigenvalues, valid));
            end
            V = real(V(:,valid));
            eigenvalues = real(eigenvalues(valid));
            selection = evp.selectModes(eigenvalues, options.nModes, evp.contextForSolver(self));
            eigenvalues = eigenvalues(selection.sortIndex);
            V = V(:,selection.sortIndex);
            modeNumber = selection.modeNumber;
            index = selection.index;
            h = evp.hFromEigenvalue(eigenvalues(:).');
            basisSet = IMBasisSet(solver=self, evp=evp, nativeModes=V, ...
                eigenvalues=eigenvalues(:).', h=h, modeNumber=modeNumber, index=index, ...
                zDomain=self.zDomain, N2Function=@(z) self.N2(z));
            basisSet = basisSet.orientModeSigns();
        end

        function [A, B] = applyEndpointLaw(self, A, B, endpointLaw, options)
            % Apply a resolved endpoint law to a matrix pair.
            %
            % Resolved endpoint laws replace the solver-native row associated
            % with their physical location.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: [A,B] = applyEndpointLaw(solver,A,B,endpointLaw,options)
            % - Parameter A: left EVP matrix
            % - Parameter B: right EVP matrix
            % - Parameter endpointLaw: resolved endpoint law
            % - Parameter options.context: framework coefficient context
            % - Returns A: left matrix with the endpoint row applied
            % - Returns B: right matrix with the endpoint row applied
            arguments
                self IMSolver
                A double
                B double
                endpointLaw IMBoundary
                options.context struct = struct()
            end

            if endpointLaw.location == ""
                error("IMSolver:UnplacedBoundaryLaw", ...
                    "Boundary law ""%s"" must be resolved before assembly.", endpointLaw.family);
            end

            index = self.boundaryIndex(endpointLaw.location);
            A(index,:) = endpointLaw.leftOperator.boundaryRow(self, endpointLaw.location, context=options.context);
            B(index,:) = endpointLaw.rightOperator.boundaryRow(self, endpointLaw.location, context=options.context);
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
        context = context(self)
        values = N2(self, z)
        values = dzLogN2(self, z)
        D = physicalDerivativeMatrix(self, derivativeOrder)
        index = boundaryIndex(self, location)
        values = evaluateNativeModes(self, nativeModes, z)
        values = evaluatePhysicalDerivative(self, nativeModes, z, derivativeOrder)
        z = innerProductGrid(self, zBounds)
        value = integrateInnerProduct(self, z, integrand, zBounds)
    end

end

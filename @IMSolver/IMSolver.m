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

        function [A, B] = applyBoundaryCondition(self, A, B, boundaryCondition, options)
            % Apply a placed boundary condition to a matrix pair.
            %
            % Active metadata-only boundary conditions do not replace matrix
            % rows. Standard placed conditions replace the solver-native row
            % associated with their physical location.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: [A,B] = applyBoundaryCondition(solver,A,B,boundaryCondition,options)
            % - Parameter A: left EVP matrix
            % - Parameter B: right EVP matrix
            % - Parameter boundaryCondition: placed boundary condition
            % - Parameter options.context: framework coefficient context
            % - Returns A: boundary-conditioned left matrix
            % - Returns B: boundary-conditioned right matrix
            arguments
                self IMSolver
                A double
                B double
                boundaryCondition IMBoundary
                options.context struct = struct()
            end

            if boundaryCondition.family == "active" || boundaryCondition.family == "partialDepthPE"
                return;
            end
            if boundaryCondition.location == ""
                error("IMSolver:UnplacedBoundaryCondition", ...
                    "Boundary condition ""%s"" must be placed before assembly.", boundaryCondition.family);
            end

            index = self.boundaryIndex(boundaryCondition.location);
            A(index,:) = boundaryCondition.leftOperator.boundaryRow(self, boundaryCondition.location, context=options.context);
            B(index,:) = boundaryCondition.rightOperator.boundaryRow(self, boundaryCondition.location, context=options.context);
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

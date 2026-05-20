classdef (Abstract) IMSolver
    % Define the shared protocol for v2 internal-mode solvers.
    %
    % `IMSolver` owns the solver-independent generalized EVP
    % workflow. Concrete subclasses provide the native grid, physical
    % derivative matrices, boundary rows, and native-mode evaluation.
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
            [eigenvalues, sortIndex] = self.sortEigenvalues(eigenvalues, evp.ordering);
            V = V(:,sortIndex);

            nRetain = min(options.nModes, length(eigenvalues));
            eigenvalues = eigenvalues(1:nRetain);
            V = V(:,1:nRetain);
            h = evp.hFromEigenvalue(eigenvalues(:).');
            index = evp.indexPolicy.classify(eigenvalues, self.context());
            basisSet = IMBasisSet(solver=self, evp=evp, nativeModes=V, ...
                eigenvalues=eigenvalues(:).', h=h, index=index, normalization=Normalization.kConstant, ...
                zDomain=self.zDomain, N2Function=@(z) self.N2(z), f0=self.f0, g=self.g);
            basisSet = basisSet.orientModeSigns();
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

    methods (Access = protected)
        function [lambdaSorted, sortIndex] = sortEigenvalues(~, eigenvalues, ordering)
            % Sort generalized-EVP eigenvalues according to an EVP policy.
            %
            % - Topic: Developer topics
            % - Declaration: [lambdaSorted,sortIndex] = sortEigenvalues(solver,eigenvalues,ordering)
            % - Parameter eigenvalues: eigenvalue vector
            % - Parameter ordering: ordering policy name
            % - Returns lambdaSorted: sorted eigenvalues
            % - Returns sortIndex: permutation indices
            % - Developer: true
            switch string(ordering)
                case "ascendingEigenvalue"
                    [lambdaSorted, sortIndex] = sort(eigenvalues, "ascend");
                case "descendingEigenvalue"
                    [lambdaSorted, sortIndex] = sort(eigenvalues, "descend");
                case "indexThenAscending"
                    [~, sortIndex] = sortrows([IMSolver.signWithZero(eigenvalues), abs(eigenvalues), eigenvalues]);
                    lambdaSorted = eigenvalues(sortIndex);
                otherwise
                    error("IMSolver:InvalidOrdering", ...
                        "Unknown EVP ordering ""%s"".", ordering);
            end
        end
    end

    methods (Static, Access = private)
        function signs = signWithZero(values)
            tolerance = 1e-10*max(1,max(abs(values)));
            signs = ones(size(values));
            signs(values < -tolerance) = -1;
            signs(abs(values) <= tolerance) = 0;
        end
    end
end

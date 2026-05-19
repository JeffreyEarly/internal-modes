classdef InternalModesSolverFiniteDifference
    % Solve physical-coordinate EVPs on a supplied finite-difference grid.
    %
    % `InternalModesSolverFiniteDifference` uses the user's physical `z`
    % grid as its native basis. It shares the v2 `solveEVP` interface with
    % the spectral solvers but evaluates modes by interpolation.
    %
    % ```matlab
    % solver = InternalModesSolverFiniteDifference(z=linspace(-1000,0,65).');
    % ```
    %
    % - Topic: Create solvers
    % - Topic: Solve EVPs
    % - Topic: Evaluate native modes
    % - Declaration: classdef InternalModesSolverFiniteDifference

    properties (SetAccess = private)
        % Number of native EVP values.
        %
        % - Topic: Inspect solvers
        nEVP

        % Physical vertical domain.
        %
        % - Topic: Inspect solvers
        zDomain

        % Native finite-difference grid.
        %
        % - Topic: Inspect solvers
        zNative

        % Native coordinate alias for API consistency.
        %
        % - Topic: Inspect solvers
        xNative

        % Native value matrix.
        %
        % - Topic: Developer topics
        % - Developer: true
        T

        % Native first-derivative matrix.
        %
        % - Topic: Developer topics
        % - Developer: true
        Tx

        % Native second-derivative matrix.
        %
        % - Topic: Developer topics
        % - Developer: true
        Txx

        % Coriolis parameter.
        %
        % - Topic: Inspect solvers
        f0

        % Gravitational acceleration.
        %
        % - Topic: Inspect solvers
        g
    end

    properties (Access = private)
        % Buoyancy frequency squared function.
        %
        % - Topic: Developer topics
        % - Developer: true
        N2Function
    end

    methods
        function self = InternalModesSolverFiniteDifference(options)
            % Create a finite-difference solver from a physical grid.
            %
            % - Topic: Create solvers
            % - Declaration: solver = InternalModesSolverFiniteDifference(options)
            % - Parameter options.z: finite-difference grid
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Returns solver: initialized finite-difference solver
            arguments
                options.z (:,1) double
                options.N2 function_handle = @(z) 1e-5*ones(size(z))
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
            end

            self.zNative = sort(options.z(:), "descend");
            self.xNative = self.zNative;
            self.zDomain = [min(self.zNative) max(self.zNative)];
            self.nEVP = length(self.zNative);
            self.N2Function = options.N2;
            self.f0 = options.f0;
            self.g = options.g;
            self.T = eye(self.nEVP);
            self.Tx = self.finiteDifferenceMatrix(1);
            self.Txx = self.finiteDifferenceMatrix(2);
        end

        function basisSet = solveEVP(self, evp, options)
            % Solve an EVP and return a native-basis solution set.
            %
            % - Topic: Solve EVPs
            % - Declaration: basisSet = solveEVP(solver,evp,options)
            % - Parameter evp: physical-coordinate EVP descriptor
            % - Parameter options.nModes: number of modes to retain
            % - Returns basisSet: solved native basis set
            arguments
                self InternalModesSolverFiniteDifference
                evp InternalModesEVP
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
            basisSet = InternalModesBasisSet(solver=self, evp=evp, nativeModes=V, ...
                eigenvalues=eigenvalues(:).', h=h, index=index, normalization=Normalization.kConstant);
        end

        function context = context(self)
            % Return coefficient functions used by physical operators.
            %
            % - Topic: Solve EVPs
            % - Developer: true
            % - Declaration: context = context(solver)
            % - Returns context: solver context structure
            context.N2 = @(z) self.N2(z);
            context.dzLogN2 = @(z) self.dzLogN2(z);
            context.f0 = self.f0;
            context.g = self.g;
            context.zDomain = self.zDomain;
            context.coordinateKind = "finiteDifference";
        end

        function values = N2(self, z)
            % Evaluate buoyancy frequency squared.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = N2(solver,z)
            % - Parameter z: physical coordinate
            % - Returns values: buoyancy frequency squared
            values = self.N2Function(z);
        end

        function values = dzLogN2(self, z)
            % Evaluate $$\partial_z\log N^2$$ numerically.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = dzLogN2(solver,z)
            % - Parameter z: physical coordinate
            % - Returns values: derivative values
            zAscending = sort(self.zNative);
            valuesAscending = gradient(log(self.N2(zAscending)), zAscending);
            values = interp1(zAscending, valuesAscending, z, "pchip");
        end

        function D = physicalDerivativeMatrix(self, derivativeOrder)
            % Return a native finite-difference matrix for a physical derivative.
            %
            % - Topic: Solve EVPs
            % - Developer: true
            % - Declaration: D = physicalDerivativeMatrix(solver,derivativeOrder)
            % - Parameter derivativeOrder: physical derivative order
            % - Returns D: derivative matrix
            switch derivativeOrder
                case 0
                    D = self.T;
                case 1
                    D = self.Tx;
                case 2
                    D = self.Txx;
                otherwise
                    error("InternalModesSolverFiniteDifference:UnsupportedDerivativeOrder", ...
                        "Derivative order %d is not supported.", derivativeOrder);
            end
        end

        function index = boundaryIndex(self, location)
            % Return the native row for a physical boundary.
            %
            % - Topic: Solve EVPs
            % - Developer: true
            % - Declaration: index = boundaryIndex(solver,location)
            % - Parameter location: `"surface"` or `"bottom"`
            % - Returns index: native row index
            switch string(location)
                case "surface"
                    index = 1;
                case "bottom"
                    index = self.nEVP;
                otherwise
                    error("InternalModesSolverFiniteDifference:InvalidBoundaryLocation", ...
                        "Boundary location must be ""surface"" or ""bottom"".");
            end
        end

        function values = evaluateNativeModes(self, nativeModes, z)
            % Interpolate native finite-difference modes to requested points.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = evaluateNativeModes(solver,nativeModes,z)
            % - Parameter nativeModes: native mode columns
            % - Parameter z: physical evaluation points
            % - Returns values: interpolated mode values
            values = interp1(flip(self.zNative), flip(nativeModes,1), z(:), "pchip");
        end

        function values = evaluatePhysicalDerivative(self, nativeModes, z, derivativeOrder)
            % Evaluate finite-difference physical derivatives at requested points.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = evaluatePhysicalDerivative(solver,nativeModes,z,derivativeOrder)
            % - Parameter nativeModes: native mode columns
            % - Parameter z: physical evaluation points
            % - Parameter derivativeOrder: physical derivative order
            % - Returns values: derivative values
            derivativeValues = self.physicalDerivativeMatrix(derivativeOrder)*nativeModes;
            values = interp1(flip(self.zNative), flip(derivativeValues,1), z(:), "pchip");
        end

        function gram = componentGramMatrix(self, basisSet, component, zBounds)
            % Integrate a component Gram matrix in physical coordinates.
            %
            % - Topic: Evaluate native modes
            % - Declaration: gram = componentGramMatrix(solver,basisSet,component,zBounds)
            % - Parameter basisSet: basis set to evaluate
            % - Parameter component: component name
            % - Parameter zBounds: physical integration bounds
            % - Returns gram: component Gram matrix
            z = linspace(min(zBounds), max(zBounds), max(256,4*self.nEVP)).';
            values = basisSet.evaluate(component, z);
            switch string(component)
                case "G"
                    weight = self.N2(z)/self.g;
                otherwise
                    weight = ones(size(z));
            end
            gram = zeros(size(values,2), size(values,2));
            for iMode = 1:size(values,2)
                for jMode = iMode:size(values,2)
                    value = trapz(z, weight(:).*values(:,iMode).*values(:,jMode));
                    gram(iMode,jMode) = value;
                    gram(jMode,iMode) = value;
                end
            end
        end

        function x = xOfZ(~, z)
            % Return the native coordinate for a physical finite-difference grid.
            %
            % - Topic: Evaluate native modes
            % - Declaration: x = xOfZ(solver,z)
            % - Parameter z: physical coordinate
            % - Returns x: native coordinate
            x = z;
        end

        function z = zOfX(~, x)
            % Return the physical coordinate for a native finite-difference grid.
            %
            % - Topic: Evaluate native modes
            % - Declaration: z = zOfX(solver,x)
            % - Parameter x: native coordinate
            % - Returns z: physical coordinate
            z = x;
        end
    end

    methods (Access = private)
        function D = finiteDifferenceMatrix(self, derivativeOrder)
            D = zeros(self.nEVP,self.nEVP);
            for iColumn = 1:self.nEVP
                e = zeros(self.nEVP,1);
                e(iColumn) = 1;
                values = e;
                for iDerivative = 1:derivativeOrder
                    values = gradient(values, self.zNative);
                end
                D(:,iColumn) = values;
            end
        end

        function [lambdaSorted, sortIndex] = sortEigenvalues(~, eigenvalues, ordering)
            switch string(ordering)
                case "ascendingEigenvalue"
                    [lambdaSorted, sortIndex] = sort(eigenvalues, "ascend");
                case "descendingEigenvalue"
                    [lambdaSorted, sortIndex] = sort(eigenvalues, "descend");
                case "indexThenAscending"
                    [~, sortIndex] = sortrows([localSignWithZero(eigenvalues), abs(eigenvalues), eigenvalues]);
                    lambdaSorted = eigenvalues(sortIndex);
                otherwise
                    error("InternalModesSolverFiniteDifference:InvalidOrdering", ...
                        "Unknown EVP ordering ""%s"".", ordering);
            end
        end
    end
end

function signs = localSignWithZero(values)
tolerance = 1e-10*max(1,max(abs(values)));
signs = ones(size(values));
signs(values < -tolerance) = -1;
signs(abs(values) <= tolerance) = 0;
end

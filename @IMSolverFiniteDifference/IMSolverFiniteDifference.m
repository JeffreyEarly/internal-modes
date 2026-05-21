classdef IMSolverFiniteDifference < IMSolver
    % Solve physical-coordinate EVPs on a supplied finite-difference grid.
    %
    % `IMSolverFiniteDifference` uses the user's physical `z`
    % grid as its native basis. It shares the v2 `solveEVP` interface with
    % the spectral solvers but evaluates modes by interpolation.
    %
    % ```matlab
    % solver = IMSolverFiniteDifference(z=linspace(-1000,0,65).');
    % ```
    %
    % - Topic: Create solvers
    % - Topic: Inspect solvers
    % - Topic: Solve EVPs
    % - Topic: Assemble EVPs
    % - Topic: Evaluate native modes
    % - Topic: Developer topics
    % - Declaration: classdef IMSolverFiniteDifference

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
    end

    properties (Access = private)
        % Buoyancy frequency squared function.
        %
        % - Topic: Developer topics
        % - Developer: true
        N2Function
    end

    methods
        function self = IMSolverFiniteDifference(options)
            % Create a finite-difference solver from a physical grid.
            %
            % - Topic: Create solvers
            % - Declaration: solver = IMSolverFiniteDifference(options)
            % - Parameter options.z: finite-difference grid
            % - Parameter options.N2: buoyancy frequency squared function
            % - Returns solver: initialized finite-difference solver
            arguments
                options.z (:,1) double
                options.N2 function_handle = @(z) 1e-5*ones(size(z))
            end

            self.zNative = sort(options.z(:), "descend");
            self.xNative = self.zNative;
            self.zDomain = [min(self.zNative) max(self.zNative)];
            self.nEVP = length(self.zNative);
            self.N2Function = options.N2;
            self.T = eye(self.nEVP);
            self.Tx = self.finiteDifferenceMatrix(1);
            self.Txx = self.finiteDifferenceMatrix(2);
        end

        function context = context(self)
            % Return the framework coefficient context.
            %
            % Solvers provide medium and discretization fields. EVPs add
            % physical constants such as `ctx.g` and `ctx.f0`.
            %
            % - Topic: Solve EVPs
            % - Developer: true
            % - Declaration: context = context(solver)
            % - Returns context: framework coefficient context
            context.N2 = @(z) self.N2(z);
            context.dzLogN2 = @(z) self.dzLogN2(z);
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
                    error("IMSolverFiniteDifference:UnsupportedDerivativeOrder", ...
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
                    error("IMSolverFiniteDifference:InvalidBoundaryLocation", ...
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

        function z = innerProductGrid(self, zBounds)
            % Return a bounded physical grid for finite-difference inner products.
            %
            % - Topic: Evaluate native modes
            % - Developer: true
            % - Declaration: z = innerProductGrid(solver,zBounds)
            % - Parameter zBounds: physical integration bounds
            % - Returns z: ascending physical integration grid
            arguments
                self IMSolverFiniteDifference
                zBounds (1,2) double
            end

            zBounds = sort(zBounds);
            zAscending = sort(self.zNative);
            interior = zAscending > zBounds(1) & zAscending < zBounds(2);
            z = unique([zBounds(1); zAscending(interior); zBounds(2)]);
        end

        function value = integrateInnerProduct(self, z, integrand, zBounds)
            % Integrate finite-difference inner-product values with trapezoids.
            %
            % - Topic: Evaluate native modes
            % - Developer: true
            % - Declaration: value = integrateInnerProduct(solver,z,integrand,zBounds)
            % - Parameter z: physical integration grid
            % - Parameter integrand: physical-coordinate integrand values
            % - Parameter zBounds: physical integration bounds
            % - Returns value: definite integral over `zBounds`
            arguments
                self IMSolverFiniteDifference
                z (:,1) double
                integrand (:,1) double
                zBounds (1,2) double
            end

            [zSorted, sortIndex] = sort(z(:));
            value = trapz(zSorted, integrand(sortIndex));
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

    end
end

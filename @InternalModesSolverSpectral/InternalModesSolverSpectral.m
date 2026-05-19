classdef InternalModesSolverSpectral
    % Solve physical-coordinate EVPs with a Chebyshev spectral discretization.
    %
    % `InternalModesSolverSpectral` owns the numerical coordinate, the
    % Chebyshev grid, and the physical-coordinate pullback rules. It does
    % not own modal normalization.
    %
    % ```matlab
    % solver = InternalModesSolverSpectral(N2=@(z) 1e-5*ones(size(z)), zDomain=[-1000 0], nEVP=64);
    % basisSet = solver.solveEVP(InternalModesEVP.hydrostaticGModes(k=1e-4));
    % ```
    %
    % - Topic: Create solvers
    % - Topic: Solve EVPs
    % - Topic: Evaluate native modes
    % - Topic: Developer topics
    % - Declaration: classdef InternalModesSolverSpectral

    properties (SetAccess = private)
        % Number of native EVP coefficients.
        %
        % - Topic: Inspect solvers
        nEVP

        % Physical vertical domain.
        %
        % - Topic: Inspect solvers
        zDomain

        % Native coordinate kind.
        %
        % - Topic: Inspect solvers
        coordinateKind

        % Native Lobatto grid.
        %
        % - Topic: Inspect solvers
        xNative

        % Physical points corresponding to `xNative`.
        %
        % - Topic: Inspect solvers
        zNative

        % Chebyshev basis matrix on the native grid.
        %
        % - Topic: Developer topics
        % - Developer: true
        T

        % Native first-derivative basis matrix.
        %
        % - Topic: Developer topics
        % - Developer: true
        Tx

        % Native second-derivative basis matrix.
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

    properties (SetAccess = private)
        % Reference physical grid for coordinate interpolation.
        %
        % - Topic: Developer topics
        % - Developer: true
        zReference

        % Reference native coordinate grid.
        %
        % - Topic: Developer topics
        % - Developer: true
        xReference

        % Reference coordinate derivative $$dx/dz$$.
        %
        % - Topic: Developer topics
        % - Developer: true
        qReference

        % Reference physical derivative of $$dx/dz$$.
        %
        % - Topic: Developer topics
        % - Developer: true
        qzReference
    end

    properties (Access = protected)
        % Buoyancy frequency squared function.
        %
        % - Topic: Developer topics
        % - Developer: true
        N2Function
    end

    methods
        function self = InternalModesSolverSpectral(options)
            % Create a coordinate-aware spectral solver.
            %
            % - Topic: Create solvers
            % - Declaration: solver = InternalModesSolverSpectral(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.nEVP: number of EVP coefficients
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.coordinateKind: native coordinate kind
            % - Returns solver: initialized spectral solver
            arguments
                options.N2 function_handle = @(z) 1e-5*ones(size(z))
                options.zDomain (1,2) double = [-1 0]
                options.nEVP (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.nEVP, 4)} = 64
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.coordinateKind {mustBeTextScalar} = "z"
            end

            self.N2Function = options.N2;
            self.zDomain = sort(options.zDomain);
            self.nEVP = options.nEVP;
            self.f0 = options.f0;
            self.g = options.g;
            self.coordinateKind = string(options.coordinateKind);
            self = self.setupCoordinate();
            self = self.setupNativeGrid();
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
                self InternalModesSolverSpectral
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
            % - Topic: Assemble EVPs
            % - Developer: true
            % - Declaration: context = context(solver)
            % - Returns context: solver context structure
            context.N2 = @(z) self.N2(z);
            context.dzLogN2 = @(z) self.dzLogN2(z);
            context.f0 = self.f0;
            context.g = self.g;
            context.zDomain = self.zDomain;
            context.coordinateKind = self.coordinateKind;
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
            N2Reference = self.N2(self.zReference);
            dzLogN2Reference = gradient(log(N2Reference), self.zReference);
            values = interp1(self.zReference, dzLogN2Reference, z, "pchip");
        end

        function D = physicalDerivativeMatrix(self, derivativeOrder)
            % Return a native matrix for a physical derivative.
            %
            % - Topic: Assemble EVPs
            % - Developer: true
            % - Declaration: D = physicalDerivativeMatrix(solver,derivativeOrder)
            % - Parameter derivativeOrder: physical derivative order
            % - Returns D: matrix mapping coefficients to derivative values
            q = self.qAtZ(self.zNative);
            qz = self.qzAtZ(self.zNative);
            switch derivativeOrder
                case 0
                    D = self.T;
                case 1
                    D = diag(q)*self.Tx;
                case 2
                    D = diag(q.*q)*self.Txx + diag(qz)*self.Tx;
                otherwise
                    error("InternalModesSolverSpectral:UnsupportedDerivativeOrder", ...
                        "Derivative order %d is not supported.", derivativeOrder);
            end
        end

        function index = boundaryIndex(self, location)
            % Return the native row for a physical boundary.
            %
            % - Topic: Assemble EVPs
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
                    error("InternalModesSolverSpectral:InvalidBoundaryLocation", ...
                        "Boundary location must be ""surface"" or ""bottom"".");
            end
        end

        function values = evaluateNativeModes(self, nativeModes, z)
            % Evaluate native Chebyshev coefficient columns at physical points.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = evaluateNativeModes(solver,nativeModes,z)
            % - Parameter nativeModes: Chebyshev coefficient columns
            % - Parameter z: physical evaluation points
            % - Returns values: mode values at `z`
            x = self.xOfZ(z(:));
            transform = InternalModesSpectral.ChebyshevTransformForGrid(self.xNative, x);
            values = transform(nativeModes);
        end

        function values = evaluatePhysicalDerivative(self, nativeModes, z, derivativeOrder)
            % Evaluate physical derivatives of native modes.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = evaluatePhysicalDerivative(solver,nativeModes,z,derivativeOrder)
            % - Parameter nativeModes: Chebyshev coefficient columns
            % - Parameter z: physical evaluation points
            % - Parameter derivativeOrder: physical derivative order
            % - Returns values: derivative values at `z`
            z = z(:);
            x = self.xOfZ(z);
            [~, TxOut, TxxOut] = InternalModesSpectral.ChebyshevPolynomialsOnGrid(x, self.nEVP);
            q = self.qAtZ(z);
            qz = self.qzAtZ(z);
            switch derivativeOrder
                case 0
                    TOut = InternalModesSpectral.ChebyshevPolynomialsOnGrid(x, self.nEVP);
                    values = TOut*nativeModes;
                case 1
                    values = diag(q)*TxOut*nativeModes;
                case 2
                    values = (diag(q.*q)*TxxOut + diag(qz)*TxOut)*nativeModes;
                otherwise
                    error("InternalModesSolverSpectral:UnsupportedDerivativeOrder", ...
                        "Derivative order %d is not supported.", derivativeOrder);
            end
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
            nGrid = max(256, 4*self.nEVP);
            z = linspace(min(zBounds), max(zBounds), nGrid).';
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

        function x = xOfZ(self, z)
            % Map physical coordinate to native coordinate.
            %
            % - Topic: Evaluate native modes
            % - Declaration: x = xOfZ(solver,z)
            % - Parameter z: physical coordinate
            % - Returns x: native coordinate
            x = interp1(self.zReference, self.xReference, z, "pchip");
        end

        function z = zOfX(self, x)
            % Map native coordinate to physical coordinate.
            %
            % - Topic: Evaluate native modes
            % - Declaration: z = zOfX(solver,x)
            % - Parameter x: native coordinate
            % - Returns z: physical coordinate
            z = interp1(self.xReference, self.zReference, x, "pchip");
        end
    end

    methods (Access = protected)
        function self = setupCoordinate(self)
            nReference = max(2001, 20*self.nEVP);
            self.zReference = linspace(self.zDomain(1), self.zDomain(2), nReference).';
            self.qReference = self.coordinateDerivative(self.zReference);
            if any(self.qReference <= 0)
                error("InternalModesSolverSpectral:InvalidCoordinate", ...
                    "The native coordinate derivative dx/dz must be positive.");
            end
            self.xReference = cumtrapz(self.zReference, self.qReference);
            self.qzReference = gradient(self.qReference, self.zReference);
        end

        function self = setupNativeGrid(self)
            xMin = min(self.xReference);
            xMax = max(self.xReference);
            self.xNative = ((xMax - xMin)/2)*(cos(((0:self.nEVP-1).')*pi/(self.nEVP-1)) + 1) + xMin;
            self.zNative = self.zOfX(self.xNative);
            self.zNative(1) = self.zDomain(2);
            self.zNative(end) = self.zDomain(1);
            [self.T, self.Tx, self.Txx] = InternalModesSpectral.ChebyshevPolynomialsOnGrid(self.xNative, self.nEVP);
        end

        function q = coordinateDerivative(self, z)
            switch self.coordinateKind
                case "z"
                    q = ones(size(z));
                case "wkb"
                    q = sqrt(self.N2(z));
                case "density"
                    q = self.N2(z);
                otherwise
                    error("InternalModesSolverSpectral:InvalidCoordinateKind", ...
                        "Unknown coordinate kind ""%s"".", self.coordinateKind);
            end
        end

        function q = qAtZ(self, z)
            q = interp1(self.zReference, self.qReference, z, "pchip");
        end

        function qz = qzAtZ(self, z)
            qz = interp1(self.zReference, self.qzReference, z, "pchip");
        end

        function [lambdaSorted, sortIndex] = sortEigenvalues(~, eigenvalues, ordering)
            switch string(ordering)
                case "ascendingEigenvalue"
                    [lambdaSorted, sortIndex] = sort(eigenvalues, "ascend");
                case "descendingEigenvalue"
                    [lambdaSorted, sortIndex] = sort(eigenvalues, "descend");
                case "indexThenAscending"
                    [~, sortIndex] = sortrows([signWithZero(eigenvalues), abs(eigenvalues), eigenvalues]);
                    lambdaSorted = eigenvalues(sortIndex);
                otherwise
                    error("InternalModesSolverSpectral:InvalidOrdering", ...
                        "Unknown EVP ordering ""%s"".", ordering);
            end
        end
    end
end

function signs = signWithZero(values)
tolerance = 1e-10*max(1,max(abs(values)));
signs = ones(size(values));
signs(values < -tolerance) = -1;
signs(abs(values) <= tolerance) = 0;
end

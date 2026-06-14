classdef IMSolverSpectral < IMSolver
    % Solve physical-coordinate EVPs with a Chebyshev spectral discretization.
    %
    % `IMSolverSpectral` owns the numerical coordinate choice, Chebyshev
    % resolution, derivative matrices, and physical-coordinate pullback
    % rules. It is configured against an EVP before solving.
    %
    % ```matlab
    % evp = IMInternalModes.waveModesAtWavenumber(N2=@(z) 1e-5*ones(size(z)), zDomain=[-1000 0], k=1e-4);
    % solver = IMSolverSpectral(nEVP=64);
    % basisSet = solver.solveEVP(evp);
    % ```
    %
    % - Topic: Create solvers
    % - Topic: Inspect solvers
    % - Topic: Solve EVPs
    % - Topic: Assemble EVPs
    % - Topic: Evaluate native modes
    % - Topic: Developer topics
    % - Declaration: classdef IMSolverSpectral

    properties (SetAccess = private)
        % Number of native EVP coefficients.
        %
        % - Topic: Inspect solvers
        nEVP

        % Physical vertical domain.
        %
        % - Topic: Inspect solvers
        zDomain = [NaN NaN]

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
        % EVP-provided buoyancy frequency squared function.
        %
        % - Topic: Developer topics
        % - Developer: true
        N2Function = []
    end

    methods
        function self = IMSolverSpectral(options)
            % Create a coordinate-aware spectral solver.
            %
            % - Topic: Create solvers
            % - Declaration: solver = IMSolverSpectral(options)
            % - Parameter options.nEVP: number of EVP coefficients
            % - Parameter options.coordinateKind: native coordinate kind
            % - Returns solver: initialized spectral solver
            arguments
                options.nEVP (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.nEVP, 4)} = 64
                options.coordinateKind {mustBeTextScalar} = "z"
            end

            self.nEVP = options.nEVP;
            self.coordinateKind = string(options.coordinateKind);
        end

        function solver = configuredForEVP(self, evp)
            % Return a spectral solver configured for an EVP.
            %
            % - Topic: Assemble EVPs
            % - Topic: Developer topics
            % - Declaration: solver = configuredForEVP(solver,evp)
            % - Parameter evp: canonical EVP descriptor
            % - Returns solver: solver with grid and coordinate matrices initialized
            % - Developer: true
            arguments
                self IMSolverSpectral
                evp IMEigenvalueProblem
            end

            solver = self;
            profile = evp.coordinateProfile(self.coordinateKind);
            solver.zDomain = evp.zDomain;
            solver.N2Function = profile.N2;
            solver = solver.setupCoordinate();
            solver = solver.setupNativeGrid();
        end

        function context = context(self)
            % Return the framework coefficient context.
            %
            % Solvers provide discretization fields. EVPs add the physical
            % domain, medium, and constants.
            %
            % - Topic: Assemble EVPs
            % - Developer: true
            % - Declaration: context = context(solver)
            % - Returns context: framework coefficient context
            context.coordinateKind = self.coordinateKind;
        end

        function values = N2(self, z)
            % Evaluate buoyancy frequency squared.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = N2(solver,z)
            % - Parameter z: physical coordinate
            % - Returns values: buoyancy frequency squared
            if isempty(self.N2Function)
                error("IMSolverSpectral:UnsupportedOperation", ...
                    "This spectral solver is not configured with an N2 function.");
            end
            values = self.N2Function(z);
        end

        function values = differentiateGridValues(self, values, derivativeOrder)
            % Differentiate values sampled on the native grid.
            %
            % - Topic: Assemble EVPs
            % - Developer: true
            % - Declaration: values = differentiateGridValues(solver,values,derivativeOrder)
            % - Parameter values: one value per native grid point
            % - Parameter derivativeOrder: physical derivative order
            % - Returns values: differentiated grid values
            values = values(:,:);
            if size(values,1) ~= self.nEVP
                error("IMSolverSpectral:InvalidGridValues", ...
                    "Grid values must have one row per native grid point.");
            end
            if derivativeOrder == 0
                return;
            end
            coefficients = self.T \ values;
            values = self.physicalDerivativeMatrix(derivativeOrder)*coefficients;
        end

        function values = dzLogN2(self, z)
            % Evaluate $$\partial_z\log N^2$$ numerically.
            %
            % - Topic: Evaluate native modes
            % - Declaration: values = dzLogN2(solver,z)
            % - Parameter z: physical coordinate
            % - Returns values: derivative values
            if isempty(self.N2Function)
                error("IMSolverSpectral:UnsupportedOperation", ...
                    "This spectral solver is not configured with an N2 function.");
            end
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
                    error("IMSolverSpectral:UnsupportedDerivativeOrder", ...
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
                    error("IMSolverSpectral:InvalidBoundaryLocation", ...
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
            z = z(:);
            if self.isNativePhysicalGrid(z)
                TOut = self.T;
            else
                x = self.clampNativeCoordinate(self.xOfZ(z));
                TOut = self.chebyshevPolynomialsAtNativePoints(x);
            end
            values = TOut*nativeModes;
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
            if self.isNativePhysicalGrid(z)
                TOut = self.T;
                TxOut = self.Tx;
                TxxOut = self.Txx;
            else
                x = self.clampNativeCoordinate(self.xOfZ(z));
                [TOut, TxOut, TxxOut] = self.chebyshevPolynomialsAtNativePoints(x);
            end
            q = self.qAtZ(z);
            qz = self.qzAtZ(z);
            switch derivativeOrder
                case 0
                    values = TOut*nativeModes;
                case 1
                    values = diag(q)*TxOut*nativeModes;
                case 2
                    values = (diag(q.*q)*TxxOut + diag(qz)*TxOut)*nativeModes;
                otherwise
                    error("IMSolverSpectral:UnsupportedDerivativeOrder", ...
                        "Derivative order %d is not supported.", derivativeOrder);
            end
        end

        function z = innerProductGrid(self, zBounds)
            % Return the native grid used for spectral inner products.
            %
            % - Topic: Evaluate native modes
            % - Developer: true
            % - Declaration: z = innerProductGrid(solver,zBounds)
            % - Parameter zBounds: physical integration bounds
            % - Returns z: physical points corresponding to the native grid
            arguments
                self IMSolverSpectral
                zBounds (1,2) double
            end

            z = self.zNative;
        end

        function value = integrateInnerProduct(self, z, integrand, zBounds)
            % Integrate inner-product values in the native Chebyshev coordinate.
            %
            % The supplied `integrand` is a physical-coordinate integrand
            % sampled at `z`. The method integrates
            % $$f(z(x))\,q^{-1}(z(x))$$ in the native coordinate $$x$$,
            % where $$q=dx/dz$$.
            %
            % - Topic: Evaluate native modes
            % - Developer: true
            % - Declaration: value = integrateInnerProduct(solver,z,integrand,zBounds)
            % - Parameter z: physical points from `innerProductGrid`
            % - Parameter integrand: physical-coordinate integrand values
            % - Parameter zBounds: physical integration bounds
            % - Returns value: definite integral over `zBounds`
            arguments
                self IMSolverSpectral
                z (:,1) double
                integrand (:,1) double
                zBounds (1,2) double
            end

            if length(z) ~= self.nEVP || max(abs(z(:) - self.zNative(:))) > self.gridTolerance()
                error("IMSolverSpectral:InvalidInnerProductGrid", ...
                    "Spectral inner products must be evaluated on the solver's native grid.");
            end
            if length(integrand) ~= self.nEVP
                error("IMSolverSpectral:InvalidIntegrandSize", ...
                    "The integrand must have one value per native grid point.");
            end

            q = self.qAtZ(self.zNative);
            integrandCheb = InternalModesSpectral.fct(integrand(:)./q(:));
            if self.boundsCoverDomain(zBounds)
                value = sum(self.chebyshevIntegrationWeights().*integrandCheb);
            else
                xBounds = self.xOfZ(sort(zBounds(:)));
                xBounds = min(max(xBounds, min(self.xNative)), max(self.xNative));
                value = InternalModesSpectral.IntegrateChebyshevVectorWithLimits( ...
                    integrandCheb, self.xNative, min(xBounds), max(xBounds));
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
                error("IMSolverSpectral:InvalidCoordinate", ...
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
            [self.T, self.Tx, self.Txx] = self.chebyshevPolynomialsAtNativePoints(self.xNative);
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
                    error("IMSolverSpectral:InvalidCoordinateKind", ...
                        "Unknown coordinate kind ""%s"".", self.coordinateKind);
            end
        end

        function q = qAtZ(self, z)
            q = self.coordinateDerivative(z);
        end

        function qz = qzAtZ(self, z)
            qz = interp1(self.zReference, self.qzReference, z, "pchip");
        end

        function x = clampNativeCoordinate(self, x)
            x = min(max(x, min(self.xNative)), max(self.xNative));
        end

        function value = isNativePhysicalGrid(self, z)
            value = length(z) == self.nEVP && max(abs(z(:) - self.zNative(:))) <= self.gridTolerance();
        end

        function [T, Tx, Txx] = chebyshevPolynomialsAtNativePoints(self, x)
            x = x(:);
            xMin = min(self.xNative);
            xMax = max(self.xNative);
            L = xMax - xMin;
            xNorm = (2/L)*(x - xMin) - 1;
            xNorm = min(max(xNorm, -1), 1);
            theta = acos(xNorm);

            T = zeros(length(x), self.nEVP);
            for iPoly = 0:(self.nEVP-1)
                T(:,iPoly+1) = cos(iPoly*theta);
            end

            if nargout < 2
                return
            end

            Tx = self.differentiateChebyshevBasis(T, L);
            if nargout < 3
                return
            end

            Txx = self.differentiateChebyshevBasis(Tx, L);
        end

        function Tx = differentiateChebyshevBasis(~, T, L)
            nPolys = size(T,2);
            Tx = zeros(size(T));
            Tx(:,2) = T(:,1);
            Tx(:,3) = 4*T(:,2);
            for j = 4:nPolys
                m = j - 1;
                Tx(:,j) = (m/(m-2))*Tx(:,j-2) + 2*m*T(:,j-1);
            end
            Tx = (2/L)*Tx;
        end

        function weights = chebyshevIntegrationWeights(self)
            np = (0:(self.nEVP-1)).';
            weights = -(1+(-1).^np)./(np.*np - 1);
            weights(2) = 0;
            weights = (max(self.xNative) - min(self.xNative))*weights/2;
        end

        function value = boundsCoverDomain(self, zBounds)
            tolerance = self.gridTolerance();
            zBounds = sort(zBounds);
            value = abs(zBounds(1) - self.zDomain(1)) <= tolerance && ...
                abs(zBounds(2) - self.zDomain(2)) <= tolerance;
        end

        function tolerance = gridTolerance(self)
            tolerance = 100*eps(max(1,max(abs(self.zDomain))));
        end

    end
end

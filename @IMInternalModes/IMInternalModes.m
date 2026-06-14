classdef IMInternalModes < IMEigenvalueProblem
    % Describe canonical EVPs with internal-mode interpretation.
    %
    % `IMInternalModes` translates standard `F` and `G` internal-mode
    % problems into the canonical scalar EVP. The solved scalar `u` is
    % either `F` or `G`; the other variable is recovered diagnostically by
    % the resulting `IMInternalModesBasis`.
    % Internal-mode EVPs own the stratification profile `N2` and physical
    % vertical domain used by solvers and basis sets.
    %
    % ```matlab
    % N2 = @(z) (5.2e-3)^2*exp(2*z/1300);
    % evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0]);
    % solver = IMSolverSpectral(nEVP=128,coordinateKind="wkb");
    % basisSet = solver.solveEVP(evp,nModes=4);
    % G = basisSet.G(linspace(-4000,0,200).');
    % ```
    %
    % - Topic: Create internal-mode EVPs
    % - Topic: Inspect internal-mode metadata
    % - Topic: Inspect internal-mode inner products
    % - Topic: Developer topics
    % - Declaration: classdef IMInternalModes < IMEigenvalueProblem

    properties (SetAccess = private)
        % Solved physical variable, `"F"` or `"G"`.
        %
        % `formulation` tells the canonical solver which variable is the
        % native unknown `u`. The complementary variable is evaluated
        % diagnostically by `IMInternalModesBasis`.
        %
        % - Topic: Inspect internal-mode metadata
        formulation = "G"

        % Buoyancy frequency squared function.
        %
        % `N2` has signature `values = N2(z)` and is owned by the EVP so
        % solvers can prepare `z`, `wkb`, or `density` coordinates from the
        % same continuous stratification.
        %
        % - Topic: Inspect internal-mode metadata
        N2

        % Coriolis parameter.
        %
        % `f0` is copied into `metadata.f0` and into coefficient contexts.
        %
        % - Topic: Inspect internal-mode metadata
        f0 = 0

        % Gravitational acceleration.
        %
        % `g` is copied into `metadata.g` and into coefficient contexts.
        %
        % - Topic: Inspect internal-mode metadata
        g = 9.81

        % Equivalent-depth conversion function.
        %
        % `hFromEigenvalue` maps retained eigenvalues to equivalent depths
        % for internal-mode basis sets. The handle has signature
        % `h = hFromEigenvalue(lambda)`, so
        % $$h_j=\texttt{hFromEigenvalue}(\lambda_j).$$
        %
        % - Topic: Inspect internal-mode metadata
        hFromEigenvalue = @(lambda) 1 ./ lambda
    end

    methods
        function self = IMInternalModes(options)
            % Create an internal-mode canonical EVP.
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes(options)
            % - Parameter options.name: short EVP name
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.formulation: solved variable, `"F"` or `"G"`
            % - Parameter options.p: canonical derivative-flux coefficient
            % - Parameter options.q: canonical left-side value coefficient
            % - Parameter options.r: canonical metric coefficient
            % - Parameter options.surfaceBoundary: surface endpoint condition
            % - Parameter options.bottomBoundary: bottom endpoint condition
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.hFromEigenvalue: equivalent-depth conversion
            % - Parameter options.metadata: additional scalar parameters
            % - Returns evp: internal-mode EVP descriptor
            arguments
                options.name {mustBeTextScalar} = "internalModes"
                options.zDomain (1,2) double
                options.N2 function_handle
                options.formulation {mustBeTextScalar} = "G"
                options.p = @(z,~) ones(size(z))
                options.q = @(z,~) zeros(size(z))
                options.r = @(z,~) ones(size(z))
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.hFromEigenvalue = @(lambda) 1 ./ lambda
                options.hasZeroMode (1,1) logical = false
                options.defaultNormalization = []
                options.normalizations struct = struct()
                options.metadata struct = struct()
            end

            formulation = IMInternalModes.validateFormulation(options.formulation);
            metadata = options.metadata;
            metadata.f0 = options.f0;
            metadata.g = options.g;
            metadata.formulation = formulation;

            self@IMEigenvalueProblem(name=options.name, p=options.p, q=options.q, r=options.r, ...
                zDomain=options.zDomain, surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                hasZeroMode=options.hasZeroMode, defaultNormalization=options.defaultNormalization, ...
                normalizations=options.normalizations, metadata=metadata);
            self.formulation = formulation;
            self.N2 = options.N2;
            self.f0 = options.f0;
            self.g = options.g;
            self.hFromEigenvalue = options.hFromEigenvalue;
        end

        function context = contextForSolver(self, solver)
            % Return the internal-mode coefficient context.
            %
            % The returned context extends the canonical context with `N2`,
            % `dzLogN2`, `f0`, `g`, and `formulation`.
            %
            % - Topic: Developer topics
            % - Declaration: context = contextForSolver(evp,solver)
            % - Parameter solver: canonical solver
            % - Returns context: coefficient context
            % - Developer: true
            context = contextForSolver@IMEigenvalueProblem(self, solver);
            context.N2 = @(z) self.N2(z);
            context.dzLogN2 = @(z) self.dzLogN2(z);
            context.f0 = self.f0;
            context.g = self.g;
            context.formulation = self.formulation;
        end

        function profile = coordinateProfile(self, coordinateKind)
            % Return internal-mode resources needed by a solver coordinate map.
            %
            % - Topic: Developer topics
            % - Declaration: profile = coordinateProfile(evp,coordinateKind)
            % - Parameter coordinateKind: `"z"`, `"wkb"`, or `"density"`
            % - Returns profile: struct with coordinate resources
            % - Developer: true
            coordinateKind = string(coordinateKind);
            if coordinateKind ~= "z" && coordinateKind ~= "wkb" && coordinateKind ~= "density"
                error("IMEigenvalueProblem:UnsupportedCoordinateKind", ...
                    "Internal-mode EVPs support coordinateKind=""z"", ""wkb"", or ""density"".");
            end
            profile.N2 = @(z) self.N2(z);
            profile.dzLogN2 = @(z) self.dzLogN2(z);
        end

        function values = dzLogN2(self, z)
            % Evaluate the vertical derivative of `log(N2)`.
            %
            % This derivative is used by coordinate mappings that need the
            % stratification slope. For more than one point it is computed
            % by finite differences on the supplied coordinate vector.
            %
            % - Topic: Inspect internal-mode metadata
            % - Declaration: values = dzLogN2(evp,z)
            % - Parameter z: physical coordinate
            % - Returns values: derivative values
            z = z(:);
            if length(z) == 1
                scale = max(1,abs(z));
                dz = sqrt(eps)*scale;
                values = (log(self.N2(z + dz)) - log(self.N2(z - dz)))/(2*dz);
            else
                values = gradient(log(self.N2(z)), z);
            end
        end

        function spec = innerProduct(self, variable)
            % Return the `F` or `G` inner-product recipe.
            %
            % For `G`, the interior weight is $$N^2/g$$. For `F`, the
            % interior weight is one. Endpoint metric terms are attached to
            % the solved formulation, because only that variable appears in
            % the canonical endpoint condition. The returned struct has
            % fields `variable`, `interiorWeight`, `surfaceWeights`,
            % `bottomWeights`, and `hasKnownBoundaryWeights`.
            %
            % - Topic: Inspect internal-mode inner products
            % - Declaration: spec = innerProduct(evp,variable)
            % - Parameter variable: `"F"` or `"G"`
            % - Returns spec: struct with interior and endpoint metric terms
            variable = IMInternalModes.validateFormulation(variable);
            spec.variable = variable;
            if variable == "G"
                spec.interiorWeight = @(z,ctx) ctx.N2(z)/ctx.g;
            else
                spec.interiorWeight = @(z,~) ones(size(z));
            end
            if variable == self.formulation
                spec.surfaceWeights = self.endpointWeights("surface");
                spec.bottomWeights = self.endpointWeights("bottom");
            else
                spec.surfaceWeights = struct("location", {}, "coefficient", {}, "c", {}, "d", {});
                spec.bottomWeights = struct("location", {}, "coefficient", {}, "c", {}, "d", {});
            end
            spec.hasKnownBoundaryWeights = true;
        end

        function basisSet = makeBasisSet(self, solver, nativeModes, eigenvalues, modeNumber, index)
            % Create an internal-mode basis set.
            %
            % - Topic: Developer topics
            % - Declaration: basisSet = makeBasisSet(evp,solver,nativeModes,eigenvalues,modeNumber,index)
            % - Returns basisSet: internal-mode basis set
            % - Developer: true
            h = self.hFromEigenvalue(eigenvalues);
            basisSet = IMInternalModesBasis(solver=solver, evp=self, nativeModes=nativeModes, ...
                eigenvalues=eigenvalues, h=h, modeNumber=modeNumber, index=index, ...
                zDomain=self.zDomain, N2Function=self.N2);
        end
    end

    methods (Static)
        function evp = hydrostaticGModes(options)
            % Create the hydrostatic `G` internal-mode EVP.
            %
            % The canonical scalar form is
            % $$-G''=\lambda N^2G/g.$$
            % The default normalization is `Normalization.geostrophic`, and
            % metadata includes `formulation`, `f0`, and `g`.
            %
            % ```matlab
            % evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0]);
            % solver = IMSolverSpectral(nEVP=128);
            % basisSet = solver.solveEVP(evp,nModes=4);
            % G = basisSet.G(z);
            % ```
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.hydrostaticGModes(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceBoundary: surface endpoint condition
            % - Parameter options.bottomBoundary: bottom endpoint condition
            % - Returns evp: hydrostatic `G` EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            normalizations = IMInternalModes.standardNormalizations();
            evp = IMInternalModes(name="hydrostaticGModes", formulation="G", N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ctx.N2(z)/ctx.g, f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                defaultNormalization=Normalization.geostrophic, normalizations=normalizations);
        end

        function evp = hydrostaticFModes(options)
            % Create the hydrostatic `F` internal-mode EVP.
            %
            % The canonical scalar form is
            % $$-\partial_z(N^{-2}F_z)=\lambda F/g.$$
            % It declares the barotropic zero mode.
            % Metadata includes `formulation` and `g`.
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.hydrostaticFModes(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceBoundary: surface endpoint condition
            % - Parameter options.bottomBoundary: bottom endpoint condition
            % - Returns evp: hydrostatic `F` EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double
                options.g (1,1) double {mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.neumann()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.neumann()
            end

            normalizations = IMInternalModes.standardNormalizations();
            evp = IMInternalModes(name="hydrostaticFModes", formulation="F", N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,ctx) 1./ctx.N2(z), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ones(size(z))/ctx.g, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                defaultNormalization=Normalization.geostrophic, normalizations=normalizations, ...
                hasZeroMode=true);
        end

        function evp = waveModesAtWavenumber(options)
            % Create the fixed-wavenumber wave-mode EVP.
            %
            % The canonical scalar form is
            % $$-G''+K^2G=\lambda(N^2-f_0^2)G/g.$$
            % The default normalization is `Normalization.kConstant`, and
            % metadata includes `k`, `formulation`, `f0`, and `g`.
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.waveModesAtWavenumber(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.k: horizontal wavenumber
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceBoundary: surface endpoint condition
            % - Parameter options.bottomBoundary: bottom endpoint condition
            % - Returns evp: fixed-wavenumber `G` EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double
                options.k (1,1) double {mustBeNonnegative}
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            k = options.k;
            metadata = struct("k", k);
            normalizations = IMInternalModes.standardNormalizations();
            evp = IMInternalModes(name="waveModesAtWavenumber", formulation="G", N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) k*k*ones(size(z)), ...
                r=@(z,ctx) (ctx.N2(z) - ctx.f0*ctx.f0)/ctx.g, ...
                f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                defaultNormalization=Normalization.kConstant, normalizations=normalizations, ...
                metadata=metadata);
        end

        function evp = waveModesAtFrequency(options)
            % Create the fixed-frequency wave-mode EVP.
            %
            % The canonical scalar form is
            % $$-G''=\lambda(N^2-\omega^2)G/g.$$
            % The default normalization is `Normalization.omegaConstant`,
            % and metadata includes `omega`, `formulation`, `f0`, and `g`.
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.waveModesAtFrequency(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.omega: wave frequency
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceBoundary: surface endpoint condition
            % - Parameter options.bottomBoundary: bottom endpoint condition
            % - Returns evp: fixed-frequency `G` EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double
                options.omega (1,1) double {mustBeNonnegative}
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            omega = options.omega;
            metadata = struct("omega", omega);
            normalizations = IMInternalModes.standardNormalizations();
            evp = IMInternalModes(name="waveModesAtFrequency", formulation="G", N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) (ctx.N2(z) - omega*omega)/ctx.g, ...
                f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                defaultNormalization=Normalization.omegaConstant, normalizations=normalizations, ...
                metadata=metadata);
        end
    end

    methods (Static, Access = private)
        function formulation = validateFormulation(formulation)
            formulation = string(formulation);
            if formulation ~= "F" && formulation ~= "G"
                error("IMInternalModes:InvalidFormulation", ...
                    "formulation must be ""F"" or ""G"".");
            end
        end

        function normalizations = standardNormalizations()
            normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor(basisSet.evp.formulation, iMode);
            normalizations.geostrophic = @(basisSet,iMode) basisSet.geostrophicNormFactor(iMode);
            normalizations.kConstant = @(basisSet,iMode) basisSet.innerProductNormFactor("G", iMode);
            normalizations.omegaConstant = @(basisSet,iMode) basisSet.innerProductNormFactor("F", iMode);
            normalizations.wMax = @(basisSet,iMode) basisSet.maxAbsFactor("G", iMode);
            normalizations.uMax = @(basisSet,iMode) basisSet.maxAbsFactor("F", iMode);
            normalizations.surfacePressure = @(basisSet,iMode) basisSet.surfacePressureNormFactor(iMode);
        end
    end
end

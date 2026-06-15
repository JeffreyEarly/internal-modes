classdef IMInternalModes < IMEigenvalueProblem
    % Describe canonical EVPs with internal-mode interpretation.
    %
    % `IMInternalModes` translates standard `F` and `G` internal-mode
    % problems into the canonical scalar EVP. The solved scalar `u` is
    % either `F` or `G`; the other variable is recovered diagnostically by
    % the relation handles `FfromGz` and `GfromFz` on the resulting
    % `IMInternalModesBasis`.
    % Internal-mode EVPs own the stratification profile `N2` and physical
    % vertical domain used by solvers and basis sets.
    %
    % ```matlab
    % N2 = @(z) (5.2e-3)^2*exp(2*z/1300);
    % evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=[-4000 0]);
    % solver = IMSolverSpectral(nEVP=128,coordinateKind="wkb");
    % basisSet = solver.solveEVP(evp,nModes=4);
    % basisSet.normalization = Normalization.geostrophic;
    % G = basisSet.G(linspace(-4000,0,200).');
    % ```
    %
    % - Topic: Create internal-mode EVPs
    % - Topic: Summarize internal-mode EVPs
    % - Topic: Inspect internal-mode configuration
    % - Topic: Inspect internal-mode inner products
    % - Topic: Developer topics
    % - Declaration: classdef IMInternalModes < IMEigenvalueProblem

    properties (SetAccess = private)
        % Solved physical variable, `"F"` or `"G"`.
        %
        % `formulation` tells the canonical solver which variable is the
        % native unknown `u`. The complementary variable is evaluated
        % diagnostically by `IMInternalModesBasis`. Coefficient handles can
        % read this value as `ctx.formulation`.
        %
        % - Topic: Inspect internal-mode configuration
        formulation = "G"

        % Buoyancy frequency squared function.
        %
        % `N2` has signature `values = N2(z)` and is owned by the EVP so
        % solvers can prepare `z`, `wkb`, or `density` coordinates from the
        % same continuous stratification.
        %
        % - Topic: Inspect internal-mode configuration
        N2

        % Coriolis parameter.
        %
        % `f0` is stored in `parameters.f0` and is available to coefficient
        % handles as `ctx.f0`.
        %
        % - Topic: Inspect internal-mode configuration
        f0 = 0

        % Gravitational acceleration.
        %
        % `g` is stored in `parameters.g` and is available to coefficient
        % handles as `ctx.g`.
        %
        % - Topic: Inspect internal-mode configuration
        g = 9.81

        % Equivalent-depth conversion function.
        %
        % `hFromEigenvalue` maps retained eigenvalues to equivalent depths
        % for internal-mode basis sets. The handle has signature
        % `h = hFromEigenvalue(lambda)`, so
        % $$h_j=\texttt{hFromEigenvalue}(\lambda_j).$$
        %
        % - Topic: Inspect internal-mode configuration
        hFromEigenvalue = @(lambda) 1 ./ lambda

        % Diagnostic relation from `G` derivative to `F`.
        %
        % `FfromGz` has signature `F = FfromGz(z,dGdz,h,ctx)`. The default
        % relation is
        % $$F_j(z)=h_j\frac{\partial G_j}{\partial z}(z).$$
        %
        % - Topic: Inspect internal-mode configuration
        FfromGz = @(z,dGdz,h,ctx) dGdz .* reshape(h, 1, [])

        % Diagnostic relation from `F` derivative to `G`.
        %
        % `GfromFz` has signature `G = GfromFz(z,dFdz,h,ctx)`. The default
        % relation is the hydrostatic inverse
        % $$G_j(z)=-\frac{g}{N^2(z)}
        % \frac{\partial F_j}{\partial z}(z).$$
        % Wave-mode factories install relation handles with the
        % appropriate wave correction factors.
        %
        % - Topic: Inspect internal-mode configuration
        GfromFz = @(z,dFdz,h,ctx) -(ctx.g./ctx.N2(z(:))).*dFdz
    end

    methods
        function self = IMInternalModes(options)
            % Create an internal-mode canonical EVP.
            %
            % The constructor copies `options.parameters`, then writes the
            % internal-mode fields `f0`, `g`, and `formulation` into the
            % parameter struct. These constructor-owned fields override
            % same-named entries in `options.parameters`.
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
            % - Parameter options.FfromGz: diagnostic relation from `G` derivative to `F`
            % - Parameter options.GfromFz: diagnostic relation from `F` derivative to `G`
            % - Parameter options.parameters: named coefficient parameters
            % - Returns evp: internal-mode EVP descriptor
            arguments
                options.name {mustBeTextScalar} = "internalModes"
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.N2 function_handle
                options.formulation {mustBeTextScalar, mustBeMember(options.formulation, ["F", "G"])} = "G"
                options.p = @(z,~) ones(size(z))
                options.q = @(z,~) zeros(size(z))
                options.r = @(z,~) ones(size(z))
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.hFromEigenvalue function_handle = @(lambda) 1 ./ lambda
                options.FfromGz function_handle = @(z,dGdz,h,ctx) dGdz .* reshape(h, 1, [])
                options.GfromFz function_handle = @(z,dFdz,h,ctx) -(ctx.g./ctx.N2(z(:))).*dFdz
                options.parameters struct = struct()
            end

            formulation = string(options.formulation);
            parameters = options.parameters;
            parameters.f0 = options.f0;
            parameters.g = options.g;
            parameters.formulation = formulation;

            self@IMEigenvalueProblem(name=options.name, p=options.p, q=options.q, r=options.r, ...
                zDomain=options.zDomain, surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                parameters=parameters);
            self.formulation = formulation;
            self.N2 = options.N2;
            self.f0 = options.f0;
            self.g = options.g;
            self.hFromEigenvalue = options.hFromEigenvalue;
            self.FfromGz = options.FfromGz;
            self.GfromFz = options.GfromFz;
        end

        summarize(self, solver)

        function context = contextForSolver(self, solver)
            % Return the internal-mode coefficient context.
            %
            % The returned context extends the canonical context with `N2`,
            % `f0`, `g`, and `formulation`.
            %
            % - Topic: Developer topics
            % - Declaration: context = contextForSolver(evp,solver)
            % - Parameter solver: canonical solver
            % - Returns context: coefficient context
            % - Developer: true
            arguments
                self IMInternalModes
                solver IMSolver
            end

            context = contextForSolver@IMEigenvalueProblem(self, solver);
            context.N2 = @(z) self.N2(z);
            context.f0 = self.f0;
            context.g = self.g;
            context.formulation = self.formulation;
        end

        function spec = innerProduct(self, variable)
            % Return the `F` or `G` inner-product recipe.
            %
            % For `G`, the interior weight is $$N^2/g$$. For `F`, the
            % interior weight is one. The returned struct has fields
            % `variable`, `interiorWeight`, `surfaceWeights`,
            % `bottomWeights`, `status`, and `reason`. `status` is
            % `"fixed"` or `"interiorOnly"` when a standalone Gram matrix is
            % available. It is `"unknown"`, `"mixed"`, or
            % `"eigenvalueDependent"` when the requested diagnostic
            % variable does not yet have an installed fixed inner-product
            % rule.
            %
            % - Topic: Inspect internal-mode inner products
            % - Declaration: spec = innerProduct(evp,variable)
            % - Parameter variable: optional variable name, `"F"` or `"G"`
            % - Returns spec: struct with interior and endpoint metric terms
            arguments
                self IMInternalModes
                variable {mustBeTextScalar, mustBeMember(variable, ["F", "G"])} = self.formulation
            end

            variable = string(variable);
            spec.variable = variable;
            if variable == "G"
                spec.interiorWeight = @(z,ctx) ctx.N2(z)/ctx.g;
            else
                spec.interiorWeight = @(z,~) ones(size(z));
            end
            if variable == self.formulation
                spec.surfaceWeights = self.endpointWeights("surface");
                spec.bottomWeights = self.endpointWeights("bottom");
                [spec.status, spec.reason] = self.solvedInnerProductStatus(spec.surfaceWeights, spec.bottomWeights);
            else
                spec.surfaceWeights = IMInternalModes.emptyEndpointWeights();
                spec.bottomWeights = IMInternalModes.emptyEndpointWeights();
                [spec.status, spec.reason] = self.diagnosticInnerProductStatus(variable);
            end
        end

        function basisSet = makeBasisSet(self, solver, nativeModes, eigenvalues, modeNumber, modeSelectionDiagnostics)
            % Create an internal-mode basis set.
            %
            % - Topic: Developer topics
            % - Declaration: basisSet = makeBasisSet(evp,solver,nativeModes,eigenvalues,modeNumber,modeSelectionDiagnostics)
            % - Returns basisSet: internal-mode basis set
            % - Developer: true
            arguments
                self IMInternalModes
                solver IMSolver
                nativeModes (:,:) double
                eigenvalues (1,:) double {mustBeReal, mustBeFinite}
                modeNumber (1,:) double {mustBeInteger}
                modeSelectionDiagnostics struct
            end

            h = self.hFromEigenvalue(eigenvalues);
            basisSet = IMInternalModesBasis(solver=solver, evp=self, nativeModes=nativeModes, ...
                eigenvalues=eigenvalues, h=h, modeNumber=modeNumber, modeSelectionDiagnostics=modeSelectionDiagnostics, ...
                zDomain=self.zDomain, N2=self.N2);
        end
    end

    methods (Static)
        function evp = hydrostaticGModes(options)
            % Create the hydrostatic `G` internal-mode EVP.
            %
            % The canonical scalar form is
            % $$-\frac{\partial^2 G}{\partial z^2}(z)
            % =\lambda\frac{N^2(z)}{g}G(z).$$
            % Solved hydrostatic basis sets install the `geostrophic`
            % normalization rule and use it by default. This factory sets
            % `parameters.formulation`, `parameters.f0`, and `parameters.g`.
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
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            evp = IMInternalModes(name="hydrostaticGModes", formulation="G", N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ctx.N2(z)/ctx.g, f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary);
        end

        function evp = hydrostaticFModes(options)
            % Create the hydrostatic `F` internal-mode EVP.
            %
            % The canonical scalar form is
            % $$-\frac{\partial}{\partial z}\left(N^{-2}(z)
            % \frac{\partial F}{\partial z}(z)\right)
            % =\lambda\frac{F(z)}{g}.$$
            % The barotropic zero mode is inferred from the canonical left
            % problem during mode selection.
            % This factory sets `parameters.formulation` and `parameters.g`;
            % `parameters.f0` is supplied by the internal-mode constructor
            % default.
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
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.neumann()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.neumann()
            end

            evp = IMInternalModes(name="hydrostaticFModes", formulation="F", N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,ctx) 1./ctx.N2(z), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ones(size(z))/ctx.g, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary);
        end

        function evp = waveModesAtWavenumber(options)
            % Create the fixed-wavenumber wave-mode EVP.
            %
            % The canonical scalar form is
            % $$-\frac{\partial^2 G}{\partial z^2}(z)+K^2G(z)
            % =\lambda\frac{N^2(z)-f_0^2}{g}G(z).$$
            % Solved fixed-wavenumber basis sets install the `kConstant`
            % normalization rule and use it by default.
            % This factory adds `parameters.k` and sets
            % `parameters.formulation`, `parameters.f0`, and `parameters.g`.
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
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.k (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            k = options.k;
            parameters = struct("k", k);
            GfromFz = @(z,dFdz,h,ctx) -(ctx.g./(ctx.N2(z(:)) - ctx.f0*ctx.f0 - ctx.g*reshape(h,1,[])*k*k)).*dFdz;
            evp = IMInternalModes(name="waveModesAtWavenumber", formulation="G", N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) k*k*ones(size(z)), ...
                r=@(z,ctx) (ctx.N2(z) - ctx.f0*ctx.f0)/ctx.g, ...
                f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                GfromFz=GfromFz, parameters=parameters);
        end

        function evp = waveModesAtFrequency(options)
            % Create the fixed-frequency wave-mode EVP.
            %
            % The canonical scalar form is
            % $$-\frac{\partial^2 G}{\partial z^2}(z)
            % =\lambda\frac{N^2(z)-\omega^2}{g}G(z).$$
            % Solved fixed-frequency basis sets use the generic `unity`
            % normalization by default. A fixed-frequency diagnostic `F`
            % inner-product normalization is deferred until the wave
            % diagnostic inner-product catalog is derived. This factory
            % adds `parameters.omega` and sets `parameters.formulation`,
            % `parameters.f0`, and `parameters.g`.
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
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.omega (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative}
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            omega = options.omega;
            parameters = struct("omega", omega);
            GfromFz = @(z,dFdz,h,ctx) -(ctx.g./(ctx.N2(z(:)) - omega*omega)).*dFdz;
            evp = IMInternalModes(name="waveModesAtFrequency", formulation="G", N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) (ctx.N2(z) - omega*omega)/ctx.g, ...
                f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                GfromFz=GfromFz, parameters=parameters);
        end
    end

    methods (Access = private)
        function [status, reason] = solvedInnerProductStatus(self, surfaceWeights, bottomWeights)
            surfaceActive = self.surfaceBoundary.isEigenvalueDependent();
            bottomActive = self.bottomBoundary.isEigenvalueDependent();
            if (surfaceActive && isempty(surfaceWeights)) || (bottomActive && isempty(bottomWeights))
                status = "unknown";
                reason = "At least one active endpoint condition has a degenerate or unavailable endpoint metric weight.";
            elseif isempty(surfaceWeights) && isempty(bottomWeights)
                status = "interiorOnly";
                reason = "The solved formulation has no endpoint metric terms, so the inner product is the interior integral only.";
            else
                status = "fixed";
                reason = "The solved formulation inner product follows the canonical scalar EVP endpoint weights.";
            end
        end

        function [status, reason] = diagnosticInnerProductStatus(self, variable)
            if string(self.name) == "hydrostaticGModes" && self.formulation == "G" && variable == "F" && self.hasDirichletEndpoints()
                status = "interiorOnly";
                reason = "For hydrostatic G modes with G=0 at both endpoints, the diagnostic F inner product is the interior F integral.";
            elseif string(self.name) == "hydrostaticFModes" && self.formulation == "F" && variable == "G" && self.hasNeumannEndpoints()
                status = "interiorOnly";
                reason = "For hydrostatic F modes with N^{-2} dF/dz=0 at both endpoints, the diagnostic G inner product is the interior N^2 G/g integral.";
            else
                status = "unknown";
                reason = "No fixed diagnostic inner-product catalog entry is installed for this EVP and boundary combination.";
            end
        end

        function tf = hasDirichletEndpoints(self)
            tf = IMInternalModes.isDirichletBoundary(self.surfaceBoundary) && IMInternalModes.isDirichletBoundary(self.bottomBoundary);
        end

        function tf = hasNeumannEndpoints(self)
            tf = IMInternalModes.isNeumannBoundary(self.surfaceBoundary) && IMInternalModes.isNeumannBoundary(self.bottomBoundary);
        end
    end

    methods (Static, Access = private)
        function weights = emptyEndpointWeights()
            weights = struct("location", {}, "coefficient", {}, "c", {}, "d", {});
        end

        function tf = isDirichletBoundary(boundary)
            tolerance = IMInternalModes.coefficientTolerance(boundary);
            tf = ~boundary.isEigenvalueDependent(tolerance) && abs(boundary.b) <= tolerance && abs(boundary.a) > tolerance;
        end

        function tf = isNeumannBoundary(boundary)
            tolerance = IMInternalModes.coefficientTolerance(boundary);
            tf = ~boundary.isEigenvalueDependent(tolerance) && abs(boundary.a) <= tolerance && abs(boundary.b) > tolerance;
        end

        function tolerance = coefficientTolerance(boundary)
            tolerance = 100*eps*max([1 abs(boundary.a) abs(boundary.b) abs(boundary.c) abs(boundary.d)]);
        end
    end
end

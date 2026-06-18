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

        % Physical mode-family declaration.
        %
        % `modeFamily` tells internal-mode utilities which physical
        % catalog and coupled normalization rules are meaningful for this
        % EVP. The default `"none"` installs only generic internal-mode
        % behavior. The `"hydrostatic"` family declares the hydrostatic
        % `F`/`G` family, enabling the generalized boundary-condition
        % catalog and the coupled `geostrophic` normalization convention.
        %
        % - Topic: Inspect internal-mode configuration
        modeFamily = "none"

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
            % - Parameter options.modeFamily: physical family, `"none"` or `"hydrostatic"`
            % - Parameter options.p: canonical derivative-flux coefficient
            % - Parameter options.q: canonical left-side value coefficient
            % - Parameter options.r: canonical metric coefficient
            % - Parameter options.surfaceBoundary: surface boundary condition
            % - Parameter options.bottomBoundary: bottom boundary condition
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
                options.modeFamily {mustBeTextScalar, mustBeMember(options.modeFamily, ["none", "hydrostatic"])} = "none"
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
            modeFamily = string(options.modeFamily);
            parameters = options.parameters;
            parameters.f0 = options.f0;
            parameters.g = options.g;
            parameters.formulation = formulation;

            self@IMEigenvalueProblem(name=options.name, p=options.p, q=options.q, r=options.r, ...
                zDomain=options.zDomain, surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                parameters=parameters);
            self.formulation = formulation;
            self.modeFamily = modeFamily;
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
            % `bottomWeights`, `endpointInnerProductTerms`,
            % `hasInnerProduct`, and `reason`. `hasInnerProduct` is true
            % when the variable has a known inner product. When it is false, Gram
            % matrices, spectra, and inner-product normalization for that
            % variable throw `IMInternalModesBasis:UnavailableInnerProduct`.
            % Diagnostic variables use the value-only hydrostatic endpoint
            % catalog only when `modeFamily` is `"hydrostatic"` and a
            % catalog row is known; other diagnostic inner products are
            % unavailable until a family catalog is added. Endpoint
            % inner-product terms from the catalog have the form
            % $$\alpha_\ell F_i(z_\ell)F_j(z_\ell)$$ or
            % $$\alpha_\ell G_i(z_\ell)G_j(z_\ell),$$
            % where $$z_\ell$$ is the bottom or surface endpoint. The
            % variable used in the endpoint term is stored as
            % `term.variable`, so a `G` inner product can contain an
            % endpoint term involving `F`, and conversely.
            %
            % - Topic: Inspect internal-mode inner products
            % - Declaration: spec = innerProduct(evp,variable)
            % - Parameter variable: optional variable name, `"F"` or `"G"`
            % - Returns spec: struct with interior and endpoint inner-product terms
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
            spec.endpointInnerProductTerms = IMHydrostaticInnerProductCatalog.emptyEndpointInnerProductTerms();
            if variable == self.formulation
                spec.surfaceWeights = self.endpointWeights("surface");
                spec.bottomWeights = self.endpointWeights("bottom");
                [spec.hasInnerProduct, spec.reason] = self.solvedInnerProductAvailability(spec.surfaceWeights, spec.bottomWeights);
            else
                spec.surfaceWeights = IMInternalModes.emptyEndpointWeights();
                spec.bottomWeights = IMInternalModes.emptyEndpointWeights();
                catalog = IMHydrostaticInnerProductCatalog.resolve(self, variable);
                spec.endpointInnerProductTerms = catalog.endpointInnerProductTerms;
                spec.hasInnerProduct = catalog.hasInnerProduct;
                spec.reason = catalog.reason;
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

            basisSet = IMInternalModesBasis(solver=solver, evp=self, nativeModes=nativeModes, eigenvalues=eigenvalues, modeNumber=modeNumber, modeSelectionDiagnostics=modeSelectionDiagnostics);
        end
    end

    methods (Static)
        function evp = hydrostaticGModes(options)
            % Create the hydrostatic `G` internal-mode EVP.
            %
            % This factory creates the hydrostatic `G`-form problem
            %
            % $$
            % -\frac{\partial^2 G_j}{\partial z^2}(z)
            % =
            % \lambda_j\frac{N^2(z)}{g}G_j(z),
            % \qquad \lambda_j=\frac{1}{h_j}.
            % $$
            %
            % At each endpoint $$z_\ell\in\{z_b,z_s\}$$, the corresponding
            % `IMBoundaryCondition(a=...,b=...,c=...,d=...)` is applied as
            %
            % $$
            % -\left[
            % a_\ell G_j(z_\ell)
            % -b_\ell\frac{\partial G_j}{\partial z}(z_\ell)
            % \right]
            % =
            % \lambda_j\left[
            % c_\ell G_j(z_\ell)
            % -d_\ell\frac{\partial G_j}{\partial z}(z_\ell)
            % \right].
            % $$
            %
            % The default surface and bottom boundary conditions are
            % `IMBoundaryCondition.dirichlet()`, giving rigid-lid and
            % rigid-bottom conditions
            %
            % $$
            % G_j(z_s)=0,\qquad G_j(z_b)=0.
            % $$
            %
            % Physical hydrostatic endpoint laws written in `F` and `G`
            % can be converted with `IMHydrostaticBoundaryCondition`
            % before they are passed to this factory:
            %
            % ```matlab
            % law = IMHydrostaticBoundaryCondition(a=A/g,b=1);
            % surfaceBoundary = law.canonicalBoundary(formulation="G",g=g);
            % evp = IMInternalModes.hydrostaticGModes(N2=N2,zDomain=zDomain,g=g,surfaceBoundary=surfaceBoundary);
            % ```
            %
            % After conversion, `innerProduct("F")` and
            % `innerProduct("G")` use the hydrostatic endpoint catalog to
            % report which bilinear forms are known.
            % Solved hydrostatic basis sets install the `geostrophic`
            % normalization rule and use it by default because they set
            % `modeFamily` to `"hydrostatic"`. This factory sets
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
            % - Parameter options.surfaceBoundary: surface boundary condition
            % - Parameter options.bottomBoundary: bottom boundary condition
            % - Returns evp: hydrostatic `G` EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            evp = IMInternalModes(name="hydrostaticGModes", formulation="G", modeFamily="hydrostatic", ...
                N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ctx.N2(z)/ctx.g, f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary);
        end

        function evp = hydrostaticFModes(options)
            % Create the hydrostatic `F` internal-mode EVP.
            %
            % This factory creates the hydrostatic `F`-form problem
            %
            % $$
            % -\frac{\partial}{\partial z}
            % \left(
            % \frac{1}{N^2(z)}
            % \frac{\partial F_j}{\partial z}(z)
            % \right)
            % =
            % \lambda_j\frac{F_j(z)}{g},
            % \qquad \lambda_j=\frac{1}{h_j}.
            % $$
            %
            % At each endpoint $$z_\ell\in\{z_b,z_s\}$$, the corresponding
            % `IMBoundaryCondition(a=...,b=...,c=...,d=...)` is applied as
            %
            % $$
            % -\left[
            % a_\ell F_j(z_\ell)
            % -b_\ell\frac{1}{N^2(z_\ell)}
            % \frac{\partial F_j}{\partial z}(z_\ell)
            % \right]
            % =
            % \lambda_j\left[
            % c_\ell F_j(z_\ell)
            % -d_\ell\frac{1}{N^2(z_\ell)}
            % \frac{\partial F_j}{\partial z}(z_\ell)
            % \right].
            % $$
            %
            % The default surface and bottom boundary conditions are
            % `IMBoundaryCondition.neumann()`, giving
            %
            % $$
            % \frac{1}{N^2(z_s)}
            % \frac{\partial F_j}{\partial z}(z_s)=0,\qquad
            % \frac{1}{N^2(z_b)}
            % \frac{\partial F_j}{\partial z}(z_b)=0.
            % $$
            %
            % Through the hydrostatic relation
            %
            % $$
            % G_j(z)=-\frac{g}{N^2(z)}
            % \frac{\partial F_j}{\partial z}(z),
            % $$
            %
            % these are the same rigid-lid and rigid-bottom conditions
            %
            % $$
            % G_j(z_s)=0,\qquad G_j(z_b)=0.
            % $$
            %
            % Physical hydrostatic endpoint laws written in `F` and `G`
            % can be converted with `IMHydrostaticBoundaryCondition`
            % before they are passed to this factory:
            %
            % ```matlab
            % law = IMHydrostaticBoundaryCondition(b=B,c=C);
            % surfaceBoundary = law.canonicalBoundary(formulation="F",g=g);
            % evp = IMInternalModes.hydrostaticFModes(N2=N2,zDomain=zDomain,g=g,surfaceBoundary=surfaceBoundary);
            % ```
            %
            % After conversion, `innerProduct("F")` and
            % `innerProduct("G")` use the hydrostatic endpoint catalog to
            % report which bilinear forms are known.
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
            % - Parameter options.surfaceBoundary: surface boundary condition
            % - Parameter options.bottomBoundary: bottom boundary condition
            % - Returns evp: hydrostatic `F` EVP
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.neumann()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.neumann()
            end

            evp = IMInternalModes(name="hydrostaticFModes", formulation="F", modeFamily="hydrostatic", ...
                N2=options.N2, zDomain=options.zDomain, ...
                p=@(z,ctx) 1./ctx.N2(z), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ones(size(z))/ctx.g, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary);
        end

        function evp = waveModesAtWavenumber(options)
            % Create the fixed-wavenumber wave-mode EVP.
            %
            % This factory creates the fixed-wavenumber `G`-form problem
            %
            % $$
            % -\frac{\partial^2 G_j}{\partial z^2}(z)
            % +k^2G_j(z)
            % =
            % \lambda_j\frac{N^2(z)-f_0^2}{g}G_j(z),
            % \qquad \lambda_j=\frac{1}{h_j}.
            % $$
            %
            % At each endpoint $$z_\ell\in\{z_b,z_s\}$$, the corresponding
            % `IMBoundaryCondition(a=...,b=...,c=...,d=...)` is applied as
            %
            % $$
            % -\left[
            % a_\ell G_j(z_\ell)
            % -b_\ell\frac{\partial G_j}{\partial z}(z_\ell)
            % \right]
            % =
            % \lambda_j\left[
            % c_\ell G_j(z_\ell)
            % -d_\ell\frac{\partial G_j}{\partial z}(z_\ell)
            % \right].
            % $$
            %
            % The default surface and bottom boundary conditions are
            % `IMBoundaryCondition.dirichlet()`, giving rigid endpoint
            % conditions
            %
            % $$
            % G_j(z_s)=0,\qquad G_j(z_b)=0.
            % $$
            %
            % A linear free-surface condition at the surface can be written
            % as
            %
            % $$
            % G_j(z_s)=h_j\frac{\partial G_j}{\partial z}(z_s),
            % \qquad \lambda_j=\frac{1}{h_j},
            % $$
            %
            % equivalently
            %
            % $$
            % \frac{\partial G_j}{\partial z}(z_s)
            % =
            % \lambda_j G_j(z_s).
            % $$
            %
            % In canonical boundary-condition coefficients this is
            % `IMBoundaryCondition(a=0,b=1,c=1,d=0)` at the surface.
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
            % - Parameter options.surfaceBoundary: surface boundary condition
            % - Parameter options.bottomBoundary: bottom boundary condition
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
            % This factory creates the fixed-frequency `G`-form problem
            %
            % $$
            % -\frac{\partial^2 G_j}{\partial z^2}(z)
            % =
            % \lambda_j\frac{N^2(z)-\omega^2}{g}G_j(z),
            % \qquad \lambda_j=\frac{1}{h_j}.
            % $$
            %
            % At each endpoint $$z_\ell\in\{z_b,z_s\}$$, the corresponding
            % `IMBoundaryCondition(a=...,b=...,c=...,d=...)` is applied as
            %
            % $$
            % -\left[
            % a_\ell G_j(z_\ell)
            % -b_\ell\frac{\partial G_j}{\partial z}(z_\ell)
            % \right]
            % =
            % \lambda_j\left[
            % c_\ell G_j(z_\ell)
            % -d_\ell\frac{\partial G_j}{\partial z}(z_\ell)
            % \right].
            % $$
            %
            % The default surface and bottom boundary conditions are
            % `IMBoundaryCondition.dirichlet()`, giving rigid endpoint
            % conditions
            %
            % $$
            % G_j(z_s)=0,\qquad G_j(z_b)=0.
            % $$
            %
            % A linear free-surface condition at the surface can be written
            % as
            %
            % $$
            % G_j(z_s)=h_j\frac{\partial G_j}{\partial z}(z_s),
            % \qquad \lambda_j=\frac{1}{h_j},
            % $$
            %
            % equivalently
            %
            % $$
            % \frac{\partial G_j}{\partial z}(z_s)
            % =
            % \lambda_j G_j(z_s).
            % $$
            %
            % In canonical boundary-condition coefficients this is
            % `IMBoundaryCondition(a=0,b=1,c=1,d=0)` at the surface.
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
            % - Parameter options.surfaceBoundary: surface boundary condition
            % - Parameter options.bottomBoundary: bottom boundary condition
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
        function [hasInnerProduct, reason] = solvedInnerProductAvailability(self, surfaceWeights, bottomWeights)
            surfaceActive = self.surfaceBoundary.isEigenvalueDependent();
            bottomActive = self.bottomBoundary.isEigenvalueDependent();
            if (surfaceActive && isempty(surfaceWeights)) || (bottomActive && isempty(bottomWeights))
                hasInnerProduct = false;
                reason = "At least one active boundary condition has a degenerate or unavailable endpoint metric weight.";
            elseif isempty(surfaceWeights) && isempty(bottomWeights)
                hasInnerProduct = true;
                reason = "The solved formulation has no endpoint metric terms, so the inner product is the interior integral only.";
            else
                hasInnerProduct = true;
                reason = "The solved formulation inner product follows the canonical scalar EVP endpoint weights.";
            end
        end
    end

    methods (Static, Access = private)
        function weights = emptyEndpointWeights()
            weights = struct("location", {}, "coefficient", {}, "c", {}, "d", {});
        end

    end
end

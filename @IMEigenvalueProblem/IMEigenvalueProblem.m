classdef IMEigenvalueProblem
    % Describe a vertical-mode generalized eigenvalue problem.
    %
    % `IMEigenvalueProblem` is the solver-independent contract for a
    % vertical-mode eigenvalue problem. A solver owns the stratification,
    % physical domain, coordinate mapping, and derivative matrices. The EVP
    % owns the physical constants, differential operators, boundary laws,
    % inner-product weights, normalization rules, equivalent-depth mapping,
    % and mode-index policy.
    %
    % Assembly combines those responsibilities in the solver native basis:
    % $$Aq_j=\lambda_jBq_j,\qquad h_j=\mathrm{hFromEigenvalue}(\lambda_j).$$
    % The retained columns become an `IMBasisSet`. The basis set evaluates the
    % solved variable and its linked diagnostic variable; `F` and `G` are not
    % independent mode families. In a `G` formulation,
    % $$F_j=h_j\partial_zG_j,$$
    % and in an `F` formulation,
    % $$G_j=-gN^{-2}\partial_zF_j.$$
    %
    % Use the static factories for standard wave and hydrostatic mode
    % problems. Use the constructor directly when defining a custom operator,
    % boundary, inner-product, or normalization contract.
    %
    % ```matlab
    % N2 = @(z) 1e-5*ones(size(z));
    % solver = IMSolverSpectral(N2=N2, zDomain=[-1000 0], nEVP=64);
    % evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=1e-4);
    % basisSet = solver.solveEVP(evp, nModes=4);
    % z = linspace(-1000, 0, 128).';
    % G = basisSet.G(z);
    % h = basisSet.h;
    % ```
    %
    % - Topic: Create standard EVPs
    % - Topic: Build custom EVPs
    % - Topic: Assemble EVPs
    % - Topic: Inspect EVP metadata
    % - Topic: Select retained modes
    % - Topic: Developer topics
    % - Declaration: classdef IMEigenvalueProblem

    properties (SetAccess = private)
        % Short EVP name.
        %
        % This is the canonical identity for factory-created EVPs, such as
        % `"waveModesAtWavenumber"` or `"hydrostaticGModes"`.
        %
        % - Topic: Inspect EVP metadata
        name = "unknown"
    end

    properties
        % Solved vertical-structure formulation.
        %
        % The formulation is either `"G"` or `"F"`. The basis set solves
        % this variable directly and obtains the other variable
        % diagnostically through
        % $$F_j=h_j\partial_zG_j$$ for `G` formulations or
        % $$G_j=-gN^{-2}\partial_zF_j$$ for `F` formulations.
        %
        % - Topic: Inspect EVP metadata
        formulation = "G"
    end

    properties (SetAccess = private)
        % Coriolis parameter.
        %
        % `f0` has units of radians per second. Solvers provide `N2` and the
        % vertical domain; the EVP provides the physical constants used by the
        % operators, boundaries, and normalization rules.
        %
        % - Topic: Inspect EVP metadata
        f0 = 0

        % Gravitational acceleration.
        %
        % `g` has units of meters per second squared. Operator coefficient and
        % diagnostic-variable functions read this value through the assembly
        % context as `ctx.g`.
        %
        % - Topic: Inspect EVP metadata
        g = 9.81
    end

    properties
        % Left differential operator.
        %
        % `leftOperator` contributes the matrix `A` in
        % $$Aq=\lambda Bq.$$
        % Coefficients are evaluated with the context returned by
        % `contextForSolver`, so coefficient functions may use solver-owned
        % fields such as `ctx.N2` and EVP-owned fields such as `ctx.g`.
        %
        % - Topic: Assemble EVPs
        leftOperator = IMOperator()

        % Right differential operator.
        %
        % `rightOperator` contributes the matrix `B` in the generalized EVP.
        % Standard factories use it for the weighted side of the strong form,
        % for example `@(z,ctx) (ctx.f0*ctx.f0 - ctx.N2(z))/ctx.g` in the
        % fixed-wavenumber wave problem.
        %
        % - Topic: Assemble EVPs
        rightOperator = IMOperator()

        % Inner-product weights for `F` and `G`.
        %
        % `innerWeights.F` and `innerWeights.G` are function handles with
        % signature `w = weight(z,ctx)`. Each handle returns the interior
        % weight in
        % $$\langle X_i,X_j\rangle_w=\int w(z)X_i(z)X_j(z)\,dz.$$
        % Boundary trace terms belong to the boundary conditions so that each
        % normalization can include the same endpoint convention as the EVP.
        %
        % - Topic: Inspect EVP metadata
        innerWeights = struct()

        % Named modal normalization rules.
        %
        % Each field stores a function handle with signature
        % `scale = rule(basisSet,iMode)`. The returned scale divides one raw
        % mode column after the solver has assembled, selected, and linked the
        % retained modes. Factory-created EVPs populate names such as
        % `unity`, `kConstant`, `omegaConstant`, `geostrophic`, `wMax`,
        % `uMax`, and `surfacePressure` when those rules are meaningful.
        %
        % - Topic: Inspect EVP metadata
        normalizations = struct()

        % Natural default normalization for this EVP.
        %
        % This is a `Normalization` value or `[]`. Empty means the EVP does
        % not declare a problem-specific default and the basis-set layer may
        % use its package fallback.
        %
        % - Topic: Inspect EVP metadata
        defaultNormalization = []

        % Placed boundary-condition array.
        %
        % The array stores location-aware `IMBoundary` values. Standard
        % factories accept location-free `surfaceBoundary` and `bottomBoundary`
        % laws and place them on the solved variable. During assembly each
        % placed boundary replaces the appropriate boundary row in `A` and `B`
        % and contributes endpoint metadata for indexing and normalization.
        %
        % - Topic: Assemble EVPs
        boundaryConditions = IMBoundary.empty(0,1)

        % Equivalent-depth conversion function.
        %
        % `hFromEigenvalue` is a function handle with signature
        % `h = hFromEigenvalue(lambda)`. Standard EVPs use
        % $$h_j=1/\lambda_j,$$
        % so `lambda` has inverse-depth units in those problems.
        %
        % - Topic: Inspect EVP metadata
        hFromEigenvalue = @(lambda) 1 ./ lambda
    end

    properties (SetAccess = private)
        % Whether the EVP declares the barotropic mode.
        %
        % When `true`, the mode-index policy expects the depth-uniform
        % zero-eigenvalue barotropic mode
        % $$F_0(z)=1,\qquad G_0(z)=0.$$
        % The barotropic mode is selected after boundary-index modes and
        % before positive baroclinic modes. Hydrostatic `F` modes declare this
        % mode; wave-mode and hydrostatic `G` EVPs do not.
        %
        % - Topic: Select retained modes
        hasBarotropicMode = false

        % Index validation behavior.
        %
        % Values are `"error"`, `"warning"`, or `"none"`. The mode-index
        % policy uses this value when the observed boundary, barotropic, or
        % interior counts differ from the counts declared by the EVP and its
        % boundary conditions.
        %
        % - Topic: Select retained modes
        indexValidationMode = "none"
    end

    properties
        % Stored factory-specific physical inputs.
        %
        % `parameters` records physical inputs supplied to a standard factory
        % but not otherwise stored as first-class EVP properties, such as
        % `parameters.k` or `parameters.omega`. The EVP identity is `name`,
        % boundary laws live in `boundaryConditions`, and physical constants
        % live in `f0` and `g`. Core EVP assembly does not consume this struct.
        %
        % - Topic: Inspect EVP metadata
        parameters = struct()
    end

    methods
        function self = IMEigenvalueProblem(options)
            % Create a physical-coordinate EVP descriptor.
            %
            % The constructor is the low-level entry point for custom EVPs.
            % Provide the solved `formulation`, the operators defining
            % $$Aq=\lambda Bq,$$
            % the placed boundary conditions for that formulation, the
            % interior inner-product weights, and any normalizations that an
            % `IMBasisSet` should expose. The standard factories are preferred
            % for the built-in wave and hydrostatic problems because they set
            % the operator, boundary, normalization, and mode-index metadata as
            % one coherent contract.
            %
            % ```matlab
            % left = IMOperator().plus(derivativeOrder=2);
            % right = IMOperator().plus(coefficient=@(z,ctx) -ctx.N2(z)/ctx.g, derivativeOrder=0);
            % evp = IMEigenvalueProblem(name="customG", formulation="G", ...
            %     leftOperator=left, rightOperator=right, ...
            %     boundaryConditions=IMBoundary.conditions(formulation="G"));
            % ```
            %
            % - Topic: Build custom EVPs
            % - Declaration: evp = IMEigenvalueProblem(options)
            % - Parameter options.name: short EVP name
            % - Parameter options.formulation: solved variable, `"F"` or `"G"`
            % - Parameter options.f0: Coriolis parameter in radians per second
            % - Parameter options.g: gravitational acceleration in meters per second squared
            % - Parameter options.leftOperator: operator that builds the generalized-EVP matrix `A`
            % - Parameter options.rightOperator: operator that builds the generalized-EVP matrix `B`
            % - Parameter options.innerWeights: interior inner-product weight handles for `F` and `G`
            % - Parameter options.normalizations: named modal normalization handles
            % - Parameter options.defaultNormalization: natural default normalization for this EVP
            % - Parameter options.boundaryConditions: placed boundary conditions
            % - Parameter options.hFromEigenvalue: equivalent-depth conversion
            % - Parameter options.hasBarotropicMode: whether the EVP declares the barotropic mode
            % - Parameter options.indexValidationMode: `"error"`, `"warning"`, or `"none"`
            % - Parameter options.parameters: stored factory-specific physical inputs
            % - Returns evp: initialized EVP descriptor
            arguments
                options.name {mustBeTextScalar} = "unknown"
                options.formulation {mustBeTextScalar} = "G"
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.leftOperator IMOperator = IMOperator()
                options.rightOperator IMOperator = IMOperator()
                options.innerWeights struct = struct()
                options.normalizations struct = struct()
                options.defaultNormalization = []
                options.boundaryConditions = IMBoundary.empty(0,1)
                options.hFromEigenvalue = @(lambda) 1 ./ lambda
                options.hasBarotropicMode (1,1) logical = false
                options.indexValidationMode {mustBeTextScalar} = "none"
                options.parameters struct = struct()
            end

            self.name = string(options.name);
            self.formulation = IMEigenvalueProblem.validateVariable(options.formulation);
            self.f0 = options.f0;
            self.g = options.g;
            self.leftOperator = options.leftOperator;
            self.rightOperator = options.rightOperator;
            self.innerWeights = IMEigenvalueProblem.resolveInnerWeights(options.innerWeights);
            self.normalizations = options.normalizations;
            self.defaultNormalization = options.defaultNormalization;
            self.boundaryConditions = options.boundaryConditions(:);
            self.validateBoundaryConditions();
            self.hFromEigenvalue = options.hFromEigenvalue;
            self.hasBarotropicMode = options.hasBarotropicMode;
            self.indexValidationMode = IMEigenvalueProblem.validateIndexValidationMode(options.indexValidationMode);
            self.parameters = options.parameters;
            metadataFieldsToRemove = intersect(fieldnames(self.parameters), {'f0'; 'g'});
            if ~isempty(metadataFieldsToRemove)
                self.parameters = rmfield(self.parameters, metadataFieldsToRemove);
            end
        end

        function [A, B] = assemble(self, solver)
            % Build generalized-EVP matrices on a solver's native basis.
            %
            % `assemble` evaluates `leftOperator` and `rightOperator` with the
            % merged solver/EVP context and returns the matrices for
            % $$Aq=\lambda Bq.$$
            % Interior rows come from the operator discretization. Boundary
            % rows are then replaced by the placed boundary conditions through
            % the solver, so a rigid `G` boundary imposes the trace row for
            % `G=0` while active or free boundaries can also declare endpoint
            % contributions used by normalization and mode indexing.
            %
            % - Topic: Assemble EVPs
            % - Declaration: [A,B] = assemble(evp,solver)
            % - Parameter solver: coordinate-aware internal-mode solver
            % - Returns A: left generalized-EVP matrix
            % - Returns B: right generalized-EVP matrix
            context = self.contextForSolver(solver);
            A = self.leftOperator.matrix(solver, context=context);
            B = self.rightOperator.matrix(solver, context=context);
            for iBoundary = 1:length(self.boundaryConditions)
                [A, B] = solver.applyBoundaryCondition(A, B, self.boundaryConditions(iBoundary), context=context);
            end
        end

        function context = contextForSolver(self, solver)
            % Return the coefficient context for this EVP and solver.
            %
            % The returned struct starts with the solver context, including
            % fields such as `N2`, `dzLogN2`, `zDomain`, and `coordinateKind`.
            % The EVP then adds physical constants as `f0` and `g`. Operator
            % coefficients, boundary rows, inner-product weights, and
            % normalization rules read this context but do not own the solver
            % discretization.
            %
            % - Topic: Assemble EVPs
            % - Declaration: context = contextForSolver(evp,solver)
            % - Parameter solver: coordinate-aware internal-mode solver
            % - Returns context: framework coefficient context
            context = solver.context();
            context.f0 = self.f0;
            context.g = self.g;
        end

        function selection = selectModes(self, eigenvalues, nModes, context)
            % Select and label retained eigenmodes.
            %
            % The index policy first classifies candidate eigenvalues using
            % boundary metadata, the optional barotropic mode, and positive
            % interior modes. Retained modes are ordered as boundary-index
            % modes, the barotropic mode when declared, then positive
            % baroclinic modes. Their labels define the `modeNumber` metadata
            % carried by the resulting `IMBasisSet`.
            %
            % - Topic: Select retained modes
            % - Topic: Developer topics
            % - Declaration: selection = selectModes(evp,eigenvalues,nModes,context)
            % - Parameter eigenvalues: candidate generalized-EVP eigenvalues
            % - Parameter nModes: number of modes to retain
            % - Parameter context: solver or analytical context
            % - Returns selection: selected-mode metadata
            % - Developer: true
            selection = self.indexPolicy().selectModes(eigenvalues, nModes, context);
        end

        function index = classifyEigenvalues(self, eigenvalues, context)
            % Classify eigenvalues using this EVP's index metadata.
            %
            % Classification reports the detected boundary-index directions,
            % the optional barotropic direction, positive interior modes, and
            % any validation mismatch. Negative boundary directions can be
            % expected by the declared boundary policy; they are not
            % automatically treated as numerical failures.
            %
            % - Topic: Select retained modes
            % - Topic: Developer topics
            % - Declaration: index = classifyEigenvalues(evp,eigenvalues,context)
            % - Parameter eigenvalues: generalized-EVP eigenvalues
            % - Parameter context: solver or analytical context
            % - Returns index: index summary structure
            % - Developer: true
            index = self.indexPolicy().classify(eigenvalues, context);
        end
    end

    methods (Access = private)
        function policy = indexPolicy(self)
            policy = IMIndexPolicy.fromBoundaryConditions(self.boundaryConditions, ...
                expectedZeroCount=double(self.hasBarotropicMode), validationMode=self.indexValidationMode);
        end

        function validateBoundaryConditions(self)
            for iBoundary = 1:length(self.boundaryConditions)
                boundaryCondition = self.boundaryConditions(iBoundary);
                if boundaryCondition.family == "active" || boundaryCondition.family == "partialDepthPE"
                    continue;
                end
                if boundaryCondition.location ~= "" && boundaryCondition.variable ~= self.formulation
                    error("IMEigenvalueProblem:BoundaryFormulationMismatch", ...
                        "Boundary condition ""%s"" targets variable ""%s"", but the EVP formulation is ""%s"".", ...
                        boundaryCondition.family, boundaryCondition.variable, self.formulation);
                end
            end
        end
    end

    methods (Static)
        function evp = waveModesAtWavenumber(options)
            % Create the wave-mode `G` EVP at fixed horizontal wavenumber.
            %
            % This factory fixes the horizontal wavenumber `K=options.k` and
            % solves the physical-coordinate strong form
            % $$G_{zz}-K^2G=\lambda(f_0^2-N^2)G/g,\qquad h=1/\lambda.$$
            % The solved variable is `G`; the linked diagnostic variable is
            % $$F=hG_z.$$
            % The factory stores `parameters.k`, uses the default
            % `Normalization.kConstant` normalization, and places the supplied
            % location-free boundary laws on the surface and bottom. Omitted
            % boundaries are rigid `G=0` boundaries.
            %
            % ```matlab
            % solver = IMSolverSpectral(N2=@(z) 1e-5*ones(size(z)), zDomain=[-500 0], nEVP=48);
            % evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=1e-4);
            % basisSet = solver.solveEVP(evp, nModes=6);
            % h = basisSet.h;
            % modeNumber = basisSet.modeNumber;
            % ```
            %
            % - Topic: Create standard EVPs
            % - Declaration: evp = IMEigenvalueProblem.waveModesAtWavenumber(options)
            % - Parameter options.k: horizontal wavenumber in radians per meter
            % - Parameter options.f0: Coriolis parameter in radians per second
            % - Parameter options.g: gravitational acceleration in meters per second squared
            % - Parameter options.surfaceBoundary: location-free surface boundary law
            % - Parameter options.bottomBoundary: location-free bottom boundary law
            % - Returns evp: fixed-wavenumber wave-mode `G` EVP
            arguments
                options.k (1,1) double {mustBeNonnegative}
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundary = IMBoundary.rigid()
                options.bottomBoundary (1,1) IMBoundary = IMBoundary.rigid()
            end

            k = options.k;
            f0 = options.f0;
            g = options.g;
            left = IMOperator().plus(derivativeOrder=2).plus(coefficient=-k*k, derivativeOrder=0);
            right = IMOperator().plus(coefficient=@(z,ctx) (ctx.f0*ctx.f0 - ctx.N2(z))/ctx.g, derivativeOrder=0);
            boundaryConditions = IMBoundary.conditions(formulation="G", surface=options.surfaceBoundary, bottom=options.bottomBoundary);
            parameters = struct("k", k);
            innerWeights.G = @(z,ctx) (ctx.N2(z) - ctx.f0*ctx.f0)/ctx.g;
            evp = IMEigenvalueProblem(name="waveModesAtWavenumber", formulation="G", ...
                f0=f0, g=g, ...
                leftOperator=left, rightOperator=right, innerWeights=innerWeights, boundaryConditions=boundaryConditions, ...
                defaultNormalization=Normalization.kConstant, ...
                hFromEigenvalue=@(lambda) 1 ./ lambda, parameters=parameters);

            evp.normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor("G", iMode);
            evp.normalizations.kConstant = @(basisSet,iMode) basisSet.innerProductNormFactor("G", iMode);
            evp.normalizations.omegaConstant = @(basisSet,iMode) basisSet.innerProductNormFactor("F", iMode);
            evp.normalizations.wMax = @(basisSet,iMode) basisSet.maxAbsFactor("G", iMode);
            evp.normalizations.uMax = @(basisSet,iMode) basisSet.maxAbsFactor("F", iMode);
            evp.normalizations.surfacePressure = @(basisSet,iMode) basisSet.surfacePressureNormFactor(iMode);
        end

        function evp = waveModesAtFrequency(options)
            % Create the wave-mode `G` EVP at fixed frequency.
            %
            % This factory fixes the frequency `omega=options.omega` and
            % solves the physical-coordinate strong form
            % $$G_{zz}=\lambda(\omega^2-N^2)G/g,\qquad h=1/\lambda.$$
            % The solved variable is `G`; the linked diagnostic variable is
            % $$F=hG_z.$$
            % The factory stores `parameters.omega`, uses the default
            % `Normalization.omegaConstant` normalization, and places the
            % supplied location-free boundary laws on the surface and bottom.
            % Omitted boundaries are rigid `G=0` boundaries. The `kConstant`
            % normalization remains available for fixed-frequency basis sets
            % and includes boundary trace terms when the boundary laws provide
            % them.
            %
            % ```matlab
            % solver = IMSolverSpectral(N2=@(z) 1e-5*ones(size(z)), zDomain=[-500 0], nEVP=48);
            % evp = IMEigenvalueProblem.waveModesAtFrequency(omega=1.2e-3, f0=1e-4);
            % basisSet = solver.solveEVP(evp, nModes=6);
            % F = basisSet.F(linspace(-500, 0, 100).');
            % ```
            %
            % - Topic: Create standard EVPs
            % - Declaration: evp = IMEigenvalueProblem.waveModesAtFrequency(options)
            % - Parameter options.omega: fixed frequency in radians per second
            % - Parameter options.f0: Coriolis parameter in radians per second
            % - Parameter options.g: gravitational acceleration in meters per second squared
            % - Parameter options.surfaceBoundary: location-free surface boundary law
            % - Parameter options.bottomBoundary: location-free bottom boundary law
            % - Returns evp: fixed-frequency wave-mode `G` EVP
            arguments
                options.omega (1,1) double {mustBeNonnegative}
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundary = IMBoundary.rigid()
                options.bottomBoundary (1,1) IMBoundary = IMBoundary.rigid()
            end

            omega = options.omega;
            f0 = options.f0;
            g = options.g;
            left = IMOperator().plus(derivativeOrder=2);
            right = IMOperator().plus(coefficient=@(z,ctx) (omega*omega - ctx.N2(z))/ctx.g, derivativeOrder=0);
            boundaryConditions = IMBoundary.conditions(formulation="G", surface=options.surfaceBoundary, bottom=options.bottomBoundary);
            parameters = struct("omega", omega);
            innerWeights.G = @(z,ctx) (ctx.N2(z) - omega*omega)/ctx.g;
            innerWeights.F = @(z,ctx) ones(size(z))/diff(ctx.zDomain);
            evp = IMEigenvalueProblem(name="waveModesAtFrequency", formulation="G", ...
                f0=f0, g=g, ...
                leftOperator=left, rightOperator=right, innerWeights=innerWeights, boundaryConditions=boundaryConditions, ...
                defaultNormalization=Normalization.omegaConstant, ...
                hFromEigenvalue=@(lambda) 1 ./ lambda, parameters=parameters);

            evp.normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor("G", iMode);
            evp.normalizations.kConstant = @(basisSet,iMode) basisSet.weightedNormFactorWithBoundaryTerms("G", iMode, @(z,ctx) (ctx.N2(z) - ctx.f0*ctx.f0)/ctx.g);
            evp.normalizations.omegaConstant = @(basisSet,iMode) basisSet.innerProductNormFactor("F", iMode);
            evp.normalizations.wMax = @(basisSet,iMode) basisSet.maxAbsFactor("G", iMode);
            evp.normalizations.uMax = @(basisSet,iMode) basisSet.maxAbsFactor("F", iMode);
            evp.normalizations.surfacePressure = @(basisSet,iMode) basisSet.surfacePressureNormFactor(iMode);
        end

        function evp = hydrostaticGModes(options)
            % Create the hydrostatic `G`-mode EVP.
            %
            % Hydrostatic `G` modes are the zero-frequency wave-mode problem
            % written directly as
            % $$G_{zz}=-\lambda N^2G/g,\qquad h=1/\lambda.$$
            % The solved variable is `G`; the linked diagnostic variable is
            % $$F=hG_z.$$
            % There is no nontrivial null `G` mode, so retained modes are the
            % boundary-index modes declared by the boundary laws followed by
            % positive interior baroclinic modes. The default normalization is
            % `Normalization.geostrophic`.
            %
            % ```matlab
            % solver = IMSolverSpectral(N2=@(z) 1e-5*ones(size(z)), zDomain=[-1000 0], nEVP=64);
            % evp = IMEigenvalueProblem.hydrostaticGModes();
            % basisSet = solver.solveEVP(evp, nModes=4);
            % G = basisSet.G(linspace(-1000, 0, 128).');
            % ```
            %
            % - Topic: Create standard EVPs
            % - Declaration: evp = IMEigenvalueProblem.hydrostaticGModes(options)
            % - Parameter options.f0: Coriolis parameter in radians per second
            % - Parameter options.g: gravitational acceleration in meters per second squared
            % - Parameter options.surfaceBoundary: location-free surface boundary law
            % - Parameter options.bottomBoundary: location-free bottom boundary law
            % - Returns evp: zero-frequency hydrostatic `G` EVP
            arguments
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundary = IMBoundary.rigid()
                options.bottomBoundary (1,1) IMBoundary = IMBoundary.rigid()
            end

            g = options.g;
            left = IMOperator().plus(derivativeOrder=2);
            right = IMOperator().plus(coefficient=@(z,ctx) -ctx.N2(z)/ctx.g, derivativeOrder=0);
            boundaryConditions = IMBoundary.conditions(formulation="G", surface=options.surfaceBoundary, bottom=options.bottomBoundary);
            innerWeights.G = @(z,ctx) ctx.N2(z)/ctx.g;
            innerWeights.F = @(z,ctx) ones(size(z));
            evp = IMEigenvalueProblem(name="hydrostaticGModes", formulation="G", ...
                f0=options.f0, g=g, ...
                leftOperator=left, rightOperator=right, innerWeights=innerWeights, boundaryConditions=boundaryConditions, ...
                defaultNormalization=Normalization.geostrophic, hFromEigenvalue=@(lambda) 1 ./ lambda);

            evp.normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor("G", iMode);
            evp.normalizations.geostrophic = @(basisSet,iMode) basisSet.geostrophicNormFactor(iMode);
            evp.normalizations.wMax = @(basisSet,iMode) basisSet.maxAbsFactor("G", iMode);
            evp.normalizations.uMax = @(basisSet,iMode) basisSet.maxAbsFactor("F", iMode);
            evp.normalizations.surfacePressure = @(basisSet,iMode) basisSet.surfacePressureNormFactor(iMode);
        end

        function evp = hydrostaticFModes(options)
            % Create the geostrophic hydrostatic `F`-mode EVP.
            %
            % The physical-coordinate strong form is
            % $$F_{zz}-(\partial_z\log N^2)F_z=-\lambda N^2F/g,\qquad h=1/\lambda.$$
            % The solved variable is `F`; the linked diagnostic variable is
            % $$G=-gN^{-2}F_z.$$
            % This EVP declares the barotropic mode,
            % $$F_0(z)=1,\qquad G_0(z)=0,$$
            % so the barotropic mode is retained before the positive
            % baroclinic modes. The default normalization is
            % `Normalization.geostrophic`.
            %
            % ```matlab
            % solver = IMSolverSpectral(N2=@(z) 1e-5*ones(size(z)), zDomain=[-1000 0], nEVP=64);
            % evp = IMEigenvalueProblem.hydrostaticFModes();
            % basisSet = solver.solveEVP(evp, nModes=4);
            % F = basisSet.F(linspace(-1000, 0, 128).');
            % ```
            %
            % - Topic: Create standard EVPs
            % - Declaration: evp = IMEigenvalueProblem.hydrostaticFModes(options)
            % - Parameter options.g: gravitational acceleration in meters per second squared
            % - Parameter options.surfaceBoundary: location-free surface boundary law
            % - Parameter options.bottomBoundary: location-free bottom boundary law
            % - Returns evp: hydrostatic `F` EVP
            arguments
                options.g (1,1) double {mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundary = IMBoundary.rigid()
                options.bottomBoundary (1,1) IMBoundary = IMBoundary.rigid()
            end

            g = options.g;
            left = IMOperator().plus(derivativeOrder=2) ...
                .plus(coefficient=@(z,ctx) -ctx.dzLogN2(z), derivativeOrder=1);
            right = IMOperator().plus(coefficient=@(z,ctx) -ctx.N2(z)/ctx.g, derivativeOrder=0);
            boundaryConditions = IMBoundary.conditions(formulation="F", surface=options.surfaceBoundary, bottom=options.bottomBoundary);
            innerWeights.F = @(z,ctx) ones(size(z));
            innerWeights.G = @(z,ctx) ctx.N2(z)/ctx.g;
            evp = IMEigenvalueProblem(name="hydrostaticFModes", formulation="F", ...
                g=g, ...
                leftOperator=left, rightOperator=right, innerWeights=innerWeights, boundaryConditions=boundaryConditions, ...
                defaultNormalization=Normalization.geostrophic, ...
                hFromEigenvalue=@(lambda) 1 ./ lambda, hasBarotropicMode=true, indexValidationMode="warning");
            evp.normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor("F", iMode);
            evp.normalizations.geostrophic = @(basisSet,iMode) basisSet.geostrophicNormFactor(iMode);
            evp.normalizations.wMax = @(basisSet,iMode) basisSet.maxAbsFactor("G", iMode);
            evp.normalizations.uMax = @(basisSet,iMode) basisSet.maxAbsFactor("F", iMode);
            evp.normalizations.surfacePressure = @(basisSet,iMode) basisSet.surfacePressureNormFactor(iMode);
        end

        function policy = partialDepthPEIndexPolicy(options)
            % Return the partial-depth potential-energy index policy.
            %
            % Partial-depth potential-energy forms can declare endpoint
            % directions through active boundary metadata. With
            % `boundarySign="positive"`, endpoint contributions use the
            % positive boundary convention. With `boundarySign="negative"`,
            % endpoint contributions use the negative boundary convention and
            % the resulting negative index directions are expected by the
            % policy rather than treated as numerical failures.
            %
            % - Topic: Select retained modes
            % - Declaration: policy = IMEigenvalueProblem.partialDepthPEIndexPolicy(options)
            % - Parameter options.boundarySign: `"positive"` or `"negative"`
            % - Parameter options.validationMode: `"error"`, `"warning"`, or `"none"`
            % - Returns policy: corresponding index policy
            arguments
                options.boundarySign {mustBeTextScalar} = "positive"
                options.validationMode {mustBeTextScalar} = "error"
            end

            boundaryConditions = IMBoundary.partialDepthPE(boundarySign=options.boundarySign);
            policy = IMIndexPolicy.fromBoundaryConditions(boundaryConditions, validationMode=options.validationMode);
        end
    end

    methods (Static, Access = private)
        function variable = validateVariable(variable)
            variable = string(variable);
            if variable ~= "F" && variable ~= "G"
                error("IMEigenvalueProblem:InvalidVariable", ...
                    "formulation must be ""F"" or ""G"".");
            end
        end

        function innerWeights = resolveInnerWeights(innerWeights)
            if ~isfield(innerWeights, "G")
                innerWeights.G = @(z,ctx) ctx.N2(z)/ctx.g;
            end
            if ~isfield(innerWeights, "F")
                innerWeights.F = @(z,ctx) ones(size(z));
            end
        end

        function validationMode = validateIndexValidationMode(validationMode)
            validationMode = string(validationMode);
            if ~ismember(validationMode, ["error", "warning", "none"])
                error("IMEigenvalueProblem:InvalidIndexValidationMode", ...
                    "indexValidationMode must be ""error"", ""warning"", or ""none"".");
            end
        end
    end

end

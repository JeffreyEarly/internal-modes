classdef IMEigenvalueProblem
    % Describe an internal-mode eigenvalue problem in physical coordinates.
    %
    % `IMEigenvalueProblem` stores the physical operators, boundary laws,
    % formulation, inner-product weights, named normalization rules,
    % default normalization, equivalent-depth interpretation, and
    % index-selection metadata needed by coordinate-aware solvers.
    %
    % ```matlab
    % evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=1e-4);
    % ```
    %
    % - Topic: Create EVPs
    % - Topic: Assemble EVPs
    % - Topic: Inspect EVP metadata
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
        % Coriolis parameter owned by this EVP.
        %
        % Solvers provide the medium and discretization; the EVP owns the
        % physical constants that define the mathematical problem.
        %
        % - Topic: Inspect EVP metadata
        f0 = 0

        % Gravitational acceleration owned by this EVP.
        %
        % - Topic: Inspect EVP metadata
        g = 9.81
    end

    properties
        % Left physical operator.
        %
        % - Topic: Assemble EVPs
        leftOperator = IMOperator()

        % Right physical operator.
        %
        % - Topic: Assemble EVPs
        rightOperator = IMOperator()

        % Inner-product weights for `F` and `G`.
        %
        % Each field stores the interior weight in
        % $$\int w(z)X_i(z)X_j(z)\,dz$$ for variable `F` or `G`. Boundary
        % trace terms are stored on boundary conditions.
        %
        % - Topic: Inspect EVP metadata
        innerWeights = struct()

        % Named modal normalization rules.
        %
        % Each field stores a function handle that returns the divisor for
        % one raw mode column.
        %
        % - Topic: Inspect EVP metadata
        normalizations = struct()

        % Natural default normalization for this EVP.
        %
        % Empty means the EVP does not declare a problem-specific default,
        % and basis sets should use the package fallback.
        %
        % - Topic: Inspect EVP metadata
        defaultNormalization = []

        % Placed boundary-condition array.
        %
        % - Topic: Assemble EVPs
        boundaryConditions = IMBoundary.empty(0,1)

        % Equivalent-depth conversion function.
        %
        % This maps generalized-EVP eigenvalues $$\lambda$$ to equivalent
        % depths $$h_j$$.
        %
        % - Topic: Inspect EVP metadata
        hFromEigenvalue = @(lambda) 1 ./ lambda
    end

    properties (SetAccess = private)
        % Number of true null modes.
        %
        % Null modes are zero-eigenvalue modes selected after boundary
        % modes and before positive interior modes. In the hydrostatic `F`
        % EVP, this is the depth-uniform mode with $$F_0(z)=1$$ and
        % $$G_0(z)=0$$.
        %
        % - Topic: Inspect EVP metadata
        nNullModes = 0

        % Index validation behavior.
        %
        % Values are `"error"`, `"warning"`, or `"none"`.
        %
        % - Topic: Inspect EVP metadata
        indexValidationMode = "none"
    end

    properties
        % Stored factory-specific physical inputs.
        %
        % `parameters` records physical inputs supplied to a factory but not
        % otherwise stored as first-class EVP properties, such as `k` or
        % `omega`. The EVP identity is `name`, boundary laws live in
        % `boundaryConditions`, and physical constants live in `f0` and `g`.
        % Core EVP assembly does not consume this struct.
        %
        % - Topic: Inspect EVP metadata
        parameters = struct()
    end

    methods
        function self = IMEigenvalueProblem(options)
            % Create a physical-coordinate EVP descriptor.
            %
            % - Topic: Create EVPs
            % - Declaration: evp = IMEigenvalueProblem(options)
            % - Parameter options.name: short EVP name
            % - Parameter options.formulation: solved variable, `"F"` or `"G"`
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.leftOperator: left physical operator
            % - Parameter options.rightOperator: right physical operator
            % - Parameter options.innerWeights: interior inner-product weights for `F` and `G`
            % - Parameter options.normalizations: named modal normalization rules
            % - Parameter options.defaultNormalization: natural default normalization for this EVP
            % - Parameter options.boundaryConditions: placed boundary conditions
            % - Parameter options.hFromEigenvalue: equivalent-depth conversion
            % - Parameter options.nNullModes: number of true null modes
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
                options.nNullModes (1,1) double {mustBeInteger, mustBeNonnegative} = 0
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
            self.nNullModes = options.nNullModes;
            self.indexValidationMode = IMEigenvalueProblem.validateIndexValidationMode(options.indexValidationMode);
            self.parameters = options.parameters;
            metadataFieldsToRemove = intersect(fieldnames(self.parameters), {'f0'; 'g'});
            if ~isempty(metadataFieldsToRemove)
                self.parameters = rmfield(self.parameters, metadataFieldsToRemove);
            end
        end

        function [A, B] = assemble(self, solver)
            % Assemble the EVP on a solver's native basis.
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
            % The solver supplies medium and discretization fields; the EVP
            % supplies physical constants.
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
            % Boundary conditions provide boundary-mode index metadata.
            % `nNullModes` provides expected zero-eigenvalue null modes. The
            % selected modes are ordered as boundary modes, null modes, then
            % positive interior modes.
            %
            % - Topic: Assemble EVPs
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
            % - Topic: Assemble EVPs
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
                expectedZeroCount=self.nNullModes, validationMode=self.indexValidationMode);
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
            % The physical-coordinate strong form is
            % $$G_{zz}-K^2G=\lambda(f_0^2-N^2)G/g$$.
            %
            % - Topic: Create EVPs
            % - Declaration: evp = IMEigenvalueProblem.waveModesAtWavenumber(options)
            % - Parameter options.k: horizontal wavenumber
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
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
        end

        function evp = waveModesAtFrequency(options)
            % Create the wave-mode `G` EVP at fixed frequency.
            %
            % The physical-coordinate strong form is
            % $$G_{zz}=\lambda(\omega^2-N^2)G/g$$.
            %
            % - Topic: Create EVPs
            % - Declaration: evp = IMEigenvalueProblem.waveModesAtFrequency(options)
            % - Parameter options.omega: fixed frequency in radians per second
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
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
        end

        function evp = hydrostaticGModes(options)
            % Create the hydrostatic `G`-mode EVP.
            %
            % Hydrostatic `G` modes satisfy
            % $$G_{zz}=-\lambda N^2G/g$$ and have no nontrivial null `G`
            % mode.
            %
            % - Topic: Create EVPs
            % - Declaration: evp = IMEigenvalueProblem.hydrostaticGModes(options)
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
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
        end

        function evp = hydrostaticFModes(options)
            % Create the geostrophic hydrostatic `F`-mode EVP.
            %
            % The physical-coordinate strong form is
            % $$F_{zz}-(\partial_z\log N^2)F_z=-\lambda N^2F/g$$.
            %
            % - Topic: Create EVPs
            % - Declaration: evp = IMEigenvalueProblem.hydrostaticFModes(options)
            % - Parameter options.g: gravitational acceleration
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
                hFromEigenvalue=@(lambda) 1 ./ lambda, nNullModes=1, indexValidationMode="warning");
            evp.normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor("F", iMode);
            evp.normalizations.geostrophic = @(basisSet,iMode) basisSet.geostrophicNormFactor(iMode);
            evp.normalizations.wMax = @(basisSet,iMode) basisSet.maxAbsFactor("G", iMode);
            evp.normalizations.uMax = @(basisSet,iMode) basisSet.maxAbsFactor("F", iMode);
        end

        function policy = partialDepthPEIndexPolicy(options)
            % Return the manuscript partial-depth PE index policy.
            %
            % - Topic: Create EVPs
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

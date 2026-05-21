classdef IMEigenvalueProblem
    % Describe an internal-mode eigenvalue problem in physical coordinates.
    %
    % `IMEigenvalueProblem` stores the physical operators, boundary laws,
    % component metadata, named normalization rules, default
    % normalization, eigenvalue interpretation, and index policy needed by
    % coordinate-aware solvers.
    %
    % ```matlab
    % evp = IMEigenvalueProblem.waveModesAtWavenumber(k=1e-4, f0=1e-4);
    % ```
    %
    % - Topic: Create EVPs
    % - Topic: Assemble EVPs
    % - Topic: Inspect EVP metadata
    % - Declaration: classdef IMEigenvalueProblem

    properties
        % Short EVP name.
        %
        % - Topic: Inspect EVP metadata
        name = "unknown"

        % Primary eigenfunction component.
        %
        % - Topic: Inspect EVP metadata
        primaryComponent = "G"
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
        leftOperator = IMOperator.strong()

        % Right physical operator.
        %
        % - Topic: Assemble EVPs
        rightOperator = IMOperator.strong()

        % Component metadata.
        %
        % - Topic: Inspect EVP metadata
        components = struct()

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

        % Resolved boundary-row array.
        %
        % - Topic: Assemble EVPs
        boundaryRows = IMBoundaryRow.empty(0,1)

        % Eigenvalue name.
        %
        % - Topic: Inspect EVP metadata
        eigenvalueName = "lambda"

        % Equivalent-depth conversion function.
        %
        % - Topic: Inspect EVP metadata
        hFromEigenvalue = @(lambda) 1 ./ lambda

        % Eigenvalue ordering rule.
        %
        % - Topic: Inspect EVP metadata
        ordering = "ascendingEigenvalue"

        % Expected index policy.
        %
        % - Topic: Inspect EVP metadata
        indexPolicy = IMIndexPolicy.none()

        % Stored problem metadata.
        %
        % `parameters` records provenance and factory-specific values for
        % inspection, tests, and analytical shortcuts. Core EVP assembly
        % uses the operators, boundary rows, and the EVP-built coefficient
        % context, not this metadata. Factory-created EVPs use
        % `problemType`, `upperBoundary`, and `lowerBoundary`, with optional
        % fields such as `k` and `omega` when they were supplied.
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
            % - Parameter options.primaryComponent: primary eigenfunction component
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.leftOperator: left physical operator
            % - Parameter options.rightOperator: right physical operator
            % - Parameter options.normalizations: named modal normalization rules
            % - Parameter options.defaultNormalization: natural default normalization for this EVP
            % - Parameter options.boundaryRows: resolved boundary rows
            % - Parameter options.eigenvalueName: eigenvalue name
            % - Parameter options.hFromEigenvalue: equivalent-depth conversion
            % - Parameter options.ordering: eigenvalue ordering rule
            % - Parameter options.indexPolicy: expected index policy
            % - Parameter options.parameters: stored problem metadata for inspection and analytical shortcuts
            % - Returns evp: initialized EVP descriptor
            arguments
                options.name {mustBeTextScalar} = "unknown"
                options.primaryComponent {mustBeTextScalar} = "G"
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.leftOperator IMOperator = IMOperator.strong()
                options.rightOperator IMOperator = IMOperator.strong()
                options.normalizations struct = struct()
                options.defaultNormalization = []
                options.boundaryRows = IMBoundaryRow.empty(0,1)
                options.eigenvalueName {mustBeTextScalar} = "lambda"
                options.hFromEigenvalue = @(lambda) 1 ./ lambda
                options.ordering {mustBeTextScalar} = "ascendingEigenvalue"
                options.indexPolicy IMIndexPolicy = IMIndexPolicy.none()
                options.parameters struct = struct()
            end

            self.name = string(options.name);
            self.primaryComponent = string(options.primaryComponent);
            self.f0 = options.f0;
            self.g = options.g;
            self.leftOperator = options.leftOperator;
            self.rightOperator = options.rightOperator;
            self.normalizations = options.normalizations;
            self.defaultNormalization = options.defaultNormalization;
            self.boundaryRows = options.boundaryRows(:);
            self.eigenvalueName = string(options.eigenvalueName);
            self.hFromEigenvalue = options.hFromEigenvalue;
            self.ordering = string(options.ordering);
            self.indexPolicy = options.indexPolicy;
            self.parameters = options.parameters;
            metadataFieldsToRemove = intersect(fieldnames(self.parameters), {'f0'; 'g'});
            if ~isempty(metadataFieldsToRemove)
                self.parameters = rmfield(self.parameters, metadataFieldsToRemove);
            end
            self.components = struct();
            self.components.(char(self.primaryComponent)).role = "eigenfunction";
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
            for iBoundary = 1:length(self.boundaryRows)
                [A, B] = self.boundaryRows(iBoundary).apply(A, B, solver, context=context);
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
            % - Parameter options.upperBoundary: location-free upper boundary law
            % - Parameter options.lowerBoundary: location-free lower boundary law
            % - Returns evp: fixed-wavenumber wave-mode `G` EVP
            arguments
                options.k (1,1) double {mustBeNonnegative}
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.upperBoundary (1,1) IMBoundary = IMBoundary.rigid()
                options.lowerBoundary (1,1) IMBoundary = IMBoundary.rigid()
            end

            k = options.k;
            f0 = options.f0;
            g = options.g;
            left = IMOperator.strong().plus(derivativeOrder=2).plus(coefficient=-k*k, derivativeOrder=0);
            right = IMOperator.strong().plus(coefficient=@(z,ctx) (ctx.f0*ctx.f0 - ctx.N2(z))/ctx.g, derivativeOrder=0);
            lowerRow = options.lowerBoundary.resolve(endpoint="lower", primaryComponent="G");
            upperRow = options.upperBoundary.resolve(endpoint="upper", primaryComponent="G");
            boundaryRows = [lowerRow; upperRow];
            parameters = struct("problemType", "waveModesAtWavenumber", "k", k, "omega", NaN, ...
                "upperBoundary", upperRow.family, "lowerBoundary", lowerRow.family);
            evp = IMEigenvalueProblem(name="waveModesAtWavenumber", primaryComponent="G", ...
                f0=f0, g=g, ...
                leftOperator=left, rightOperator=right, boundaryRows=boundaryRows, ...
                defaultNormalization=Normalization.kConstant, ...
                eigenvalueName="inverseEquivalentDepth", hFromEigenvalue=@(lambda) 1 ./ lambda, ...
                ordering="ascendingEigenvalue", ...
                indexPolicy=IMIndexPolicy.fromBoundaryRows(boundaryRows, validationMode="none"), parameters=parameters);

            evp.components.G.innerWeight = @(z,ctx) (ctx.N2(z) - ctx.f0*ctx.f0)/ctx.g;
            evp.components.F.role = "diagnostic";
            evp.components.F.from = "G";
            evp.components.F.operator = IMOperator.strong().plus(derivativeOrder=1);
            evp.components.F.modalScale = "h";

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
            % - Parameter options.upperBoundary: location-free upper boundary law
            % - Parameter options.lowerBoundary: location-free lower boundary law
            % - Returns evp: fixed-frequency wave-mode `G` EVP
            arguments
                options.omega (1,1) double {mustBeNonnegative}
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.upperBoundary (1,1) IMBoundary = IMBoundary.rigid()
                options.lowerBoundary (1,1) IMBoundary = IMBoundary.rigid()
            end

            omega = options.omega;
            f0 = options.f0;
            g = options.g;
            left = IMOperator.strong().plus(derivativeOrder=2);
            right = IMOperator.strong().plus(coefficient=@(z,ctx) (omega*omega - ctx.N2(z))/ctx.g, derivativeOrder=0);
            lowerRow = options.lowerBoundary.resolve(endpoint="lower", primaryComponent="G");
            upperRow = options.upperBoundary.resolve(endpoint="upper", primaryComponent="G");
            boundaryRows = [lowerRow; upperRow];
            parameters = struct("problemType", "waveModesAtFrequency", "k", NaN, "omega", omega, ...
                "upperBoundary", upperRow.family, "lowerBoundary", lowerRow.family);
            evp = IMEigenvalueProblem(name="waveModesAtFrequency", primaryComponent="G", ...
                f0=f0, g=g, ...
                leftOperator=left, rightOperator=right, boundaryRows=boundaryRows, ...
                defaultNormalization=Normalization.omegaConstant, ...
                eigenvalueName="inverseEquivalentDepth", hFromEigenvalue=@(lambda) 1 ./ lambda, ...
                ordering="ascendingEigenvalue", ...
                indexPolicy=IMIndexPolicy.fromBoundaryRows(boundaryRows, validationMode="none"), parameters=parameters);
            evp.components.G.innerWeight = @(z,ctx) (ctx.N2(z) - omega*omega)/ctx.g;
            evp.components.F.role = "diagnostic";
            evp.components.F.from = "G";
            evp.components.F.operator = IMOperator.strong().plus(derivativeOrder=1);
            evp.components.F.modalScale = "h";
            evp.components.F.innerWeight = @(z,ctx) ones(size(z))/diff(ctx.zDomain);

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
            % $$G_{zz}=-\lambda N^2G/g$$ and have no nontrivial
            % barotropic `G` mode.
            %
            % - Topic: Create EVPs
            % - Declaration: evp = IMEigenvalueProblem.hydrostaticGModes(options)
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.upperBoundary: location-free upper boundary law
            % - Parameter options.lowerBoundary: location-free lower boundary law
            % - Returns evp: zero-frequency hydrostatic `G` EVP
            arguments
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.upperBoundary (1,1) IMBoundary = IMBoundary.rigid()
                options.lowerBoundary (1,1) IMBoundary = IMBoundary.rigid()
            end

            f0 = options.f0;
            g = options.g;
            left = IMOperator.strong().plus(derivativeOrder=2);
            right = IMOperator.strong().plus(coefficient=@(z,ctx) -ctx.N2(z)/ctx.g, derivativeOrder=0);
            lowerRow = options.lowerBoundary.resolve(endpoint="lower", primaryComponent="G");
            upperRow = options.upperBoundary.resolve(endpoint="upper", primaryComponent="G");
            boundaryRows = [lowerRow; upperRow];
            parameters = struct("problemType", "hydrostaticGModes", "omega", 0, ...
                "upperBoundary", upperRow.family, "lowerBoundary", lowerRow.family);
            evp = IMEigenvalueProblem(name="hydrostaticGModes", primaryComponent="G", ...
                f0=f0, g=g, ...
                leftOperator=left, rightOperator=right, boundaryRows=boundaryRows, ...
                defaultNormalization=Normalization.geostrophic, ...
                eigenvalueName="inverseEquivalentDepth", hFromEigenvalue=@(lambda) 1 ./ lambda, ...
                ordering="indexPolicyThenAscending", ...
                indexPolicy=IMIndexPolicy.fromBoundaryRows(boundaryRows, validationMode="none"), parameters=parameters);
            evp.components.G.innerWeight = @(z,ctx) ctx.N2(z)/ctx.g;
            evp.components.F.role = "diagnostic";
            evp.components.F.from = "G";
            evp.components.F.operator = IMOperator.strong().plus(derivativeOrder=1);
            evp.components.F.modalScale = "h";
            evp.components.F.innerWeight = @(z,ctx) ones(size(z));

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
            % - Parameter options.upperBoundary: location-free upper boundary law
            % - Parameter options.lowerBoundary: location-free lower boundary law
            % - Returns evp: hydrostatic `F` EVP
            arguments
                options.g (1,1) double {mustBePositive} = 9.81
                options.upperBoundary (1,1) IMBoundary = IMBoundary.rigid()
                options.lowerBoundary (1,1) IMBoundary = IMBoundary.rigid()
            end

            g = options.g;
            left = IMOperator.strong().plus(derivativeOrder=2) ...
                .plus(coefficient=@(z,ctx) -ctx.dzLogN2(z), derivativeOrder=1);
            right = IMOperator.strong().plus(coefficient=@(z,ctx) -ctx.N2(z)/ctx.g, derivativeOrder=0);
            lowerRow = options.lowerBoundary.resolve(endpoint="lower", primaryComponent="F");
            upperRow = options.upperBoundary.resolve(endpoint="upper", primaryComponent="F");
            boundaryRows = [lowerRow; upperRow];
            evp = IMEigenvalueProblem(name="hydrostaticFModes", primaryComponent="F", ...
                g=g, ...
                leftOperator=left, rightOperator=right, boundaryRows=boundaryRows, ...
                defaultNormalization=Normalization.geostrophic, ...
                eigenvalueName="inverseEquivalentDepth", hFromEigenvalue=@(lambda) 1 ./ lambda, ...
                ordering="indexPolicyThenAscending", ...
                indexPolicy=IMIndexPolicy.fromBoundaryRows(boundaryRows, expectedZeroCount=1, validationMode="warning"), ...
                parameters=struct("problemType", "hydrostaticFModes", "omega", NaN, ...
                "upperBoundary", upperRow.family, "lowerBoundary", lowerRow.family));
            evp.components.F.innerWeight = @(z,ctx) ones(size(z));
            evp.components.G.role = "diagnostic";
            evp.components.G.from = "F";
            evp.components.G.operator = IMOperator.strong().plus(coefficient=@(z,ctx) -ctx.g./ctx.N2(z), derivativeOrder=1);
            evp.components.G.innerWeight = @(z,ctx) ctx.N2(z)/ctx.g;
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

            boundaryRows = IMBoundaryRow.partialDepthPE(boundarySign=options.boundarySign);
            policy = IMIndexPolicy.fromBoundaryRows(boundaryRows, validationMode=options.validationMode);
        end
    end

end

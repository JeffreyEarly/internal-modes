classdef IMEigenvalueProblem
    % Describe an internal-mode eigenvalue problem in physical coordinates.
    %
    % `IMEigenvalueProblem` stores the physical operators, boundary laws,
    % component metadata, named normalization rules, eigenvalue
    % interpretation, and index policy needed by coordinate-aware solvers.
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

        % Physical parameters used to create the EVP.
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
            % - Parameter options.leftOperator: left physical operator
            % - Parameter options.rightOperator: right physical operator
            % - Parameter options.normalizations: named modal normalization rules
            % - Parameter options.boundaryRows: resolved boundary rows
            % - Parameter options.eigenvalueName: eigenvalue name
            % - Parameter options.hFromEigenvalue: equivalent-depth conversion
            % - Parameter options.ordering: eigenvalue ordering rule
            % - Parameter options.indexPolicy: expected index policy
            % - Parameter options.parameters: physical parameters used to create the EVP
            % - Returns evp: initialized EVP descriptor
            arguments
                options.name {mustBeTextScalar} = "unknown"
                options.primaryComponent {mustBeTextScalar} = "G"
                options.leftOperator IMOperator = IMOperator.strong()
                options.rightOperator IMOperator = IMOperator.strong()
                options.normalizations struct = struct()
                options.boundaryRows = IMBoundaryRow.empty(0,1)
                options.eigenvalueName {mustBeTextScalar} = "lambda"
                options.hFromEigenvalue = @(lambda) 1 ./ lambda
                options.ordering {mustBeTextScalar} = "ascendingEigenvalue"
                options.indexPolicy IMIndexPolicy = IMIndexPolicy.none()
                options.parameters struct = struct()
            end

            self.name = string(options.name);
            self.primaryComponent = string(options.primaryComponent);
            self.leftOperator = options.leftOperator;
            self.rightOperator = options.rightOperator;
            self.normalizations = options.normalizations;
            self.boundaryRows = options.boundaryRows(:);
            self.eigenvalueName = string(options.eigenvalueName);
            self.hFromEigenvalue = options.hFromEigenvalue;
            self.ordering = string(options.ordering);
            self.indexPolicy = options.indexPolicy;
            self.parameters = options.parameters;
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
            A = self.leftOperator.matrix(solver);
            B = self.rightOperator.matrix(solver);
            for iBoundary = 1:length(self.boundaryRows)
                [A, B] = self.boundaryRows(iBoundary).apply(A, B, solver);
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
            right = IMOperator.strong().plus(coefficient=@(z,ctx) (f0*f0 - ctx.N2(z))/g, derivativeOrder=0);
            context = struct("problemType", "waveModesAtWavenumber", "primaryComponent", "G", ...
                "k", k, "omega", NaN, "f0", f0, "g", g);
            lowerRow = options.lowerBoundary.resolve(endpoint="lower", context=context);
            upperRow = options.upperBoundary.resolve(endpoint="upper", context=context);
            boundaryRows = [lowerRow; upperRow];
            parameters = struct("problemType", "waveModesAtWavenumber", "k", k, "omega", NaN, ...
                "f0", f0, "g", g, ...
                "upperBoundary", upperRow.family, "lowerBoundary", lowerRow.family);
            evp = IMEigenvalueProblem(name="waveModesAtWavenumber", primaryComponent="G", ...
                leftOperator=left, rightOperator=right, boundaryRows=boundaryRows, ...
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
            right = IMOperator.strong().plus(coefficient=@(z,ctx) (omega*omega - ctx.N2(z))/g, derivativeOrder=0);
            context = struct("problemType", "waveModesAtFrequency", "primaryComponent", "G", ...
                "k", NaN, "omega", omega, "f0", f0, "g", g);
            lowerRow = options.lowerBoundary.resolve(endpoint="lower", context=context);
            upperRow = options.upperBoundary.resolve(endpoint="upper", context=context);
            boundaryRows = [lowerRow; upperRow];
            parameters = struct("problemType", "waveModesAtFrequency", "k", NaN, "omega", omega, ...
                "f0", f0, "g", g, ...
                "upperBoundary", upperRow.family, "lowerBoundary", lowerRow.family);
            evp = IMEigenvalueProblem(name="waveModesAtFrequency", primaryComponent="G", ...
                leftOperator=left, rightOperator=right, boundaryRows=boundaryRows, ...
                eigenvalueName="inverseEquivalentDepth", hFromEigenvalue=@(lambda) 1 ./ lambda, ...
                ordering="ascendingEigenvalue", ...
                indexPolicy=IMIndexPolicy.fromBoundaryRows(boundaryRows, validationMode="none"), parameters=parameters);
            evp.components.G.innerWeight = @(z,ctx) (ctx.N2(z) - ctx.omega*ctx.omega)/ctx.g;
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
            % Hydrostatic `G` modes are the zero-frequency limit of the
            % fixed-frequency wave-mode EVP.
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

            evp = IMEigenvalueProblem.waveModesAtFrequency(omega=0, f0=options.f0, g=options.g, upperBoundary=options.upperBoundary, lowerBoundary=options.lowerBoundary);
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
            right = IMOperator.strong().plus(coefficient=@(z,ctx) -ctx.N2(z)/g, derivativeOrder=0);
            context = struct("problemType", "hydrostaticFModes", "primaryComponent", "F", "omega", NaN, "g", g);
            lowerRow = options.lowerBoundary.resolve(endpoint="lower", context=context);
            upperRow = options.upperBoundary.resolve(endpoint="upper", context=context);
            boundaryRows = [lowerRow; upperRow];
            evp = IMEigenvalueProblem(name="hydrostaticFModes", primaryComponent="F", ...
                leftOperator=left, rightOperator=right, boundaryRows=boundaryRows, ...
                eigenvalueName="inverseEquivalentDepth", hFromEigenvalue=@(lambda) 1 ./ lambda, ...
                ordering="indexThenAscending", ...
                indexPolicy=IMIndexPolicy.fromBoundaryRows(boundaryRows, expectedZeroCount=1, validationMode="warning"), ...
                parameters=struct("problemType", "hydrostaticFModes", "omega", NaN, "g", g, ...
                "upperBoundary", upperRow.family, "lowerBoundary", lowerRow.family));
            evp.components.F.innerWeight = @(z,ctx) ones(size(z))/diff(ctx.zDomain);
            evp.components.G.role = "diagnostic";
            evp.components.G.from = "F";
            evp.components.G.operator = IMOperator.strong().plus( coefficient=@(z,ctx) -ctx.g./ctx.N2(z), derivativeOrder=1);
            evp.normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor("F", iMode);
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

classdef InternalModesEVP
    % Describe an internal-mode eigenvalue problem in physical coordinates.
    %
    % `InternalModesEVP` stores the physical operators, boundary conditions,
    % component metadata, eigenvalue interpretation, and index policy needed
    % by coordinate-aware solvers.
    %
    % ```matlab
    % evp = InternalModesEVP.hydrostaticGModes(k=1e-4, f0=1e-4);
    % ```
    %
    % - Topic: Create EVPs
    % - Topic: Assemble EVPs
    % - Topic: Inspect EVP metadata
    % - Declaration: classdef InternalModesEVP

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
        leftOperator = InternalModesOperator.strong()

        % Right physical operator.
        %
        % - Topic: Assemble EVPs
        rightOperator = InternalModesOperator.strong()

        % Component metadata.
        %
        % - Topic: Inspect EVP metadata
        components = struct()

        % Boundary-condition array.
        %
        % - Topic: Assemble EVPs
        boundaryConditions = InternalModesBoundaryCondition.empty(0,1)

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
        indexPolicy = InternalModesIndexPolicy.none()
    end

    methods
        function self = InternalModesEVP(options)
            % Create a physical-coordinate EVP descriptor.
            %
            % - Topic: Create EVPs
            % - Declaration: evp = InternalModesEVP(options)
            % - Parameter options.name: short EVP name
            % - Parameter options.primaryComponent: primary eigenfunction component
            % - Parameter options.leftOperator: left physical operator
            % - Parameter options.rightOperator: right physical operator
            % - Parameter options.boundaryConditions: boundary conditions
            % - Parameter options.eigenvalueName: eigenvalue name
            % - Parameter options.hFromEigenvalue: equivalent-depth conversion
            % - Parameter options.ordering: eigenvalue ordering rule
            % - Parameter options.indexPolicy: expected index policy
            % - Returns evp: initialized EVP descriptor
            arguments
                options.name {mustBeTextScalar} = "unknown"
                options.primaryComponent {mustBeTextScalar} = "G"
                options.leftOperator InternalModesOperator = InternalModesOperator.strong()
                options.rightOperator InternalModesOperator = InternalModesOperator.strong()
                options.boundaryConditions = InternalModesBoundaryCondition.empty(0,1)
                options.eigenvalueName {mustBeTextScalar} = "lambda"
                options.hFromEigenvalue = @(lambda) 1 ./ lambda
                options.ordering {mustBeTextScalar} = "ascendingEigenvalue"
                options.indexPolicy InternalModesIndexPolicy = InternalModesIndexPolicy.none()
            end

            self.name = string(options.name);
            self.primaryComponent = string(options.primaryComponent);
            self.leftOperator = options.leftOperator;
            self.rightOperator = options.rightOperator;
            self.boundaryConditions = options.boundaryConditions(:);
            self.eigenvalueName = string(options.eigenvalueName);
            self.hFromEigenvalue = options.hFromEigenvalue;
            self.ordering = string(options.ordering);
            self.indexPolicy = options.indexPolicy;
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
            for iCondition = 1:length(self.boundaryConditions)
                [A, B] = self.boundaryConditions(iCondition).apply(A, B, solver);
            end
        end
    end

    methods (Static)
        function evp = hydrostaticGModes(options)
            % Create the hydrostatic `G`-mode EVP at fixed horizontal wavenumber.
            %
            % The physical-coordinate strong form is
            % $$G_{zz}-K^2G=\lambda(f_0^2-N^2)G/g$$.
            %
            % - Topic: Create EVPs
            % - Declaration: evp = InternalModesEVP.hydrostaticGModes(options)
            % - Parameter options.k: horizontal wavenumber
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.upperBoundary: upper boundary condition
            % - Parameter options.lowerBoundary: lower boundary condition
            % - Returns evp: hydrostatic `G` EVP
            arguments
                options.k (1,1) double {mustBeNonnegative}
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.upperBoundary {mustBeTextScalar} = "rigidLid"
                options.lowerBoundary {mustBeTextScalar} = "rigid"
            end

            k = options.k;
            f0 = options.f0;
            g = options.g;
            left = InternalModesOperator.strong().plus(derivativeOrder=2).plus(coefficient=-k*k, derivativeOrder=0);
            right = InternalModesOperator.strong().plus(coefficient=@(z,ctx) (f0*f0 - ctx.N2(z))/g, derivativeOrder=0);
            boundaryConditions = InternalModesEVP.hydrostaticGBoundaries(options.upperBoundary, options.lowerBoundary);
            evp = InternalModesEVP(name="hydrostaticGModes", primaryComponent="G", ...
                leftOperator=left, rightOperator=right, boundaryConditions=boundaryConditions, ...
                eigenvalueName="inverseEquivalentDepth", hFromEigenvalue=@(lambda) 1 ./ lambda, ...
                ordering="ascendingEigenvalue", indexPolicy=InternalModesIndexPolicy.none());
            evp.components.F.role = "diagnostic";
            evp.components.F.from = "G";
            evp.components.F.operator = InternalModesOperator.strong().plus(derivativeOrder=1);
        end

        function evp = hydrostaticFModes(options)
            % Create the geostrophic hydrostatic `F`-mode EVP.
            %
            % The physical-coordinate strong form is
            % $$F_{zz}-(\partial_z\log N^2)F_z=-\lambda N^2F/g$$.
            %
            % - Topic: Create EVPs
            % - Declaration: evp = InternalModesEVP.hydrostaticFModes(options)
            % - Parameter options.k: horizontal wavenumber
            % - Parameter options.g: gravitational acceleration
            % - Returns evp: hydrostatic `F` EVP
            arguments
                options.k (1,1) double {mustBeNonnegative} = 0
                options.g (1,1) double {mustBePositive} = 9.81
            end

            g = options.g;
            left = InternalModesOperator.strong().plus(derivativeOrder=2) ...
                .plus(coefficient=@(z,ctx) -ctx.dzLogN2(z), derivativeOrder=1);
            right = InternalModesOperator.strong().plus(coefficient=@(z,ctx) -ctx.N2(z)/g, derivativeOrder=0);
            boundaryConditions = [
                InternalModesBoundaryCondition.neumann(location="bottom", component="F")
                InternalModesBoundaryCondition.neumann(location="surface", component="F")
            ];
            evp = InternalModesEVP(name="hydrostaticFModes", primaryComponent="F", ...
                leftOperator=left, rightOperator=right, boundaryConditions=boundaryConditions, ...
                eigenvalueName="inverseEquivalentDepth", hFromEigenvalue=@(lambda) 1 ./ lambda, ...
                ordering="indexThenAscending", indexPolicy=InternalModesIndexPolicy.fixed(expectedZeroCount=1, validationMode="warning"));
            evp.components.G.role = "diagnostic";
            evp.components.G.from = "F";
        end

        function policy = partialDepthPEIndexPolicy(options)
            % Return the manuscript partial-depth PE index policy.
            %
            % - Topic: Create EVPs
            % - Declaration: policy = InternalModesEVP.partialDepthPEIndexPolicy(options)
            % - Parameter options.boundarySign: `"positive"` or `"negative"`
            % - Parameter options.validationMode: `"error"`, `"warning"`, or `"none"`
            % - Returns policy: corresponding index policy
            arguments
                options.boundarySign {mustBeTextScalar} = "positive"
                options.validationMode {mustBeTextScalar} = "error"
            end

            switch string(options.boundarySign)
                case "positive"
                    signs = [1; 1];
                case "negative"
                    signs = [-1; -1];
                otherwise
                    error("InternalModesEVP:InvalidBoundarySign", ...
                        "boundarySign must be ""positive"" or ""negative"".");
            end
            policy = InternalModesIndexPolicy.fromBoundarySigns(signs, validationMode=options.validationMode);
        end
    end

    methods (Static, Access = private)
        function boundaryConditions = hydrostaticGBoundaries(upperBoundary, lowerBoundary)
            lowerBoundary = string(lowerBoundary);
            upperBoundary = string(upperBoundary);
            switch lowerBoundary
                case {"rigid", "rigidLid", "freeSlip"}
                    lower = InternalModesBoundaryCondition.dirichlet(location="bottom", component="G");
                otherwise
                    error("InternalModesEVP:UnsupportedBoundary", ...
                        "Unsupported hydrostatic G lower boundary ""%s"".", lowerBoundary);
            end

            switch upperBoundary
                case {"rigid", "rigidLid"}
                    upper = InternalModesBoundaryCondition.dirichlet(location="surface", component="G");
                case "freeSurface"
                    left = InternalModesOperator.strong().plus(derivativeOrder=1);
                    right = InternalModesOperator.strong().plus(derivativeOrder=0);
                    upper = InternalModesBoundaryCondition.robin(location="surface", component="G", ...
                        operator=left, eigenvalueOperator=right);
                otherwise
                    error("InternalModesEVP:UnsupportedBoundary", ...
                        "Unsupported hydrostatic G upper boundary ""%s"".", upperBoundary);
            end

            boundaryConditions = [lower; upper];
        end
    end
end

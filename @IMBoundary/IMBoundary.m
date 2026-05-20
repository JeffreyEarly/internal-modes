classdef IMBoundary
    % Describe a location-free boundary law for a v2 internal-mode EVP.
    %
    % `IMBoundary` stores the physical law a user wants to impose. EVP
    % factories place that law at the upper or lower endpoint by calling
    % `resolve`, which returns an `IMBoundaryRow` with explicit row
    % operators, endpoint trace terms, and index metadata.
    %
    % ```matlab
    % evp = IMEigenvalueProblem.waveModesAtWavenumber( ...
    %     k=1e-4, upperBoundary=IMBoundary.free(), lowerBoundary=IMBoundary.rigid());
    % ```
    %
    % - Topic: Create boundary laws
    % - Topic: Resolve boundary laws
    % - Topic: Inspect boundary laws
    % - Topic: Developer topics
    % - Declaration: classdef IMBoundary

    properties
        % Short boundary-law name.
        %
        % - Topic: Inspect boundary laws
        name = "boundary"

        % Boundary family name.
        %
        % Examples include `"rigid"`, `"noSlip"`, `"free"`, `"linearF"`,
        % and `"linearG"`.
        %
        % - Topic: Inspect boundary laws
        family = "custom"

        % Component used by explicit developer-level laws.
        %
        % Named physical laws such as `rigid`, `noSlip`, and `free` choose
        % their component when resolved for a `G` or `F` EVP.
        %
        % - Topic: Inspect boundary laws
        component = ""

        % Left-side boundary functional for explicit laws.
        %
        % - Topic: Developer topics
        % - Developer: true
        leftOperator = IMOperator.strong()

        % Right-side eigenvalue boundary functional for explicit laws.
        %
        % - Topic: Developer topics
        % - Developer: true
        rightOperator = IMOperator.strong()

        % Endpoint trace terms for explicit or linear laws before placement.
        %
        % Terms whose location is empty are placed at the resolved endpoint.
        %
        % - Topic: Developer topics
        % - Developer: true
        endpointTerms = IMBoundaryRow.emptyEndpointTerms()

        % Status of the associated orthogonality metadata.
        %
        % Values are `"complete"` or `"unresolved"`. Unresolved linear
        % families can still assemble EVP rows, but their endpoint inner
        % products and index metadata are incomplete.
        %
        % - Topic: Inspect boundary laws
        orthogonalityStatus = "complete"

        % Family-specific coefficients for location-free laws.
        %
        % `linearF` stores fields `a`, `b`, `c`, and `d`. `linearG` stores
        % fields `a`, `b`, `c`, and `e`.
        %
        % - Topic: Developer topics
        % - Developer: true
        coefficients = struct()
    end

    methods
        function self = IMBoundary(options)
            % Create a location-free boundary law.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary(options)
            % - Parameter options.name: short boundary-law name
            % - Parameter options.family: boundary family name
            % - Parameter options.component: constrained component for explicit laws
            % - Parameter options.leftOperator: left-side boundary functional
            % - Parameter options.rightOperator: right-side eigenvalue functional
            % - Parameter options.endpointTerms: endpoint trace-pair terms
            % - Parameter options.orthogonalityStatus: endpoint-form status
            % - Parameter options.coefficients: family-specific coefficient metadata
            % - Returns boundary: initialized boundary law
            arguments
                options.name {mustBeTextScalar} = "boundary"
                options.family {mustBeTextScalar} = "custom"
                options.component {mustBeTextScalar} = ""
                options.leftOperator IMOperator = IMOperator.strong()
                options.rightOperator IMOperator = IMOperator.strong()
                options.endpointTerms struct = IMBoundaryRow.emptyEndpointTerms()
                options.orthogonalityStatus {mustBeTextScalar} = "complete"
                options.coefficients struct = struct()
            end

            self.name = string(options.name);
            self.family = string(options.family);
            self.component = string(options.component);
            self.leftOperator = options.leftOperator;
            self.rightOperator = options.rightOperator;
            self.endpointTerms = options.endpointTerms(:);
            self.orthogonalityStatus = string(options.orthogonalityStatus);
            self.coefficients = options.coefficients;
        end

        function row = resolve(self, options)
            % Resolve this boundary law at one endpoint.
            %
            % The returned `IMBoundaryRow` is the explicit object consumed
            % by EVP assembly, basis-set endpoint inner products, and index
            % policies. `context.primaryComponent` determines how named
            % physical laws resolve for `G` and `F` EVPs.
            %
            % - Topic: Resolve boundary laws
            % - Declaration: row = resolve(boundary,options)
            % - Parameter options.endpoint: `"upper"` or `"lower"`
            % - Parameter options.context: EVP construction context
            % - Returns row: resolved boundary row
            arguments
                self IMBoundary
                options.endpoint {mustBeTextScalar}
                options.context struct = struct()
            end

            endpoint = string(options.endpoint);
            location = IMBoundaryRow.locationForEndpoint(endpoint);
            primaryComponent = IMBoundary.contextString(options.context, "primaryComponent", "G");

            switch self.family
                case "rigid"
                    switch primaryComponent
                        case "G"
                            left = IMOperator.strong().plus(derivativeOrder=0);
                            row = IMBoundaryRow(family="rigid", endpoint=endpoint, component="G", leftOperator=left);
                        case "F"
                            left = IMOperator.strong().plus(derivativeOrder=1);
                            row = IMBoundaryRow(family="rigid", endpoint=endpoint, component="F", leftOperator=left);
                        otherwise
                            IMBoundary.unsupportedResolution(self.family, primaryComponent);
                    end
                case "noSlip"
                    switch primaryComponent
                        case "G"
                            left = IMOperator.strong().plus(derivativeOrder=1);
                            row = IMBoundaryRow(family="noSlip", endpoint=endpoint, component="G", leftOperator=left);
                        case "F"
                            left = IMOperator.strong().plus(derivativeOrder=0);
                            row = IMBoundaryRow(family="noSlip", endpoint=endpoint, component="F", leftOperator=left);
                        otherwise
                            IMBoundary.unsupportedResolution(self.family, primaryComponent);
                    end
                case "free"
                    endpointTerm = IMBoundaryRow.endpointTerm("G", location, 1, ...
                        IMBoundaryRow.trace("G"), IMBoundaryRow.trace("G"));
                    switch primaryComponent
                        case "G"
                            left = IMOperator.strong().plus(derivativeOrder=1);
                            right = IMOperator.strong().plus(derivativeOrder=0);
                            row = IMBoundaryRow(family="free", endpoint=endpoint, component="G", ...
                                leftOperator=left, rightOperator=right, endpointTerms=endpointTerm);
                        case "F"
                            left = IMOperator.strong() ...
                                .plus(derivativeOrder=0) ...
                                .plus(coefficient=@(z,ctx) ctx.g./ctx.N2(z), derivativeOrder=1);
                            row = IMBoundaryRow(family="free", endpoint=endpoint, component="F", ...
                                leftOperator=left, endpointTerms=endpointTerm);
                        otherwise
                            IMBoundary.unsupportedResolution(self.family, primaryComponent);
                    end
                case "dirichlet"
                    left = IMOperator.strong().plus(derivativeOrder=0);
                    row = IMBoundaryRow(family="dirichlet", endpoint=endpoint, component=self.component, leftOperator=left);
                case "neumann"
                    left = IMOperator.strong().plus(derivativeOrder=1);
                    row = IMBoundaryRow(family="neumann", endpoint=endpoint, component=self.component, leftOperator=left);
                case "robin"
                    endpointTerms = IMBoundary.resolveEndpointLocations(self.endpointTerms, location);
                    row = IMBoundaryRow(family="robin", endpoint=endpoint, component=self.component, ...
                        leftOperator=self.leftOperator, rightOperator=self.rightOperator, endpointTerms=endpointTerms);
                case "linearF"
                    if primaryComponent ~= "F"
                        IMBoundary.unsupportedResolution(self.family, primaryComponent);
                    end
                    c = self.coefficients;
                    left = IMOperator.strong() ...
                        .plus(coefficient=-c.a, derivativeOrder=0) ...
                        .plus(coefficient=@(z,ctx) c.b./ctx.N2(z), derivativeOrder=1);
                    right = IMOperator.strong() ...
                        .plus(coefficient=@(~,ctx) c.c/ctx.g, derivativeOrder=0) ...
                        .plus(coefficient=@(z,ctx) -c.d./(ctx.g*ctx.N2(z)), derivativeOrder=1);
                    endpointTerms = IMBoundary.resolveEndpointLocations(self.endpointTerms, location);
                    row = IMBoundaryRow(family="linearF", endpoint=endpoint, component="F", ...
                        leftOperator=left, rightOperator=right, endpointTerms=endpointTerms, ...
                        orthogonalityStatus=self.orthogonalityStatus);
                case "linearG"
                    if primaryComponent ~= "G"
                        IMBoundary.unsupportedResolution(self.family, primaryComponent);
                    end
                    c = self.coefficients;
                    left = IMOperator.strong() ...
                        .plus(coefficient=-c.e, derivativeOrder=0) ...
                        .plus(coefficient=c.a, derivativeOrder=1);
                    right = IMOperator.strong() ...
                        .plus(coefficient=@(~,ctx) c.b/ctx.g, derivativeOrder=0) ...
                        .plus(coefficient=@(~,ctx) -c.c/ctx.g, derivativeOrder=1);
                    endpointTerms = IMBoundary.resolveEndpointLocations(self.endpointTerms, location);
                    row = IMBoundaryRow(family="linearG", endpoint=endpoint, component="G", ...
                        leftOperator=left, rightOperator=right, endpointTerms=endpointTerms, ...
                        orthogonalityStatus=self.orthogonalityStatus);
                case "custom"
                    endpointTerms = IMBoundary.resolveEndpointLocations(self.endpointTerms, location);
                    row = IMBoundaryRow(family="custom", endpoint=endpoint, component=self.component, ...
                        leftOperator=self.leftOperator, rightOperator=self.rightOperator, endpointTerms=endpointTerms, ...
                        orthogonalityStatus=self.orthogonalityStatus);
                otherwise
                    error("IMBoundary:UnsupportedResolution", ...
                        "Boundary family ""%s"" cannot be resolved by the v2 boundary API.", self.family);
            end
        end
    end

    methods (Static)
        function boundary = dirichlet(options)
            % Create a location-free homogeneous Dirichlet boundary law.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.dirichlet(options)
            % - Parameter options.component: constrained component
            % - Returns boundary: initialized boundary law
            arguments
                options.component {mustBeTextScalar} = "G"
            end

            boundary = IMBoundary(name=string(options.component) + "Dirichlet", family="dirichlet", ...
                component=options.component);
        end

        function boundary = neumann(options)
            % Create a location-free homogeneous Neumann boundary law.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.neumann(options)
            % - Parameter options.component: constrained component
            % - Returns boundary: initialized boundary law
            arguments
                options.component {mustBeTextScalar} = "G"
            end

            boundary = IMBoundary(name=string(options.component) + "Neumann", family="neumann", ...
                component=options.component);
        end

        function boundary = robin(options)
            % Create a location-free homogeneous Robin boundary law.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.robin(options)
            % - Parameter options.component: constrained component
            % - Parameter options.leftOperator: left-side boundary functional
            % - Parameter options.rightOperator: right-side eigenvalue functional
            % - Parameter options.endpointTerms: endpoint trace-pair terms
            % - Returns boundary: initialized boundary law
            arguments
                options.component {mustBeTextScalar} = "G"
                options.leftOperator IMOperator
                options.rightOperator IMOperator = IMOperator.strong()
                options.endpointTerms struct = IMBoundaryRow.emptyEndpointTerms()
            end

            boundary = IMBoundary(name=string(options.component) + "Robin", family="robin", ...
                component=options.component, leftOperator=options.leftOperator, ...
                rightOperator=options.rightOperator, endpointTerms=options.endpointTerms);
        end

        function boundary = rigid()
            % Create a location-free rigid boundary law.
            %
            % In a `G` EVP, `rigid` resolves to $$G=0$$. In an `F` EVP, it
            % resolves to $$F_z=0$$.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.rigid()
            % - Returns boundary: initialized rigid boundary law
            boundary = IMBoundary(name="rigid", family="rigid");
        end

        function boundary = noSlip()
            % Create a location-free no-slip boundary law.
            %
            % In a `G` EVP, `noSlip` resolves to $$G_z=0$$. In an `F` EVP,
            % it resolves to $$F=0$$.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.noSlip()
            % - Returns boundary: initialized no-slip boundary law
            boundary = IMBoundary(name="noSlip", family="noSlip");
        end

        function boundary = free()
            % Create a location-free free boundary law.
            %
            % In a `G` EVP, `free` resolves to $$G_z=\lambda G$$. In an `F`
            % EVP, it resolves to $$F+gF_z/N^2=0$$. The law may be placed at
            % either endpoint by an EVP factory.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.free()
            % - Returns boundary: initialized free boundary law
            boundary = IMBoundary(name="free", family="free");
        end

        function boundary = linearF(options)
            % Create a location-free linear `F` boundary law.
            %
            % When resolved, the assembled boundary row represents
            % $$-(aF-bF_z/N^2)=\lambda(cF-dF_z/N^2)/g$$. Trusted pairwise
            % cases include endpoint terms; unresolved cases warn and leave
            % orthogonality metadata incomplete.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.linearF(options)
            % - Parameter options.a: coefficient multiplying `F`
            % - Parameter options.b: coefficient multiplying `F_z/N2`
            % - Parameter options.c: eigenvalue coefficient multiplying `F`
            % - Parameter options.d: eigenvalue coefficient multiplying `F_z/N2`
            % - Returns boundary: initialized linear `F` boundary law
            arguments
                options.a (1,1) double = 0
                options.b (1,1) double = 0
                options.c (1,1) double = 0
                options.d (1,1) double = 0
            end

            [endpointTerms, status] = IMBoundary.linearFEndpointTerms("", options.a, options.b, options.c, options.d);
            coefficients = struct("a", options.a, "b", options.b, "c", options.c, "d", options.d);
            boundary = IMBoundary(name="linearF", family="linearF", component="F", ...
                endpointTerms=endpointTerms, orthogonalityStatus=status, coefficients=coefficients);
        end

        function boundary = linearG(options)
            % Create a location-free linear `G` boundary law.
            %
            % When resolved, the assembled boundary row represents
            % $$-(eG-aG_z)=\lambda(bG-cG_z)/g$$. Trusted pairwise cases
            % include endpoint terms; unresolved cases warn and leave
            % orthogonality metadata incomplete.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.linearG(options)
            % - Parameter options.a: coefficient multiplying `G_z`
            % - Parameter options.b: eigenvalue coefficient multiplying `G`
            % - Parameter options.c: eigenvalue coefficient multiplying `G_z`
            % - Parameter options.e: coefficient multiplying `G`
            % - Returns boundary: initialized linear `G` boundary law
            arguments
                options.a (1,1) double = 0
                options.b (1,1) double = 0
                options.c (1,1) double = 0
                options.e (1,1) double = 0
            end

            [endpointTerms, status] = IMBoundary.linearGEndpointTerms("", options.a, options.b, options.c, options.e);
            coefficients = struct("a", options.a, "b", options.b, "c", options.c, "e", options.e);
            boundary = IMBoundary(name="linearG", family="linearG", component="G", ...
                endpointTerms=endpointTerms, orthogonalityStatus=status, coefficients=coefficients);
        end
    end

    methods (Static, Access = private)
        function [terms, status] = linearFEndpointTerms(location, a, b, c, d)
            nonzero = abs([a b c d]) > 0;
            terms = IMBoundaryRow.emptyEndpointTerms();
            status = "complete";
            if nnz(nonzero) > 2
                IMBoundary.warnUnresolved("linearF");
                status = "unresolved";
                return;
            end

            g = @(ctx) ctx.g;
            if ~nonzero(1) && ~nonzero(2) && nonzero(3) && nonzero(4)
                terms = IMBoundaryRow.endpointTerm("G", location, -d/c, IMBoundaryRow.trace("G"), IMBoundaryRow.trace("G"));
            elseif nonzero(1) && nonzero(2) && ~nonzero(3) && ~nonzero(4)
                terms = IMBoundaryRow.endpointTerm("G", location, @(ctx) -b/(g(ctx)*a), ...
                    IMBoundaryRow.trace("G"), IMBoundaryRow.trace("G"));
            elseif ~nonzero(1) && nonzero(2) && nonzero(3) && ~nonzero(4)
                terms = IMBoundaryRow.endpointTerm("F", location, c/b, IMBoundaryRow.trace("F"), IMBoundaryRow.trace("F"));
            elseif nonzero(1) && ~nonzero(2) && ~nonzero(3) && nonzero(4)
                terms = [
                    IMBoundaryRow.endpointTerm("F", location, @(ctx) -d/(g(ctx)*g(ctx)*a), ...
                        IMBoundaryRow.trace("G"), IMBoundaryRow.trace("G"))
                    IMBoundaryRow.endpointTerm("G", location, 1, IMBoundaryRow.trace("G"), IMBoundaryRow.trace("F"))
                    IMBoundaryRow.endpointTerm("G", location, 1, IMBoundaryRow.trace("F"), IMBoundaryRow.trace("G"))
                ];
            elseif nnz(nonzero) <= 2
            else
                IMBoundary.warnUnresolved("linearF");
                status = "unresolved";
            end
        end

        function [terms, status] = linearGEndpointTerms(location, a, b, c, e)
            nonzero = abs([a b c e]) > 0;
            terms = IMBoundaryRow.emptyEndpointTerms();
            status = "complete";
            if nnz(nonzero) > 2
                IMBoundary.warnUnresolved("linearG");
                status = "unresolved";
                return;
            end

            if ~nonzero(1) && ~nonzero(2) && nonzero(3) && nonzero(4)
                terms = IMBoundaryRow.endpointTerm("G", location, -c/e, ...
                    IMBoundaryRow.trace("G", derivativeOrder=1), IMBoundaryRow.trace("G", derivativeOrder=1));
            elseif nonzero(1) && nonzero(2) && ~nonzero(3) && ~nonzero(4)
                terms = IMBoundaryRow.endpointTerm("G", location, b/a, IMBoundaryRow.trace("G"), IMBoundaryRow.trace("G"));
            elseif nonzero(1) && ~nonzero(2) && ~nonzero(3) && nonzero(4)
                terms = IMBoundaryRow.endpointTerm("F", location, -a/e, IMBoundaryRow.trace("F"), IMBoundaryRow.trace("F"));
            elseif ~nonzero(1) && nonzero(2) && nonzero(3) && ~nonzero(4)
                terms = IMBoundaryRow.endpointTerm("F", location, -c/b, IMBoundaryRow.trace("F"), IMBoundaryRow.trace("F"));
            elseif nnz(nonzero) <= 2
            else
                IMBoundary.warnUnresolved("linearG");
                status = "unresolved";
            end
        end

        function terms = resolveEndpointLocations(terms, location)
            for iTerm = 1:length(terms)
                if string(terms(iTerm).location) == ""
                    terms(iTerm).location = string(location);
                end
            end
        end

        function value = contextString(context, name, defaultValue)
            if isfield(context, name)
                value = string(context.(name));
            else
                value = string(defaultValue);
            end
        end

        function unsupportedResolution(family, primaryComponent)
            error("IMBoundary:UnsupportedResolution", ...
                "Boundary family ""%s"" cannot resolve for primary component ""%s"".", ...
                string(family), string(primaryComponent));
        end

        function warnUnresolved(family)
            warning("IMBoundary:UnsupportedLinearFamily", ...
                "%s can assemble this boundary row, but its orthogonality and index metadata are unresolved.", family);
        end
    end
end

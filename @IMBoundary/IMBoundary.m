classdef IMBoundary
    % Describe a location-free boundary law for a v2 internal-mode EVP.
    %
    % `IMBoundary` stores the physical law a user wants to impose, not the
    % placed matrix row. Users normally pass boundary laws to EVP factories:
    %
    % ```matlab
    % evp = IMEigenvalueProblem.waveModesAtWavenumber( ...
    %     k=1e-4, upperBoundary=IMBoundary.free(), lowerBoundary=IMBoundary.rigid());
    % ```
    %
    % The factory places the law at the upper or lower endpoint and resolves
    % it for the EVP primary component. For example, `IMBoundary.dirichlet()`
    % acts on the primary component, while `IMBoundary.dirichlet(on="F")`
    % is an explicit advanced override. The resolved developer object is an
    % `IMBoundaryRow`, which owns the matrix row, boundary inner-product
    % terms, and index metadata.
    %
    % Boundary inner-product terms are not boundary conditions. They are the
    % boundary trace products added to modal inner products so the EVP
    % remains orthogonal under the chosen boundary law. Examples include
    % scalar products like $$G_iG_j$$ and mixed products like
    % $$G_iF_j+F_iG_j$$.
    %
    % - Topic: Create boundary laws
    % - Topic: Resolve boundary laws
    % - Topic: Inspect boundary laws
    % - Topic: Developer topics
    % - Declaration: classdef IMBoundary

    properties (SetAccess = private)
        % Boundary family identifier.
        %
        % Examples include `"rigid"`, `"noSlip"`, `"free"`, `"linearF"`,
        % and `"linearG"`. This is the stable public identifier recorded in
        % EVP metadata and resolved boundary rows.
        %
        % - Topic: Inspect boundary laws
        family = "operator"
    end

    properties (Access = private)
        on = "primary"
        leftOperator = IMOperator.strong()
        rightOperator = IMOperator.strong()
        innerProductTerms = IMBoundaryRow.emptyInnerProductTerms()
        hasKnownInnerProduct = true
        coefficients = struct()
    end

    methods (Hidden)
        function row = resolve(self, options)
            % Resolve this boundary law at one endpoint.
            %
            % This is a developer-facing method used by EVP factories and
            % tests. `endpoint` determines the placement, and
            % `primaryComponent` determines how `"primary"` component laws
            % such as `dirichlet()` and `neumann()` are interpreted.
            %
            % ```matlab
            % row = IMBoundary.dirichlet().resolve(endpoint="upper", primaryComponent="G");
            % row.component  % "G"
            % ```
            %
            % - Topic: Resolve boundary laws
            % - Developer: true
            % - Declaration: row = resolve(boundary,options)
            % - Parameter options.endpoint: `"upper"` or `"lower"`
            % - Parameter options.primaryComponent: primary EVP component, `"G"` or `"F"`
            % - Returns row: resolved boundary row
            arguments
                self IMBoundary
                options.endpoint {mustBeTextScalar}
                options.primaryComponent {mustBeTextScalar} = "G"
            end

            endpoint = string(options.endpoint);
            location = IMBoundaryRow.locationForEndpoint(endpoint);
            primaryComponent = IMBoundary.validateComponent(options.primaryComponent);

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
                    innerProductTerm = IMBoundaryRow.innerProductTerm("G", location, 1, ...
                        IMBoundaryRow.trace("G"), IMBoundaryRow.trace("G"));
                    switch primaryComponent
                        case "G"
                            left = IMOperator.strong().plus(derivativeOrder=1);
                            right = IMOperator.strong().plus(derivativeOrder=0);
                            row = IMBoundaryRow(family="free", endpoint=endpoint, component="G", ...
                                leftOperator=left, rightOperator=right, innerProductTerms=innerProductTerm);
                        case "F"
                            left = IMOperator.strong() ...
                                .plus(derivativeOrder=0) ...
                                .plus(coefficient=@(z,ctx) ctx.g./ctx.N2(z), derivativeOrder=1);
                            row = IMBoundaryRow(family="free", endpoint=endpoint, component="F", ...
                                leftOperator=left, innerProductTerms=innerProductTerm);
                        otherwise
                            IMBoundary.unsupportedResolution(self.family, primaryComponent);
                    end
                case "dirichlet"
                    left = IMOperator.strong().plus(derivativeOrder=0);
                    component = self.resolvedComponent(primaryComponent);
                    row = IMBoundaryRow(family="dirichlet", endpoint=endpoint, component=component, leftOperator=left);
                case "neumann"
                    left = IMOperator.strong().plus(derivativeOrder=1);
                    component = self.resolvedComponent(primaryComponent);
                    row = IMBoundaryRow(family="neumann", endpoint=endpoint, component=component, leftOperator=left);
                case "operator"
                    component = self.resolvedComponent(primaryComponent);
                    innerProductTerms = IMBoundary.placeInnerProductTerms(self.innerProductTerms, location);
                    row = IMBoundaryRow(family="operator", endpoint=endpoint, component=component, ...
                        leftOperator=self.leftOperator, rightOperator=self.rightOperator, ...
                        innerProductTerms=innerProductTerms, hasKnownInnerProduct=self.hasKnownInnerProduct);
                case "linearF"
                    if primaryComponent ~= "F"
                        IMBoundary.unsupportedResolution(self.family, primaryComponent);
                    end
                    row = self.resolveLinearF(endpoint, location);
                case "linearG"
                    if primaryComponent ~= "G"
                        IMBoundary.unsupportedResolution(self.family, primaryComponent);
                    end
                    row = self.resolveLinearG(endpoint, location);
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
            % By default this constrains the EVP primary component. Use
            % `on="F"` or `on="G"` for advanced overrides.
            %
            % ```matlab
            % IMBoundary.dirichlet()        % primary component
            % IMBoundary.dirichlet(on="F")  % explicit F condition
            % ```
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.dirichlet(options)
            % - Parameter options.on: component target, `"primary"`, `"F"`, or `"G"`
            % - Returns boundary: initialized boundary law
            arguments
                options.on {mustBeTextScalar} = "primary"
            end

            boundary = IMBoundary(family="dirichlet", on=options.on);
        end

        function boundary = neumann(options)
            % Create a location-free homogeneous Neumann boundary law.
            %
            % By default this constrains the EVP primary component. Use
            % `on="F"` or `on="G"` for advanced overrides.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.neumann(options)
            % - Parameter options.on: component target, `"primary"`, `"F"`, or `"G"`
            % - Returns boundary: initialized boundary law
            arguments
                options.on {mustBeTextScalar} = "primary"
            end

            boundary = IMBoundary(family="neumann", on=options.on);
        end

        function boundary = operator(options)
            % Create a manual location-free operator boundary law.
            %
            % The law applies `left = lambda right` to the resolved component
            % at whichever endpoint the EVP factory places it. By default it
            % acts on the EVP primary component.
            %
            % ```matlab
            % boundary = IMBoundary.operator(left=IMOperator.strong().plus(derivativeOrder=1));
            % ```
            %
            % `innerProductTerms` are optional boundary trace products
            % associated with the law. If a term is passed with an empty
            % location, it is placed at the resolved endpoint.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.operator(options)
            % - Parameter options.left: left-side boundary functional
            % - Parameter options.right: right-side eigenvalue functional
            % - Parameter options.innerProductTerms: boundary inner-product terms
            % - Parameter options.hasKnownInnerProduct: true when the compatible boundary inner product is known
            % - Parameter options.on: component target, `"primary"`, `"F"`, or `"G"`
            % - Returns boundary: initialized boundary law
            arguments
                options.left IMOperator
                options.right IMOperator = IMOperator.strong()
                options.innerProductTerms struct = IMBoundaryRow.emptyInnerProductTerms()
                options.hasKnownInnerProduct (1,1) logical = true
                options.on {mustBeTextScalar} = "primary"
            end

            boundary = IMBoundary(family="operator", on=options.on, leftOperator=options.left, ...
                rightOperator=options.right, innerProductTerms=options.innerProductTerms, ...
                hasKnownInnerProduct=options.hasKnownInnerProduct);
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
            boundary = IMBoundary(family="rigid");
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
            boundary = IMBoundary(family="noSlip");
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
            boundary = IMBoundary(family="free");
        end

        function boundary = linearF(options)
            % Create a location-free linear `F` boundary law.
            %
            % When resolved, the assembled boundary row represents
            % $$-(aF-bF_z/N^2)=\lambda(cF-dF_z/N^2)/g$$. Supported pairwise
            % cases also declare the compatible boundary inner-product
            % terms. General unresolved coefficient patterns can still be
            % solved, but resolving them warns because their inner-product
            % contribution is unknown.
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

            coefficients = struct("a", options.a, "b", options.b, "c", options.c, "d", options.d);
            boundary = IMBoundary(family="linearF", coefficients=coefficients);
        end

        function boundary = linearG(options)
            % Create a location-free linear `G` boundary law.
            %
            % When resolved, the assembled boundary row represents
            % $$-(eG-aG_z)=\lambda(bG-cG_z)/g$$. Supported pairwise cases
            % also declare the compatible boundary inner-product terms.
            % General unresolved coefficient patterns can still be solved,
            % but resolving them warns because their inner-product
            % contribution is unknown.
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

            coefficients = struct("a", options.a, "b", options.b, "c", options.c, "e", options.e);
            boundary = IMBoundary(family="linearG", coefficients=coefficients);
        end
    end

    methods (Access = private)
        function self = IMBoundary(options)
            arguments
                options.family {mustBeTextScalar} = "operator"
                options.on {mustBeTextScalar} = "primary"
                options.leftOperator IMOperator = IMOperator.strong()
                options.rightOperator IMOperator = IMOperator.strong()
                options.innerProductTerms struct = IMBoundaryRow.emptyInnerProductTerms()
                options.hasKnownInnerProduct (1,1) logical = true
                options.coefficients struct = struct()
            end

            self.family = string(options.family);
            self.on = IMBoundary.validateTarget(options.on);
            self.leftOperator = options.leftOperator;
            self.rightOperator = options.rightOperator;
            self.innerProductTerms = options.innerProductTerms(:);
            self.hasKnownInnerProduct = options.hasKnownInnerProduct;
            self.coefficients = options.coefficients;
        end

        function row = resolveLinearF(self, endpoint, location)
            c = self.coefficients;
            left = IMOperator.strong() ...
                .plus(coefficient=-c.a, derivativeOrder=0) ...
                .plus(coefficient=@(z,ctx) c.b./ctx.N2(z), derivativeOrder=1);
            right = IMOperator.strong() ...
                .plus(coefficient=@(~,ctx) c.c/ctx.g, derivativeOrder=0) ...
                .plus(coefficient=@(z,ctx) -c.d./(ctx.g*ctx.N2(z)), derivativeOrder=1);
            [innerProductTerms, hasKnownInnerProduct] = IMBoundary.linearFInnerProductTerms(location, c.a, c.b, c.c, c.d);
            row = IMBoundaryRow(family="linearF", endpoint=endpoint, component="F", ...
                leftOperator=left, rightOperator=right, innerProductTerms=innerProductTerms, ...
                hasKnownInnerProduct=hasKnownInnerProduct);
        end

        function row = resolveLinearG(self, endpoint, location)
            c = self.coefficients;
            left = IMOperator.strong() ...
                .plus(coefficient=-c.e, derivativeOrder=0) ...
                .plus(coefficient=c.a, derivativeOrder=1);
            right = IMOperator.strong() ...
                .plus(coefficient=@(~,ctx) c.b/ctx.g, derivativeOrder=0) ...
                .plus(coefficient=@(~,ctx) -c.c/ctx.g, derivativeOrder=1);
            [innerProductTerms, hasKnownInnerProduct] = IMBoundary.linearGInnerProductTerms(location, c.a, c.b, c.c, c.e);
            row = IMBoundaryRow(family="linearG", endpoint=endpoint, component="G", ...
                leftOperator=left, rightOperator=right, innerProductTerms=innerProductTerms, ...
                hasKnownInnerProduct=hasKnownInnerProduct);
        end

        function component = resolvedComponent(self, primaryComponent)
            if self.on == "primary"
                component = primaryComponent;
            else
                component = self.on;
            end
        end
    end

    methods (Static, Access = private)
        function [terms, hasKnownInnerProduct] = linearFInnerProductTerms(location, a, b, c, d)
            nonzero = abs([a b c d]) > 0;
            terms = IMBoundaryRow.emptyInnerProductTerms();
            hasKnownInnerProduct = true;
            if nnz(nonzero) > 2
                IMBoundary.warnUnknownInnerProduct("linearF");
                hasKnownInnerProduct = false;
                return;
            end

            g = @(ctx) ctx.g;
            if ~nonzero(1) && ~nonzero(2) && nonzero(3) && nonzero(4)
                terms = IMBoundaryRow.innerProductTerm("G", location, -d/c, IMBoundaryRow.trace("G"), IMBoundaryRow.trace("G"));
            elseif nonzero(1) && nonzero(2) && ~nonzero(3) && ~nonzero(4)
                terms = IMBoundaryRow.innerProductTerm("G", location, @(ctx) -b/(g(ctx)*a), ...
                    IMBoundaryRow.trace("G"), IMBoundaryRow.trace("G"));
            elseif ~nonzero(1) && nonzero(2) && nonzero(3) && ~nonzero(4)
                terms = IMBoundaryRow.innerProductTerm("F", location, c/b, IMBoundaryRow.trace("F"), IMBoundaryRow.trace("F"));
            elseif nonzero(1) && ~nonzero(2) && ~nonzero(3) && nonzero(4)
                terms = [
                    IMBoundaryRow.innerProductTerm("F", location, @(ctx) -d/(g(ctx)*g(ctx)*a), ...
                        IMBoundaryRow.trace("G"), IMBoundaryRow.trace("G"))
                    IMBoundaryRow.innerProductTerm("G", location, 1, IMBoundaryRow.trace("G"), IMBoundaryRow.trace("F"))
                    IMBoundaryRow.innerProductTerm("G", location, 1, IMBoundaryRow.trace("F"), IMBoundaryRow.trace("G"))
                ];
            end
        end

        function [terms, hasKnownInnerProduct] = linearGInnerProductTerms(location, a, b, c, e)
            nonzero = abs([a b c e]) > 0;
            terms = IMBoundaryRow.emptyInnerProductTerms();
            hasKnownInnerProduct = true;
            if nnz(nonzero) > 2
                IMBoundary.warnUnknownInnerProduct("linearG");
                hasKnownInnerProduct = false;
                return;
            end

            if ~nonzero(1) && ~nonzero(2) && nonzero(3) && nonzero(4)
                terms = IMBoundaryRow.innerProductTerm("G", location, -c/e, ...
                    IMBoundaryRow.trace("G", derivativeOrder=1), IMBoundaryRow.trace("G", derivativeOrder=1));
            elseif nonzero(1) && nonzero(2) && ~nonzero(3) && ~nonzero(4)
                terms = IMBoundaryRow.innerProductTerm("G", location, b/a, IMBoundaryRow.trace("G"), IMBoundaryRow.trace("G"));
            elseif nonzero(1) && ~nonzero(2) && ~nonzero(3) && nonzero(4)
                terms = IMBoundaryRow.innerProductTerm("F", location, -a/e, IMBoundaryRow.trace("F"), IMBoundaryRow.trace("F"));
            elseif ~nonzero(1) && nonzero(2) && nonzero(3) && ~nonzero(4)
                terms = IMBoundaryRow.innerProductTerm("F", location, -c/b, IMBoundaryRow.trace("F"), IMBoundaryRow.trace("F"));
            end
        end

        function terms = placeInnerProductTerms(terms, location)
            for iTerm = 1:length(terms)
                if string(terms(iTerm).location) == ""
                    terms(iTerm).location = string(location);
                end
            end
        end

        function component = validateComponent(component)
            component = string(component);
            if component ~= "F" && component ~= "G"
                error("IMBoundary:InvalidComponent", ...
                    "component must be ""F"" or ""G"".");
            end
        end

        function target = validateTarget(target)
            target = string(target);
            if target ~= "primary" && target ~= "F" && target ~= "G"
                error("IMBoundary:InvalidTarget", ...
                    "on must be ""primary"", ""F"", or ""G"".");
            end
        end

        function unsupportedResolution(family, primaryComponent)
            error("IMBoundary:UnsupportedResolution", ...
                "Boundary family ""%s"" cannot resolve for primary component ""%s"".", ...
                string(family), string(primaryComponent));
        end

        function warnUnknownInnerProduct(family)
            warning("IMBoundary:UnknownInnerProduct", ...
                "%s can assemble this boundary row, but its compatible boundary inner-product contribution is unknown.", ...
                family);
        end
    end
end

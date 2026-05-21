classdef IMBoundary
    % Describe an internal-mode boundary condition.
    %
    % `IMBoundary` stores a mathematical boundary condition. A boundary can
    % be location-free, such as `IMBoundary.free()`, or placed at a physical
    % endpoint, such as `IMBoundary.free().at("surface",formulation="G")`.
    % EVP factories normally accept location-free surface and bottom laws
    % and place them through `IMBoundary.conditions(...)`.
    %
    % Boundary operators use the physical coordinate derivative
    % $$\partial_z$$ at both endpoints. A boundary functional has the form
    % $$\sum_p a_p(z_\ell)\partial_z^p u(z_\ell)$$ at endpoint
    % $$z_\ell$$; no outward-normal derivative sign is inserted at the
    % bottom.
    %
    % Boundary inner-product terms are not boundary conditions. They are the
    % boundary trace products added to modal inner products so the EVP
    % remains orthogonal under the chosen boundary law. Endpoint orientation
    % enters through Green's identity:
    % $$s_\ell=+1$$ at the surface and $$s_\ell=-1$$ at the bottom, so a
    % generated term with coefficient $$c$$ contributes
    % $$s_\ell c\,u_i(z_\ell)v_j(z_\ell)$$.
    %
    % ```matlab
    % boundaryConditions = IMBoundary.conditions(formulation="G", ...
    %     surface=IMBoundary.free(), bottom=IMBoundary.rigid());
    % ```
    %
    % - Topic: Create boundary conditions
    % - Topic: Place boundary conditions
    % - Topic: Inspect boundary conditions
    % - Topic: Developer topics
    % - Declaration: classdef IMBoundary

    properties (SetAccess = private)
        % Boundary family identifier.
        %
        % Examples include `"rigid"`, `"noSlip"`, `"free"`, `"linearF"`,
        % `"linearG"`, and `"custom"`.
        %
        % - Topic: Inspect boundary conditions
        family = "custom"

        % Physical endpoint location.
        %
        % Empty for location-free boundary conditions. Placed boundary
        % conditions use `"surface"` or `"bottom"`.
        %
        % - Topic: Inspect boundary conditions
        location = ""

        % Variable constrained by this boundary condition.
        %
        % Location-free boundary conditions may use `"formulation"` as a
        % target placeholder. Placed boundary conditions use `"F"` or `"G"`.
        %
        % - Topic: Inspect boundary conditions
        variable = "formulation"

        % Left-side boundary functional.
        %
        % The assembled condition is `leftOperator = lambda rightOperator`
        % at the boundary location.
        %
        % - Topic: Inspect boundary conditions
        leftOperator = IMOperator()

        % Right-side eigenvalue boundary functional.
        %
        % - Topic: Inspect boundary conditions
        rightOperator = IMOperator()

        % Boundary terms implied by this condition for modal inner products.
        %
        % Each term has an inner-product variable, physical location,
        % coefficient, left trace, and right trace.
        %
        % - Topic: Inspect boundary conditions
        innerProductTerms = IMBoundary.emptyInnerProductTerms()

        % True when the compatible boundary inner-product terms are known.
        %
        % A condition can still assemble when this is false, but Gram
        % matrices, normalization, and index metadata may be incomplete for
        % boundary modes.
        %
        % - Topic: Inspect boundary conditions
        hasKnownInnerProductTerms = true

        % Sign of the boundary-mode eigenvalue contribution.
        %
        % This sign is used when selecting a declared endpoint boundary
        % branch. Negative signs also contribute to the expected negative
        % eigenvalue index.
        %
        % - Topic: Inspect boundary conditions
        indexSign = 0

        % Rank of the boundary-mode eigenvalue contribution.
        %
        % - Topic: Inspect boundary conditions
        indexRank = 0

        % Physical mode number for an endpoint boundary branch.
        %
        % Boundary-mode numbers are physical labels, not eigenvalue signs.
        % The surface branch is always `-1`, the bottom branch is always
        % `-2`, and `NaN` means this condition does not declare a boundary
        % branch.
        %
        % - Topic: Inspect boundary conditions
        boundaryModeNumber = NaN
    end

    properties (Access = private)
        coefficients = struct()
    end

    methods
        function boundary = at(self, location, options)
            % Place a boundary condition at a physical endpoint.
            %
            % Location-free families such as `rigid`, `free`, and `noSlip`
            % resolve their variable-dependent formulas when placed. The
            % assembled operator still uses coordinate $$\partial_z$$ at
            % both endpoints.
            %
            % ```matlab
            % condition = IMBoundary.free().at("surface", formulation="G");
            % ```
            %
            % - Topic: Place boundary conditions
            % - Declaration: boundary = at(boundary,location,options)
            % - Parameter location: `"surface"` or `"bottom"`
            % - Parameter options.formulation: EVP formulation, `"F"` or `"G"`
            % - Returns boundary: placed boundary condition
            arguments
                self IMBoundary
                location {mustBeTextScalar}
                options.formulation {mustBeTextScalar} = "formulation"
            end

            location = IMBoundary.validateLocation(location);
            if self.location ~= ""
                if self.location ~= location
                    error("IMBoundary:AlreadyPlaced", ...
                        "Boundary condition is already placed at ""%s"".", self.location);
                end
                boundary = self;
                return;
            end

            requestedVariable = IMBoundary.validateTarget(options.formulation);
            if requestedVariable == "formulation"
                requestedVariable = self.variable;
            end
            if requestedVariable == "formulation"
                error("IMBoundary:MissingVariable", ...
                    "A location-free boundary condition must be placed with formulation=""F"" or formulation=""G"".");
            end
            formulation = IMBoundary.validateVariable(requestedVariable);
            variable = self.resolvedVariable(formulation);
            endpointSign = IMBoundary.endpointSign(location);

            switch self.family
                case "rigid"
                    switch variable
                        case "G"
                            left = IMOperator().plus(derivativeOrder=0);
                            boundary = IMBoundary(family="rigid", location=location, variable="G", leftOperator=left);
                        case "F"
                            left = IMOperator().plus(derivativeOrder=1);
                            boundary = IMBoundary(family="rigid", location=location, variable="F", leftOperator=left);
                    end
                case "noSlip"
                    switch variable
                        case "G"
                            left = IMOperator().plus(derivativeOrder=1);
                            boundary = IMBoundary(family="noSlip", location=location, variable="G", leftOperator=left);
                        case "F"
                            left = IMOperator().plus(derivativeOrder=0);
                            boundary = IMBoundary(family="noSlip", location=location, variable="F", leftOperator=left);
                    end
                case "free"
                    innerProductTerm = IMBoundary.innerProductTerm("G", location, endpointSign, ...
                        IMBoundary.trace("G"), IMBoundary.trace("G"));
                    switch variable
                        case "G"
                            left = IMOperator().plus(derivativeOrder=1);
                            right = IMOperator().plus(derivativeOrder=0);
                            boundary = IMBoundary(family="free", location=location, variable="G", ...
                                leftOperator=left, rightOperator=right, innerProductTerms=innerProductTerm, ...
                                indexSign=endpointSign, indexRank=1, ...
                                boundaryModeNumber=IMBoundary.boundaryModeNumberForLocation(location));
                        case "F"
                            left = IMOperator() ...
                                .plus(derivativeOrder=0) ...
                                .plus(coefficient=@(z,ctx) ctx.g./ctx.N2(z), derivativeOrder=1);
                            boundary = IMBoundary(family="free", location=location, variable="F", ...
                                leftOperator=left, innerProductTerms=innerProductTerm);
                    end
                case "dirichlet"
                    left = IMOperator().plus(derivativeOrder=0);
                    boundary = IMBoundary(family="dirichlet", location=location, variable=variable, leftOperator=left);
                case "neumann"
                    left = IMOperator().plus(derivativeOrder=1);
                    boundary = IMBoundary(family="neumann", location=location, variable=variable, leftOperator=left);
                case "custom"
                    innerProductTerms = IMBoundary.placeInnerProductTerms(self.innerProductTerms, location);
                    boundary = IMBoundary(family="custom", location=location, variable=variable, ...
                        leftOperator=self.leftOperator, rightOperator=self.rightOperator, ...
                        innerProductTerms=innerProductTerms, hasKnownInnerProductTerms=self.hasKnownInnerProductTerms, ...
                        indexSign=self.indexSign, indexRank=self.indexRank, boundaryModeNumber=self.boundaryModeNumber);
                case "linearF"
                    if formulation ~= "F"
                        IMBoundary.unsupportedPlacement(self.family, formulation);
                    end
                    boundary = self.placeLinearF(location);
                case "linearG"
                    if formulation ~= "G"
                        IMBoundary.unsupportedPlacement(self.family, formulation);
                    end
                    boundary = self.placeLinearG(location);
                otherwise
                    error("IMBoundary:UnsupportedPlacement", ...
                        "Boundary family ""%s"" cannot be placed by the boundary API.", self.family);
            end
        end

        function tf = hasInnerProductTerms(self)
            % Return true when this condition contributes inner-product terms.
            %
            % - Topic: Inspect boundary conditions
            % - Declaration: tf = hasInnerProductTerms(boundary)
            % - Returns tf: true when boundary inner-product terms are present
            tf = ~isempty(self.innerProductTerms);
        end

        function count = expectedNegativeCount(self)
            % Return the negative-index contribution from this condition.
            %
            % Only conditions that declare an endpoint boundary mode
            % contribute expected index counts.
            %
            % - Topic: Inspect boundary conditions
            % - Declaration: count = expectedNegativeCount(boundary)
            % - Returns count: expected negative count contribution
            if ~self.hasKnownInnerProductTerms || isnan(self.boundaryModeNumber)
                count = 0;
            elseif self.indexSign < 0
                count = self.indexRank;
            else
                count = 0;
            end
        end

        function count = expectedZeroCount(self)
            % Return the zero-index contribution from this condition.
            %
            % Only conditions that declare an endpoint boundary mode
            % contribute expected index counts.
            %
            % - Topic: Inspect boundary conditions
            % - Declaration: count = expectedZeroCount(boundary)
            % - Returns count: expected zero count contribution
            if ~self.hasKnownInnerProductTerms || isnan(self.boundaryModeNumber)
                count = 0;
            elseif self.indexSign == 0
                count = self.indexRank;
            else
                count = 0;
            end
        end

        function descriptors = boundaryModeDescriptors(self)
            % Return declared endpoint boundary-mode metadata.
            %
            % Boundary-mode numbers are endpoint labels used for retained
            % modes. They are selected using `indexSign`, but the sign of the
            % eigenvalue and the physical mode number are intentionally
            % separate concepts.
            %
            % - Topic: Inspect boundary conditions
            % - Declaration: descriptors = boundaryModeDescriptors(boundary)
            % - Returns descriptors: structure array with `modeNumber` and `indexSign`
            descriptors = struct("modeNumber", {}, "indexSign", {});
            if ~self.hasKnownInnerProductTerms || self.indexRank == 0 || isnan(self.boundaryModeNumber)
                return;
            end
            descriptors = struct("modeNumber", self.boundaryModeNumber, "indexSign", self.indexSign);
        end
    end

    methods (Static)
        function boundaryConditions = conditions(options)
            % Place bottom and surface boundary conditions for one variable.
            %
            % The returned array is ordered bottom first and surface second,
            % matching the usual vertical domain order.
            %
            % - Topic: Place boundary conditions
            % - Declaration: boundaryConditions = IMBoundary.conditions(options)
            % - Parameter options.formulation: EVP formulation, `"F"` or `"G"`
            % - Parameter options.surface: location-free surface boundary law
            % - Parameter options.bottom: location-free bottom boundary law
            % - Returns boundaryConditions: placed boundary-condition array
            arguments
                options.formulation {mustBeTextScalar}
                options.surface (1,1) IMBoundary
                options.bottom (1,1) IMBoundary
            end

            variable = IMBoundary.validateVariable(options.formulation);
            boundaryConditions = [
                options.bottom.at("bottom", formulation=variable)
                options.surface.at("surface", formulation=variable)
            ];
        end

        function boundary = dirichlet()
            % Create a location-free homogeneous Dirichlet boundary law.
            %
            % This constrains the EVP formulation when placed.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundary.dirichlet()
            % - Returns boundary: initialized boundary condition
            boundary = IMBoundary(family="dirichlet");
        end

        function boundary = neumann()
            % Create a location-free homogeneous Neumann boundary law.
            %
            % This constrains the EVP formulation when placed.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundary.neumann()
            % - Returns boundary: initialized boundary condition
            boundary = IMBoundary(family="neumann");
        end

        function boundary = custom(options)
            % Create a custom location-free operator boundary law.
            %
            % The law applies `left = lambda right` at whichever endpoint
            % the EVP factory places it. Both operators are written with the
            % physical coordinate derivative $$\partial_z$$ at the endpoint.
            %
            % `innerProductTerms` are optional boundary trace products
            % associated with the law. If a term is passed with an empty
            % location, it is placed at the resolved endpoint and receives
            % the endpoint orientation sign from Green's identity. Explicitly
            % located terms are treated as final and are not reoriented.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundary.custom(options)
            % - Parameter options.left: left-side boundary functional
            % - Parameter options.right: right-side eigenvalue functional
            % - Parameter options.innerProductTerms: boundary inner-product terms
            % - Parameter options.hasKnownInnerProductTerms: true when the compatible boundary inner-product terms are known
            % - Parameter options.variable: variable target, `"formulation"`, `"F"`, or `"G"`
            % - Parameter options.indexSign: expected boundary-mode eigenvalue sign, `-1`, `0`, or `1`
            % - Parameter options.indexRank: number of boundary-mode directions; currently must be `1` when `boundaryModeNumber` is supplied
            % - Parameter options.boundaryModeNumber: explicit endpoint mode number, `-1` for surface or `-2` for bottom
            % - Returns boundary: initialized boundary condition
            arguments
                options.left IMOperator
                options.right IMOperator = IMOperator()
                options.innerProductTerms struct = IMBoundary.emptyInnerProductTerms()
                options.hasKnownInnerProductTerms (1,1) logical = true
                options.variable {mustBeTextScalar} = "formulation"
                options.indexSign (1,1) double {mustBeMember(options.indexSign, [-1 0 1])} = 0
                options.indexRank (1,1) double {mustBeInteger, mustBeNonnegative} = 0
                options.boundaryModeNumber (1,1) double = NaN
            end

            boundary = IMBoundary(family="custom", variable=options.variable, ...
                leftOperator=options.left, rightOperator=options.right, ...
                innerProductTerms=options.innerProductTerms, ...
                hasKnownInnerProductTerms=options.hasKnownInnerProductTerms, ...
                indexSign=options.indexSign, indexRank=options.indexRank, ...
                boundaryModeNumber=options.boundaryModeNumber);
        end

        function boundary = rigid()
            % Create a location-free rigid boundary law.
            %
            % In a `G` EVP, `rigid` resolves to $$G=0$$. In an `F` EVP, it
            % resolves to $$F_z=0$$.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundary.rigid()
            % - Returns boundary: initialized rigid boundary condition
            boundary = IMBoundary(family="rigid");
        end

        function boundary = noSlip()
            % Create a location-free no-slip boundary law.
            %
            % In a `G` EVP, `noSlip` resolves to $$G_z=0$$. In an `F` EVP,
            % it resolves to $$F=0$$.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundary.noSlip()
            % - Returns boundary: initialized no-slip boundary condition
            boundary = IMBoundary(family="noSlip");
        end

        function boundary = free()
            % Create a location-free free boundary law.
            %
            % In a `G` EVP, `free` resolves to $$G_z=\lambda G$$. In an `F`
            % EVP, it resolves to $$F+gF_z/N^2=0$$. The derivative $$G_z$$
            % or $$F_z$$ always means the coordinate derivative
            % $$\partial_z$$. A placed `G` free boundary declares an endpoint
            % boundary branch with mode number `-1` at the surface or `-2`
            % at the bottom.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundary.free()
            % - Returns boundary: initialized free boundary condition
            boundary = IMBoundary(family="free");
        end

        function boundary = linearF(options)
            % Create a location-free linear `F` boundary law.
            %
            % When placed, the assembled boundary condition represents
            % $$-(aF-bF_z/N^2)=\lambda(cF-dF_z/N^2)/g$$. Supported pairwise
            % cases also declare compatible boundary inner-product terms.
            % General unresolved coefficient patterns can still be solved,
            % but placing them warns because their inner-product contribution
            % is unknown. The derivative $$F_z$$ always means
            % $$\partial_z F$$ at either endpoint.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundary.linearF(options)
            % - Parameter options.a: coefficient multiplying `F`
            % - Parameter options.b: coefficient multiplying `F_z/N2`
            % - Parameter options.c: eigenvalue coefficient multiplying `F`
            % - Parameter options.d: eigenvalue coefficient multiplying `F_z/N2`
            % - Returns boundary: initialized linear `F` boundary condition
            arguments
                options.a (1,1) double = 0
                options.b (1,1) double = 0
                options.c (1,1) double = 0
                options.d (1,1) double = 0
            end

            coefficients = struct("a", options.a, "b", options.b, "c", options.c, "d", options.d);
            boundary = IMBoundary(family="linearF", variable="F", coefficients=coefficients);
        end

        function boundary = linearG(options)
            % Create a location-free linear `G` boundary law.
            %
            % When placed, the assembled boundary condition represents
            % $$-(eG-aG_z)=\lambda(bG-cG_z)/g$$. Supported pairwise cases
            % also declare compatible boundary inner-product terms. General
            % unresolved coefficient patterns can still be solved, but
            % placing them warns because their inner-product contribution is
            % unknown. The derivative $$G_z$$ always means
            % $$\partial_z G$$ at either endpoint.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundary.linearG(options)
            % - Parameter options.a: coefficient multiplying `G_z`
            % - Parameter options.b: eigenvalue coefficient multiplying `G`
            % - Parameter options.c: eigenvalue coefficient multiplying `G_z`
            % - Parameter options.e: coefficient multiplying `G`
            % - Returns boundary: initialized linear `G` boundary condition
            arguments
                options.a (1,1) double = 0
                options.b (1,1) double = 0
                options.c (1,1) double = 0
                options.e (1,1) double = 0
            end

            coefficients = struct("a", options.a, "b", options.b, "c", options.c, "e", options.e);
            boundary = IMBoundary(family="linearG", variable="G", coefficients=coefficients);
        end

        function boundary = active(options)
            % Create an active metadata-only boundary condition.
            %
            % Active conditions are already placed because their trace terms
            % refer to an endpoint of a partial-depth interval. They declare
            % endpoint boundary-mode numbers using `-1` at the surface and
            % `-2` at the bottom.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundary = IMBoundary.active(options)
            % - Parameter options.location: `"surface"` or `"bottom"`
            % - Parameter options.variable: active variable
            % - Parameter options.indexSign: active-boundary sign
            % - Parameter options.indexRank: number of active directions
            % - Parameter options.innerProductTerms: boundary inner-product terms
            % - Returns boundary: initialized active boundary condition
            arguments
                options.location {mustBeTextScalar}
                options.variable {mustBeTextScalar} = "G"
                options.indexSign (1,1) double {mustBeMember(options.indexSign, [-1 1])}
                options.indexRank (1,1) double {mustBeInteger, mustBePositive} = 1
                options.innerProductTerms struct = IMBoundary.emptyInnerProductTerms()
            end

            boundary = IMBoundary(family="active", location=options.location, variable=options.variable, ...
                innerProductTerms=options.innerProductTerms, indexSign=options.indexSign, indexRank=options.indexRank, ...
                boundaryModeNumber=IMBoundary.boundaryModeNumberForLocation(options.location));
        end

        function boundaryConditions = partialDepthPE(options)
            % Create partial-depth potential-energy active boundary conditions.
            %
            % Positive boundary signs add no negative index directions.
            % Negative boundary signs add one negative direction at each
            % window endpoint.
            %
            % - Topic: Create boundary conditions
            % - Declaration: boundaryConditions = IMBoundary.partialDepthPE(options)
            % - Parameter options.boundarySign: `"positive"` or `"negative"`
            % - Returns boundaryConditions: bottom and surface active conditions
            arguments
                options.boundarySign {mustBeTextScalar} = "positive"
            end

            switch string(options.boundarySign)
                case "positive"
                    indexSign = 1;
                case "negative"
                    indexSign = -1;
                otherwise
                    error("IMBoundary:InvalidBoundarySign", ...
                        "boundarySign must be ""positive"" or ""negative"".");
            end
            boundaryConditions = [
                IMBoundary(family="partialDepthPE", location="bottom", variable="G", ...
                    indexSign=indexSign, indexRank=1, boundaryModeNumber=-2)
                IMBoundary(family="partialDepthPE", location="surface", variable="G", ...
                    indexSign=indexSign, indexRank=1, boundaryModeNumber=-1)
            ];
        end

        function trace = trace(variable, options)
            % Create an endpoint trace descriptor.
            %
            % A trace describes which variable value or first derivative is
            % evaluated at a boundary endpoint.
            %
            % - Topic: Create boundary conditions
            % - Declaration: trace = IMBoundary.trace(variable,options)
            % - Parameter variable: variable name
            % - Parameter options.derivativeOrder: physical derivative order
            % - Returns trace: endpoint trace descriptor
            arguments
                variable {mustBeTextScalar}
                options.derivativeOrder (1,1) double {mustBeInteger, mustBeNonnegative} = 0
            end

            trace = struct("variable", IMBoundary.validateVariable(variable), "derivativeOrder", options.derivativeOrder);
        end

        function term = innerProductTerm(innerProductVariable, location, coefficient, leftTrace, rightTrace)
            % Create a boundary inner-product trace-pair term.
            %
            % The term contributes `coefficient*leftTrace_i*rightTrace_j`
            % at a boundary endpoint to the named variable's inner product.
            %
            % - Topic: Create boundary conditions
            % - Declaration: term = IMBoundary.innerProductTerm(innerProductVariable,location,coefficient,leftTrace,rightTrace)
            % - Parameter innerProductVariable: variable whose inner product receives the term
            % - Parameter location: boundary location
            % - Parameter coefficient: scalar or context function handle
            % - Parameter leftTrace: trace evaluated for the left mode
            % - Parameter rightTrace: trace evaluated for the right mode
            % - Returns term: boundary inner-product term
            location = string(location);
            if location ~= ""
                location = IMBoundary.validateLocation(location);
            end
            term = struct("innerProductVariable", IMBoundary.validateVariable(innerProductVariable), ...
                "location", location, "coefficient", coefficient, "leftTrace", leftTrace, "rightTrace", rightTrace);
        end

        function terms = emptyInnerProductTerms()
            % Return an empty boundary inner-product-term structure.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: terms = IMBoundary.emptyInnerProductTerms()
            % - Returns terms: empty boundary inner-product-term structure
            terms = struct("innerProductVariable", {}, "location", {}, "coefficient", {}, ...
                "leftTrace", {}, "rightTrace", {});
        end
    end

    methods (Access = private)
        function self = IMBoundary(options)
            arguments
                options.family {mustBeTextScalar} = "custom"
                options.location {mustBeTextScalar} = ""
                options.variable {mustBeTextScalar} = "formulation"
                options.leftOperator IMOperator = IMOperator()
                options.rightOperator IMOperator = IMOperator()
                options.innerProductTerms struct = IMBoundary.emptyInnerProductTerms()
                options.hasKnownInnerProductTerms (1,1) logical = true
                options.indexSign (1,1) double = 0
                options.indexRank (1,1) double {mustBeInteger, mustBeNonnegative} = 0
                options.boundaryModeNumber (1,1) double = NaN
                options.coefficients struct = struct()
            end

            self.family = string(options.family);
            self.location = string(options.location);
            if self.location ~= ""
                self.location = IMBoundary.validateLocation(self.location);
            end
            if self.location == ""
                self.variable = IMBoundary.validateTarget(options.variable);
            else
                self.variable = IMBoundary.validateVariable(options.variable);
            end
            self.leftOperator = options.leftOperator;
            self.rightOperator = options.rightOperator;
            self.innerProductTerms = options.innerProductTerms(:);
            self.hasKnownInnerProductTerms = options.hasKnownInnerProductTerms;
            self.indexSign = sign(options.indexSign);
            self.indexRank = options.indexRank;
            self.boundaryModeNumber = IMBoundary.validateBoundaryModeNumber(options.boundaryModeNumber);
            if ~isnan(self.boundaryModeNumber) && self.indexRank ~= 1
                error("IMBoundary:UnsupportedBoundaryModeRank", ...
                    "Boundary-mode labels support exactly one branch at each endpoint.");
            end
            self.coefficients = options.coefficients;
        end

        function boundary = placeLinearF(self, location)
            c = self.coefficients;
            left = IMOperator() ...
                .plus(coefficient=-c.a, derivativeOrder=0) ...
                .plus(coefficient=@(z,ctx) c.b./ctx.N2(z), derivativeOrder=1);
            right = IMOperator() ...
                .plus(coefficient=@(~,ctx) c.c/ctx.g, derivativeOrder=0) ...
                .plus(coefficient=@(z,ctx) -c.d./(ctx.g*ctx.N2(z)), derivativeOrder=1);
            [innerProductTerms, hasKnownInnerProductTerms] = IMBoundary.linearFInnerProductTerms(location, c.a, c.b, c.c, c.d);
            boundary = IMBoundary(family="linearF", location=location, variable="F", ...
                leftOperator=left, rightOperator=right, innerProductTerms=innerProductTerms, ...
                hasKnownInnerProductTerms=hasKnownInnerProductTerms);
        end

        function boundary = placeLinearG(self, location)
            c = self.coefficients;
            left = IMOperator() ...
                .plus(coefficient=-c.e, derivativeOrder=0) ...
                .plus(coefficient=c.a, derivativeOrder=1);
            right = IMOperator() ...
                .plus(coefficient=@(~,ctx) c.b/ctx.g, derivativeOrder=0) ...
                .plus(coefficient=@(~,ctx) -c.c/ctx.g, derivativeOrder=1);
            [innerProductTerms, hasKnownInnerProductTerms] = IMBoundary.linearGInnerProductTerms(location, c.a, c.b, c.c, c.e);
            boundary = IMBoundary(family="linearG", location=location, variable="G", ...
                leftOperator=left, rightOperator=right, innerProductTerms=innerProductTerms, ...
                hasKnownInnerProductTerms=hasKnownInnerProductTerms);
        end

        function variable = resolvedVariable(self, formulation)
            if self.variable == "formulation"
                variable = formulation;
            else
                variable = self.variable;
            end
        end
    end

    methods (Static, Access = private)
        function [terms, hasKnownInnerProductTerms] = linearFInnerProductTerms(location, a, b, c, d)
            nonzero = abs([a b c d]) > 0;
            terms = IMBoundary.emptyInnerProductTerms();
            hasKnownInnerProductTerms = true;
            if nnz(nonzero) > 2
                IMBoundary.warnUnknownInnerProduct("linearF");
                hasKnownInnerProductTerms = false;
                return;
            end

            endpointSign = IMBoundary.endpointSign(location);
            g = @(ctx) ctx.g;
            if ~nonzero(1) && ~nonzero(2) && nonzero(3) && nonzero(4)
                terms = IMBoundary.innerProductTerm("G", location, endpointSign*(-d/c), ...
                    IMBoundary.trace("G"), IMBoundary.trace("G"));
            elseif nonzero(1) && nonzero(2) && ~nonzero(3) && ~nonzero(4)
                terms = IMBoundary.innerProductTerm("G", location, @(ctx) endpointSign*(-b/(g(ctx)*a)), ...
                    IMBoundary.trace("G"), IMBoundary.trace("G"));
            elseif ~nonzero(1) && nonzero(2) && nonzero(3) && ~nonzero(4)
                terms = IMBoundary.innerProductTerm("F", location, endpointSign*(c/b), ...
                    IMBoundary.trace("F"), IMBoundary.trace("F"));
            elseif nonzero(1) && ~nonzero(2) && ~nonzero(3) && nonzero(4)
                terms = [
                    IMBoundary.innerProductTerm("F", location, @(ctx) endpointSign*(-d/(g(ctx)*g(ctx)*a)), ...
                        IMBoundary.trace("G"), IMBoundary.trace("G"))
                    IMBoundary.innerProductTerm("G", location, endpointSign, IMBoundary.trace("G"), IMBoundary.trace("F"))
                    IMBoundary.innerProductTerm("G", location, endpointSign, IMBoundary.trace("F"), IMBoundary.trace("G"))
                ];
            end
        end

        function [terms, hasKnownInnerProductTerms] = linearGInnerProductTerms(location, a, b, c, e)
            nonzero = abs([a b c e]) > 0;
            terms = IMBoundary.emptyInnerProductTerms();
            hasKnownInnerProductTerms = true;
            if nnz(nonzero) > 2
                IMBoundary.warnUnknownInnerProduct("linearG");
                hasKnownInnerProductTerms = false;
                return;
            end

            endpointSign = IMBoundary.endpointSign(location);
            if ~nonzero(1) && ~nonzero(2) && nonzero(3) && nonzero(4)
                terms = IMBoundary.innerProductTerm("G", location, endpointSign*(-c/e), ...
                    IMBoundary.trace("G", derivativeOrder=1), IMBoundary.trace("G", derivativeOrder=1));
            elseif nonzero(1) && nonzero(2) && ~nonzero(3) && ~nonzero(4)
                terms = IMBoundary.innerProductTerm("G", location, endpointSign*(b/a), ...
                    IMBoundary.trace("G"), IMBoundary.trace("G"));
            elseif nonzero(1) && ~nonzero(2) && ~nonzero(3) && nonzero(4)
                terms = IMBoundary.innerProductTerm("F", location, endpointSign*(-a/e), ...
                    IMBoundary.trace("F"), IMBoundary.trace("F"));
            elseif ~nonzero(1) && nonzero(2) && nonzero(3) && ~nonzero(4)
                terms = IMBoundary.innerProductTerm("F", location, endpointSign*(-c/b), ...
                    IMBoundary.trace("F"), IMBoundary.trace("F"));
            end
        end

        function terms = placeInnerProductTerms(terms, location)
            for iTerm = 1:length(terms)
                if string(terms(iTerm).location) == ""
                    terms(iTerm).location = string(location);
                    if string(location) == "bottom"
                        coefficient = terms(iTerm).coefficient;
                        if isa(coefficient, "function_handle")
                            terms(iTerm).coefficient = @(varargin) -coefficient(varargin{:});
                        else
                            terms(iTerm).coefficient = -coefficient;
                        end
                    end
                end
            end
        end

        function signValue = endpointSign(location)
            location = IMBoundary.validateLocation(location);
            if location == "bottom"
                signValue = -1;
            else
                signValue = 1;
            end
        end

        function modeNumber = boundaryModeNumberForLocation(location)
            location = IMBoundary.validateLocation(location);
            if location == "surface"
                modeNumber = -1;
            else
                modeNumber = -2;
            end
        end

        function modeNumber = validateBoundaryModeNumber(modeNumber)
            if isnan(modeNumber)
                return;
            end
            if modeNumber ~= -1 && modeNumber ~= -2
                error("IMBoundary:InvalidBoundaryModeNumber", ...
                    "boundaryModeNumber must be -1 for the surface branch, -2 for the bottom branch, or NaN.");
            end
        end

        function variable = validateVariable(variable)
            variable = string(variable);
            if variable ~= "F" && variable ~= "G"
                error("IMBoundary:InvalidVariable", ...
                    "variable must be ""F"" or ""G"".");
            end
        end

        function target = validateTarget(target)
            target = string(target);
            if target ~= "formulation" && target ~= "F" && target ~= "G"
                error("IMBoundary:InvalidTarget", ...
                    "variable must be ""formulation"", ""F"", or ""G"".");
            end
        end

        function location = validateLocation(location)
            location = string(location);
            if location ~= "surface" && location ~= "bottom"
                error("IMBoundary:InvalidLocation", ...
                    "location must be ""surface"" or ""bottom"".");
            end
        end

        function unsupportedPlacement(family, variable)
            error("IMBoundary:UnsupportedPlacement", ...
                "Boundary family ""%s"" cannot place for variable ""%s"".", ...
                string(family), string(variable));
        end

        function warnUnknownInnerProduct(family)
            warning("IMBoundary:UnknownInnerProduct", ...
                "%s can assemble this boundary condition, but its compatible boundary inner-product contribution is unknown.", ...
                family);
        end
    end
end

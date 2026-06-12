classdef IMBoundary
    % Describe an internal-mode boundary law.
    %
    % `IMBoundary` stores the physical boundary law supplied to an
    % `IMEigenvalueProblem`. Standard EVPs keep one law at the surface and
    % one law at the bottom. When an EVP is assembled, each law is resolved
    % with the EVP formulation and endpoint to produce the operator row,
    % endpoint weights, and mode-index metadata needed by the solver and
    % basis-set layers.
    %
    % Boundary operators use the physical coordinate derivative
    % $$\partial_z$$. A boundary functional has the form
    % $$\sum_p a_p(z_\ell)\partial_z^p u(z_\ell)$$ at endpoint
    % $$z_\ell$$; no outward-normal derivative sign is inserted at the
    % bottom.
    %
    % Boundary weights are not boundary laws. They are endpoint products
    % added to modal inner products so the EVP remains orthogonal under the
    % chosen boundary law. Endpoint orientation enters through Green's
    % identity:
    % $$s_\ell=+1$$ at the surface and $$s_\ell=-1$$ at the bottom, so a
    % generated term with coefficient $$c$$ contributes
    % $$s_\ell c\,u_i(z_\ell)v_j(z_\ell)$$.
    %
    % ```matlab
    % evp = IMEigenvalueProblem.hydrostaticGModes( ...
    %     surfaceBoundary=IMBoundary.free(), bottomBoundary=IMBoundary.rigid());
    % ```
    %
    % - Topic: Create boundary laws
    % - Topic: Inspect boundary laws
    % - Topic: Developer topics
    % - Declaration: classdef IMBoundary

    properties (SetAccess = private)
        % Boundary family identifier.
        %
        % Examples include `"rigid"`, `"noSlip"`, `"free"`, `"linearF"`,
        % `"linearG"`, and `"custom"`.
        %
        % - Topic: Inspect boundary laws
        family = "custom"

        % Physical endpoint location.
        %
        % Empty for location-free boundary laws. Internally resolved laws
        % use `"surface"` or `"bottom"`.
        %
        % - Topic: Inspect boundary laws
        location = ""

        % Variable constrained by this boundary law.
        %
        % Location-free laws may use `"formulation"` as a target placeholder.
        % Resolved laws use `"F"` or `"G"`.
        %
        % - Topic: Inspect boundary laws
        variable = "formulation"

        % Left-side boundary functional.
        %
        % The assembled condition is `leftOperator = lambda rightOperator`
        % at the boundary location.
        %
        % - Topic: Inspect boundary laws
        leftOperator = IMOperator()

        % Right-side eigenvalue boundary functional.
        %
        % - Topic: Inspect boundary laws
        rightOperator = IMOperator()

        % Endpoint weights implied by this law for modal inner products.
        %
        % Resolved endpoint laws store `IMBoundaryWeight` arrays. Each
        % weight selects the inner product that receives the contribution
        % and the endpoint factors used in the bilinear product.
        %
        % - Topic: Inspect boundary laws
        boundaryWeights = IMBoundaryWeight.empty(0,1)

        % True when the compatible boundary weights are known.
        %
        % A law can still assemble when this is false, but Gram
        % matrices, normalization, and index metadata may be incomplete for
        % boundary modes.
        %
        % - Topic: Inspect boundary laws
        hasKnownBoundaryWeights = true

        % Sign of the boundary-mode eigenvalue contribution.
        %
        % This sign is used when selecting a declared endpoint boundary
        % branch. Negative signs also contribute to the expected negative
        % eigenvalue index.
        %
        % - Topic: Inspect boundary laws
        indexSign = 0

        % Rank of the boundary-mode eigenvalue contribution.
        %
        % - Topic: Inspect boundary laws
        indexRank = 0

        % Physical mode number for an endpoint boundary branch.
        %
        % Boundary-mode numbers are physical labels, not eigenvalue signs.
        % The surface branch is always `-1`, the bottom branch is always
        % `-2`, and `NaN` means this condition does not declare a boundary
        % branch.
        %
        % - Topic: Inspect boundary laws
        boundaryModeNumber = NaN
    end

    properties (Access = private)
        coefficients = struct()
    end

    methods (Hidden)
        function boundary = at(self, location, options)
            % Resolve a boundary law at a physical endpoint.
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
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: boundary = at(boundary,location,options)
            % - Parameter location: `"surface"` or `"bottom"`
            % - Parameter options.formulation: EVP formulation, `"F"` or `"G"`
            % - Returns boundary: resolved boundary law
            arguments
                self IMBoundary
                location {mustBeTextScalar}
                options.formulation {mustBeTextScalar} = "formulation"
            end

            location = IMBoundary.validateLocation(location);
            if self.location ~= ""
                if self.location ~= location
                    error("IMBoundary:AlreadyPlaced", ...
                        "Boundary law is already placed at ""%s"".", self.location);
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
                    "A location-free boundary law must be placed with formulation=""F"" or formulation=""G"".");
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
                    boundaryWeight = IMBoundaryWeight(innerProduct="G", location=location, ...
                        coefficient=endpointSign, leftVariable="G", rightVariable="G");
                    switch variable
                        case "G"
                            left = IMOperator().plus(derivativeOrder=1);
                            right = IMOperator().plus(derivativeOrder=0);
                            boundary = IMBoundary(family="free", location=location, variable="G", ...
                                leftOperator=left, rightOperator=right, boundaryWeights=boundaryWeight, ...
                                indexSign=endpointSign, indexRank=1, ...
                                boundaryModeNumber=IMBoundary.boundaryModeNumberForLocation(location));
                        case "F"
                            left = IMOperator() ...
                                .plus(derivativeOrder=0) ...
                                .plus(coefficient=@(z,ctx) ctx.g./ctx.N2(z), derivativeOrder=1);
                            boundary = IMBoundary(family="free", location=location, variable="F", ...
                                leftOperator=left, boundaryWeights=boundaryWeight);
                    end
                case "dirichlet"
                    left = IMOperator().plus(derivativeOrder=0);
                    boundary = IMBoundary(family="dirichlet", location=location, variable=variable, leftOperator=left);
                case "neumann"
                    left = IMOperator().plus(derivativeOrder=1);
                    boundary = IMBoundary(family="neumann", location=location, variable=variable, leftOperator=left);
                case "custom"
                    boundaryWeights = IMBoundary.placeBoundaryWeights(self.boundaryWeights, location);
                    boundary = IMBoundary(family="custom", location=location, variable=variable, ...
                        leftOperator=self.leftOperator, rightOperator=self.rightOperator, ...
                        boundaryWeights=boundaryWeights, hasKnownBoundaryWeights=self.hasKnownBoundaryWeights, ...
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

        function tf = hasBoundaryWeights(self)
            % Return true when this law contributes boundary weights.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: tf = hasBoundaryWeights(boundary)
            % - Returns tf: true when boundary weights are present
            tf = ~isempty(self.boundaryWeights);
        end

        function count = expectedNegativeCount(self)
            % Return the negative-index contribution from this condition.
            %
            % Only conditions that declare an endpoint boundary mode
            % contribute expected index counts.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: count = expectedNegativeCount(boundary)
            % - Returns count: expected negative count contribution
            if ~self.hasKnownBoundaryWeights || isnan(self.boundaryModeNumber)
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
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: count = expectedZeroCount(boundary)
            % - Returns count: expected zero count contribution
            if ~self.hasKnownBoundaryWeights || isnan(self.boundaryModeNumber)
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
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: descriptors = boundaryModeDescriptors(boundary)
            % - Returns descriptors: structure array with `modeNumber` and `indexSign`
            descriptors = struct("modeNumber", {}, "indexSign", {});
            if ~self.hasKnownBoundaryWeights || self.indexRank == 0 || isnan(self.boundaryModeNumber)
                return;
            end
            descriptors = struct("modeNumber", self.boundaryModeNumber, "indexSign", self.indexSign);
        end
    end

    methods (Static)
        function boundary = dirichlet()
            % Create a location-free homogeneous Dirichlet boundary law.
            %
            % This constrains the EVP formulation when placed.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.dirichlet()
            % - Returns boundary: initialized boundary law
            boundary = IMBoundary(family="dirichlet");
        end

        function boundary = neumann()
            % Create a location-free homogeneous Neumann boundary law.
            %
            % This constrains the EVP formulation when placed.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.neumann()
            % - Returns boundary: initialized boundary law
            boundary = IMBoundary(family="neumann");
        end

        function boundary = custom(options)
            % Create a custom location-free operator boundary law.
            %
            % The law applies `left = lambda right` at whichever endpoint
            % the EVP factory places it. Both operators are written with the
            % physical coordinate derivative $$\partial_z$$ at the endpoint.
            %
            % `boundaryWeights` are optional endpoint contributions associated
            % with the law. Location-free weights are placed at the endpoint
            % and receive the endpoint orientation sign from Green's identity.
            % Explicitly located weights must match the endpoint where the law
            % is resolved.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.custom(options)
            % - Parameter options.left: left-side boundary functional
            % - Parameter options.right: right-side eigenvalue functional
            % - Parameter options.boundaryWeights: endpoint weights implied by this law
            % - Parameter options.hasKnownBoundaryWeights: true when the compatible boundary weights are known
            % - Parameter options.variable: variable target, `"formulation"`, `"F"`, or `"G"`
            % - Parameter options.indexSign: expected boundary-mode eigenvalue sign, `-1`, `0`, or `1`
            % - Parameter options.indexRank: number of boundary-mode directions; currently must be `1` when `boundaryModeNumber` is supplied
            % - Parameter options.boundaryModeNumber: explicit endpoint mode number, `-1` for surface or `-2` for bottom
            % - Returns boundary: initialized boundary law
            arguments
                options.left IMOperator
                options.right IMOperator = IMOperator()
                options.boundaryWeights IMBoundaryWeight = IMBoundaryWeight.empty(0,1)
                options.hasKnownBoundaryWeights (1,1) logical = true
                options.variable {mustBeTextScalar} = "formulation"
                options.indexSign (1,1) double {mustBeMember(options.indexSign, [-1 0 1])} = 0
                options.indexRank (1,1) double {mustBeInteger, mustBeNonnegative} = 0
                options.boundaryModeNumber (1,1) double = NaN
            end

            boundary = IMBoundary(family="custom", variable=options.variable, ...
                leftOperator=options.left, rightOperator=options.right, ...
                boundaryWeights=options.boundaryWeights, ...
                hasKnownBoundaryWeights=options.hasKnownBoundaryWeights, ...
                indexSign=options.indexSign, indexRank=options.indexRank, ...
                boundaryModeNumber=options.boundaryModeNumber);
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
            % EVP, it resolves to $$F+gF_z/N^2=0$$. The derivative $$G_z$$
            % or $$F_z$$ always means the coordinate derivative
            % $$\partial_z$$. A placed `G` free boundary declares an endpoint
            % boundary branch with mode number `-1` at the surface or `-2`
            % at the bottom.
            %
            % - Topic: Create boundary laws
            % - Declaration: boundary = IMBoundary.free()
            % - Returns boundary: initialized free boundary law
            boundary = IMBoundary(family="free");
        end

        function boundary = linearF(options)
            % Create a location-free linear `F` boundary law.
            %
            % When placed, the assembled boundary law represents
            % $$-(aF-bF_z/N^2)=\lambda(cF-dF_z/N^2)/g$$. Supported pairwise
            % cases also declare compatible boundary weights.
            % General unresolved coefficient patterns can still be solved,
            % but placing them warns because their endpoint contribution
            % is unknown. The derivative $$F_z$$ always means
            % $$\partial_z F$$ at either endpoint.
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
            boundary = IMBoundary(family="linearF", variable="F", coefficients=coefficients);
        end

        function boundary = linearG(options)
            % Create a location-free linear `G` boundary law.
            %
            % When placed, the assembled boundary law represents
            % $$-(eG-aG_z)=\lambda(bG-cG_z)/g$$. Supported pairwise cases
            % also declare compatible boundary weights. General
            % unresolved coefficient patterns can still be solved, but
            % placing them warns because their endpoint contribution is
            % unknown. The derivative $$G_z$$ always means
            % $$\partial_z G$$ at either endpoint.
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
            boundary = IMBoundary(family="linearG", variable="G", coefficients=coefficients);
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
                options.boundaryWeights IMBoundaryWeight = IMBoundaryWeight.empty(0,1)
                options.hasKnownBoundaryWeights (1,1) logical = true
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
            self.boundaryWeights = options.boundaryWeights(:);
            self.hasKnownBoundaryWeights = options.hasKnownBoundaryWeights;
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
            [boundaryWeights, hasKnownBoundaryWeights] = IMBoundary.linearFBoundaryWeights(location, c.a, c.b, c.c, c.d);
            boundary = IMBoundary(family="linearF", location=location, variable="F", ...
                leftOperator=left, rightOperator=right, boundaryWeights=boundaryWeights, ...
                hasKnownBoundaryWeights=hasKnownBoundaryWeights);
        end

        function boundary = placeLinearG(self, location)
            c = self.coefficients;
            left = IMOperator() ...
                .plus(coefficient=-c.e, derivativeOrder=0) ...
                .plus(coefficient=c.a, derivativeOrder=1);
            right = IMOperator() ...
                .plus(coefficient=@(~,ctx) c.b/ctx.g, derivativeOrder=0) ...
                .plus(coefficient=@(~,ctx) -c.c/ctx.g, derivativeOrder=1);
            [boundaryWeights, hasKnownBoundaryWeights] = IMBoundary.linearGBoundaryWeights(location, c.a, c.b, c.c, c.e);
            boundary = IMBoundary(family="linearG", location=location, variable="G", ...
                leftOperator=left, rightOperator=right, boundaryWeights=boundaryWeights, ...
                hasKnownBoundaryWeights=hasKnownBoundaryWeights);
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
        function [weights, hasKnownBoundaryWeights] = linearFBoundaryWeights(location, a, b, c, d)
            nonzero = abs([a b c d]) > 0;
            weights = IMBoundaryWeight.empty(0,1);
            hasKnownBoundaryWeights = true;
            if nnz(nonzero) > 2
                IMBoundary.warnUnknownBoundaryWeights("linearF");
                hasKnownBoundaryWeights = false;
                return;
            end

            endpointSign = IMBoundary.endpointSign(location);
            g = @(ctx) ctx.g;
            if ~nonzero(1) && ~nonzero(2) && nonzero(3) && nonzero(4)
                weights = IMBoundaryWeight(innerProduct="G", location=location, coefficient=endpointSign*(-d/c), ...
                    leftVariable="G", rightVariable="G");
            elseif nonzero(1) && nonzero(2) && ~nonzero(3) && ~nonzero(4)
                weights = IMBoundaryWeight(innerProduct="G", location=location, coefficient=@(ctx) endpointSign*(-b/(g(ctx)*a)), ...
                    leftVariable="G", rightVariable="G");
            elseif ~nonzero(1) && nonzero(2) && nonzero(3) && ~nonzero(4)
                weights = IMBoundaryWeight(innerProduct="F", location=location, coefficient=endpointSign*(c/b), ...
                    leftVariable="F", rightVariable="F");
            elseif nonzero(1) && ~nonzero(2) && ~nonzero(3) && nonzero(4)
                weights = [
                    IMBoundaryWeight(innerProduct="F", location=location, coefficient=@(ctx) endpointSign*(-d/(g(ctx)*g(ctx)*a)), ...
                        leftVariable="G", rightVariable="G")
                    IMBoundaryWeight(innerProduct="G", location=location, coefficient=endpointSign, ...
                        leftVariable="G", rightVariable="F")
                    IMBoundaryWeight(innerProduct="G", location=location, coefficient=endpointSign, ...
                        leftVariable="F", rightVariable="G")
                ];
            end
        end

        function [weights, hasKnownBoundaryWeights] = linearGBoundaryWeights(location, a, b, c, e)
            nonzero = abs([a b c e]) > 0;
            weights = IMBoundaryWeight.empty(0,1);
            hasKnownBoundaryWeights = true;
            if nnz(nonzero) > 2
                IMBoundary.warnUnknownBoundaryWeights("linearG");
                hasKnownBoundaryWeights = false;
                return;
            end

            endpointSign = IMBoundary.endpointSign(location);
            if ~nonzero(1) && ~nonzero(2) && nonzero(3) && nonzero(4)
                weights = IMBoundaryWeight(innerProduct="G", location=location, coefficient=endpointSign*(-c/e), ...
                    leftVariable="G", leftDerivativeOrder=1, rightVariable="G", rightDerivativeOrder=1);
            elseif nonzero(1) && nonzero(2) && ~nonzero(3) && ~nonzero(4)
                weights = IMBoundaryWeight(innerProduct="G", location=location, coefficient=endpointSign*(b/a), ...
                    leftVariable="G", rightVariable="G");
            elseif nonzero(1) && ~nonzero(2) && ~nonzero(3) && nonzero(4)
                weights = IMBoundaryWeight(innerProduct="F", location=location, coefficient=endpointSign*(-a/e), ...
                    leftVariable="F", rightVariable="F");
            elseif ~nonzero(1) && nonzero(2) && nonzero(3) && ~nonzero(4)
                weights = IMBoundaryWeight(innerProduct="F", location=location, coefficient=endpointSign*(-c/b), ...
                    leftVariable="F", rightVariable="F");
            end
        end

        function weights = placeBoundaryWeights(weights, location)
            weights = weights.at(location);
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

        function warnUnknownBoundaryWeights(family)
            warning("IMBoundary:UnknownBoundaryWeights", ...
                "%s can assemble this boundary law, but its compatible boundary weights are unknown.", ...
                family);
        end
    end
end

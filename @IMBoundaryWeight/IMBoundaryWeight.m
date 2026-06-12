classdef IMBoundaryWeight
    % Describe one endpoint contribution to a modal inner product.
    %
    % `IMBoundaryWeight` is an atomic bilinear endpoint term. The
    % `innerProduct` property names the modal inner product receiving the
    % contribution, while the left and right endpoint factors identify which
    % modal variables and physical derivative orders are evaluated.
    %
    % A resolved weight contributes
    % $$c\,X_i^{(p)}(z_\ell)Y_j^{(q)}(z_\ell)$$
    % at `location`, where `coefficient` is $$c$$, `leftVariable` is $$X$$,
    % `leftDerivativeOrder` is $$p$$, `rightVariable` is $$Y$$, and
    % `rightDerivativeOrder` is $$q$$. Location-free weights may be stored on
    % a boundary law and are placed when the law is resolved at the surface
    % or bottom endpoint.
    %
    % ```matlab
    % weight = IMBoundaryWeight(innerProduct="G", location="surface", ...
    %     coefficient=1, leftVariable="G", rightVariable="G");
    % ```
    %
    % - Topic: Create boundary weights
    % - Topic: Inspect boundary weights
    % - Topic: Developer topics
    % - Declaration: classdef IMBoundaryWeight

    properties (SetAccess = private)
        % Inner product receiving this endpoint contribution.
        %
        % The value is `"F"` or `"G"`. It selects the modal inner product
        % that receives this term; it does not restrict the endpoint factors
        % to the same variable.
        %
        % - Topic: Inspect boundary weights
        innerProduct = "G"

        % Physical endpoint location.
        %
        % Resolved weights use `"surface"` or `"bottom"`. Boundary laws may
        % store a location-free template with `location=""`; resolving that
        % template places it at the endpoint and applies the bottom
        % orientation sign when needed.
        %
        % - Topic: Inspect boundary weights
        location = ""

        % Endpoint coefficient.
        %
        % This is a scalar or a function handle evaluated with the basis-set
        % context. Function handles may accept `ctx`, `(z,ctx)`, or `z`.
        %
        % - Topic: Inspect boundary weights
        coefficient = 0

        % Variable evaluated for the left mode factor.
        %
        % - Topic: Inspect boundary weights
        leftVariable = "G"

        % Physical derivative order of the left mode factor.
        %
        % - Topic: Inspect boundary weights
        leftDerivativeOrder = 0

        % Variable evaluated for the right mode factor.
        %
        % - Topic: Inspect boundary weights
        rightVariable = "G"

        % Physical derivative order of the right mode factor.
        %
        % - Topic: Inspect boundary weights
        rightDerivativeOrder = 0
    end

    methods
        function self = IMBoundaryWeight(options)
            % Create an endpoint inner-product weight.
            %
            % - Topic: Create boundary weights
            % - Declaration: weight = IMBoundaryWeight(options)
            % - Parameter options.innerProduct: inner product receiving the term, `"F"` or `"G"`
            % - Parameter options.location: `"surface"`, `"bottom"`, or `""` for a location-free template
            % - Parameter options.coefficient: scalar or context function handle
            % - Parameter options.leftVariable: left endpoint variable, `"F"` or `"G"`
            % - Parameter options.leftDerivativeOrder: left physical derivative order
            % - Parameter options.rightVariable: right endpoint variable, `"F"` or `"G"`
            % - Parameter options.rightDerivativeOrder: right physical derivative order
            % - Returns weight: initialized boundary weight
            arguments
                options.innerProduct {mustBeTextScalar} = "G"
                options.location {mustBeTextScalar} = ""
                options.coefficient = 0
                options.leftVariable {mustBeTextScalar} = "G"
                options.leftDerivativeOrder (1,1) double {mustBeInteger, mustBeNonnegative} = 0
                options.rightVariable {mustBeTextScalar} = "G"
                options.rightDerivativeOrder (1,1) double {mustBeInteger, mustBeNonnegative} = 0
            end

            self.innerProduct = IMBoundaryWeight.validateVariable(options.innerProduct);
            self.location = IMBoundaryWeight.validateLocation(options.location, allowEmpty=true);
            self.coefficient = options.coefficient;
            self.leftVariable = IMBoundaryWeight.validateVariable(options.leftVariable);
            self.leftDerivativeOrder = options.leftDerivativeOrder;
            self.rightVariable = IMBoundaryWeight.validateVariable(options.rightVariable);
            self.rightDerivativeOrder = options.rightDerivativeOrder;
        end

        function weights = at(self, location)
            % Place location-free weights at a physical endpoint.
            %
            % Explicitly located weights must already match the requested
            % endpoint. Location-free weights receive the endpoint location;
            % bottom placement flips the coefficient sign.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: weights = at(weights,location)
            % - Parameter location: `"surface"` or `"bottom"`
            % - Returns weights: placed boundary weights
            location = IMBoundaryWeight.validateLocation(location);
            weights = self;
            for iWeight = 1:numel(weights)
                if weights(iWeight).location == ""
                    weights(iWeight).location = location;
                    if location == "bottom"
                        weights(iWeight).coefficient = IMBoundaryWeight.negateCoefficient(weights(iWeight).coefficient);
                    end
                elseif weights(iWeight).location ~= location
                    error("IMBoundaryWeight:LocationMismatch", ...
                        "A boundary weight located at ""%s"" cannot be placed at ""%s"".", ...
                        weights(iWeight).location, location);
                end
            end
        end

        function tf = isScalarEndpointValue(self, variable)
            % Return true for scalar same-variable endpoint-value weights.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: tf = isScalarEndpointValue(weight,variable)
            % - Parameter variable: variable name
            % - Returns tf: true for a scalar value-value contribution
            variable = IMBoundaryWeight.validateVariable(variable);
            tf = isnumeric(self.coefficient) && isscalar(self.coefficient) ...
                && self.leftVariable == variable && self.rightVariable == variable ...
                && self.leftDerivativeOrder == 0 && self.rightDerivativeOrder == 0;
        end
    end

    methods (Static, Access = private)
        function variable = validateVariable(variable)
            variable = string(variable);
            if variable ~= "F" && variable ~= "G"
                error("IMBoundaryWeight:InvalidVariable", ...
                    "variable must be ""F"" or ""G"".");
            end
        end

        function location = validateLocation(location, options)
            arguments
                location {mustBeTextScalar}
                options.allowEmpty (1,1) logical = false
            end

            location = string(location);
            if options.allowEmpty && location == ""
                return;
            end
            if location ~= "surface" && location ~= "bottom"
                error("IMBoundaryWeight:InvalidLocation", ...
                    "location must be ""surface"" or ""bottom"".");
            end
        end

        function coefficient = negateCoefficient(coefficient)
            if isa(coefficient, "function_handle")
                coefficient = @(varargin) -coefficient(varargin{:});
            else
                coefficient = -coefficient;
            end
        end
    end
end

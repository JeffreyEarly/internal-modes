classdef IMBoundaryRow
    % Store a resolved endpoint row for a v2 internal-mode EVP.
    %
    % `IMBoundaryRow` is the placed form of an `IMBoundary` law. It owns the
    % matrix row contribution, endpoint trace terms used in inner products,
    % and boundary index metadata for one endpoint.
    %
    % ```matlab
    % row = IMBoundary.rigid().resolve(endpoint="upper", context=struct("primaryComponent","G"));
    % ```
    %
    % - Topic: Create boundary rows
    % - Topic: Inspect boundary rows
    % - Topic: Apply boundary rows
    % - Topic: Developer topics
    % - Declaration: classdef IMBoundaryRow

    properties
        % Boundary family name.
        %
        % - Topic: Inspect boundary rows
        family = "custom"

        % Resolved endpoint name.
        %
        % Supported values are `"upper"` and `"lower"`.
        %
        % - Topic: Inspect boundary rows
        endpoint = "upper"

        % Component constrained by this boundary row.
        %
        % - Topic: Inspect boundary rows
        component = ""

        % Left-side boundary functional.
        %
        % - Topic: Apply boundary rows
        leftOperator = IMOperator.strong()

        % Right-side eigenvalue boundary functional.
        %
        % - Topic: Apply boundary rows
        rightOperator = IMOperator.strong()

        % Endpoint trace terms implied by this boundary row.
        %
        % Each term has an inner-product component, physical location,
        % coefficient, left trace, and right trace.
        %
        % - Topic: Inspect boundary rows
        endpointTerms = IMBoundaryRow.emptyEndpointTerms()

        % Sign of the boundary index contribution.
        %
        % Negative signs contribute to the expected negative index.
        %
        % - Topic: Inspect boundary rows
        indexSign = 0

        % Rank of the boundary index contribution.
        %
        % - Topic: Inspect boundary rows
        indexRank = 0

        % Status of the associated orthogonality metadata.
        %
        % Values are `"complete"` or `"unresolved"`.
        %
        % - Topic: Inspect boundary rows
        orthogonalityStatus = "complete"
    end

    methods
        function self = IMBoundaryRow(options)
            % Create a resolved boundary row.
            %
            % - Topic: Create boundary rows
            % - Declaration: row = IMBoundaryRow(options)
            % - Parameter options.family: boundary family name
            % - Parameter options.endpoint: `"upper"` or `"lower"`
            % - Parameter options.component: constrained component
            % - Parameter options.leftOperator: left-side boundary functional
            % - Parameter options.rightOperator: right-side eigenvalue functional
            % - Parameter options.endpointTerms: endpoint trace-pair terms
            % - Parameter options.indexSign: index contribution sign
            % - Parameter options.indexRank: index contribution rank
            % - Parameter options.orthogonalityStatus: endpoint-form status
            % - Returns row: initialized boundary row
            arguments
                options.family {mustBeTextScalar} = "custom"
                options.endpoint {mustBeTextScalar} = "upper"
                options.component {mustBeTextScalar} = ""
                options.leftOperator IMOperator = IMOperator.strong()
                options.rightOperator IMOperator = IMOperator.strong()
                options.endpointTerms struct = IMBoundaryRow.emptyEndpointTerms()
                options.indexSign (1,1) double = 0
                options.indexRank (1,1) double {mustBeInteger, mustBeNonnegative} = 0
                options.orthogonalityStatus {mustBeTextScalar} = "complete"
            end

            self.family = string(options.family);
            self.endpoint = IMBoundaryRow.validateEndpoint(options.endpoint);
            self.component = string(options.component);
            self.leftOperator = options.leftOperator;
            self.rightOperator = options.rightOperator;
            self.endpointTerms = options.endpointTerms(:);
            self.indexSign = sign(options.indexSign);
            self.indexRank = options.indexRank;
            self.orthogonalityStatus = string(options.orthogonalityStatus);
        end

        function [A, B] = apply(self, A, B, solver)
            % Apply the boundary row to a matrix pair.
            %
            % Active metadata-only rows do not replace matrix rows.
            %
            % - Topic: Apply boundary rows
            % - Declaration: [A,B] = apply(row,A,B,solver)
            % - Parameter A: left EVP matrix
            % - Parameter B: right EVP matrix
            % - Parameter solver: coordinate-aware internal-mode solver
            % - Returns A: boundary-conditioned left matrix
            % - Returns B: boundary-conditioned right matrix
            if self.family == "active" || self.family == "partialDepthPE"
                return;
            end

            location = IMBoundaryRow.locationForEndpoint(self.endpoint);
            index = solver.boundaryIndex(location);
            A(index,:) = self.leftOperator.boundaryRow(solver, location);
            B(index,:) = self.rightOperator.boundaryRow(solver, location);
        end

        function tf = hasEndpointTerms(self)
            % Return true when this row contributes endpoint terms.
            %
            % - Topic: Inspect boundary rows
            % - Declaration: tf = hasEndpointTerms(row)
            % - Returns tf: true when endpoint terms are present
            tf = ~isempty(self.endpointTerms);
        end

        function count = expectedNegativeCount(self)
            % Return the negative-index contribution from this row.
            %
            % - Topic: Inspect boundary rows
            % - Declaration: count = expectedNegativeCount(row)
            % - Returns count: expected negative count contribution
            if self.orthogonalityStatus == "unresolved"
                count = 0;
            elseif self.indexSign < 0
                count = self.indexRank;
            else
                count = 0;
            end
        end

        function count = expectedZeroCount(self)
            % Return the zero-index contribution from this row.
            %
            % - Topic: Inspect boundary rows
            % - Declaration: count = expectedZeroCount(row)
            % - Returns count: expected zero count contribution
            if self.orthogonalityStatus == "unresolved"
                count = 0;
            elseif self.indexSign == 0
                count = self.indexRank;
            else
                count = 0;
            end
        end
    end

    methods (Static)
        function row = active(options)
            % Create an active metadata-only boundary row.
            %
            % Active rows are already placed because their trace terms refer
            % to an endpoint of a partial-depth interval.
            %
            % - Topic: Create boundary rows
            % - Declaration: row = IMBoundaryRow.active(options)
            % - Parameter options.endpoint: `"upper"` or `"lower"`
            % - Parameter options.component: active component
            % - Parameter options.indexSign: active-boundary sign
            % - Parameter options.indexRank: number of active directions
            % - Parameter options.endpointTerms: endpoint trace-pair terms
            % - Returns row: initialized active boundary row
            arguments
                options.endpoint {mustBeTextScalar}
                options.component {mustBeTextScalar} = "G"
                options.indexSign (1,1) double {mustBeMember(options.indexSign, [-1 1])}
                options.indexRank (1,1) double {mustBeInteger, mustBePositive} = 1
                options.endpointTerms struct = IMBoundaryRow.emptyEndpointTerms()
            end

            row = IMBoundaryRow(family="active", endpoint=options.endpoint, component=options.component, ...
                endpointTerms=options.endpointTerms, indexSign=options.indexSign, indexRank=options.indexRank);
        end

        function rows = partialDepthPE(options)
            % Create partial-depth potential-energy active boundary rows.
            %
            % Positive boundary signs add no negative index directions.
            % Negative boundary signs add one negative direction at each
            % window endpoint.
            %
            % - Topic: Create boundary rows
            % - Declaration: rows = IMBoundaryRow.partialDepthPE(options)
            % - Parameter options.boundarySign: `"positive"` or `"negative"`
            % - Returns rows: lower and upper active boundary rows
            arguments
                options.boundarySign {mustBeTextScalar} = "positive"
            end

            switch string(options.boundarySign)
                case "positive"
                    indexSign = 1;
                case "negative"
                    indexSign = -1;
                otherwise
                    error("IMBoundaryRow:InvalidBoundarySign", ...
                        "boundarySign must be ""positive"" or ""negative"".");
            end
            rows = [
                IMBoundaryRow(family="partialDepthPE", endpoint="lower", indexSign=indexSign, indexRank=1)
                IMBoundaryRow(family="partialDepthPE", endpoint="upper", indexSign=indexSign, indexRank=1)
            ];
        end

        function trace = trace(component, options)
            % Create an endpoint trace descriptor.
            %
            % A trace describes which component value or first derivative is
            % evaluated at a boundary endpoint.
            %
            % - Topic: Create boundary rows
            % - Declaration: trace = IMBoundaryRow.trace(component,options)
            % - Parameter component: component name
            % - Parameter options.derivativeOrder: physical derivative order
            % - Returns trace: endpoint trace descriptor
            arguments
                component {mustBeTextScalar}
                options.derivativeOrder (1,1) double {mustBeInteger, mustBeNonnegative} = 0
            end

            trace = struct("component", string(component), "derivativeOrder", options.derivativeOrder);
        end

        function term = endpointTerm(innerProductComponent, location, coefficient, leftTrace, rightTrace)
            % Create an endpoint trace-pair term.
            %
            % The term contributes `coefficient*leftTrace_i*rightTrace_j`
            % to the named component's inner product.
            %
            % - Topic: Create boundary rows
            % - Declaration: term = IMBoundaryRow.endpointTerm(innerProductComponent,location,coefficient,leftTrace,rightTrace)
            % - Parameter innerProductComponent: component whose inner product receives the term
            % - Parameter location: boundary location
            % - Parameter coefficient: scalar or context function handle
            % - Parameter leftTrace: trace evaluated for the left mode
            % - Parameter rightTrace: trace evaluated for the right mode
            % - Returns term: endpoint trace-pair term
            term = struct("innerProductComponent", string(innerProductComponent), ...
                "location", string(location), "coefficient", coefficient, ...
                "leftTrace", leftTrace, "rightTrace", rightTrace);
        end

        function terms = emptyEndpointTerms()
            % Return an empty endpoint-term structure.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: terms = IMBoundaryRow.emptyEndpointTerms()
            % - Returns terms: empty endpoint-term structure
            terms = struct("innerProductComponent", {}, "location", {}, "coefficient", {}, ...
                "leftTrace", {}, "rightTrace", {});
        end

        function location = locationForEndpoint(endpoint)
            % Return the physical location associated with an endpoint.
            %
            % - Topic: Developer topics
            % - Developer: true
            % - Declaration: location = IMBoundaryRow.locationForEndpoint(endpoint)
            % - Parameter endpoint: `"upper"` or `"lower"`
            % - Returns location: `"surface"` or `"bottom"`
            endpoint = IMBoundaryRow.validateEndpoint(endpoint);
            switch endpoint
                case "upper"
                    location = "surface";
                case "lower"
                    location = "bottom";
            end
        end
    end

    methods (Static, Access = private)
        function endpoint = validateEndpoint(endpoint)
            endpoint = string(endpoint);
            if endpoint ~= "upper" && endpoint ~= "lower"
                error("IMBoundaryRow:InvalidEndpoint", ...
                    "endpoint must be ""upper"" or ""lower"".");
            end
        end
    end
end

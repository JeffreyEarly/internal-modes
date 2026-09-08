classdef IMProjectionRecipe
    % Bind continuous basis mathematics to fixed-grid numerical projection.
    %
    % The signed sample metric defines projection; the positive majorant
    % defines error magnitudes. This value object snapshots its basis and
    % normalization. It does not choose a mode count or certify solve accuracy.
    %
    % ```matlab
    % recipe = basis.projectionRecipe(variable="F");
    % projection = recipe.projection(z,weights);
    % ```
    %
    % - Topic: Create projection recipes
    % - Topic: Inspect projection recipes
    % - Topic: Evaluate projection recipes
    % - Declaration: classdef IMProjectionRecipe
    properties (SetAccess = private)
        % Whether the continuous projection pairing is available.
        % - Topic: Inspect projection recipes
        available
        % Explanation when the continuous pairing is unavailable.
        % - Topic: Inspect projection recipes
        reason
        % Scalar variable evaluated by this recipe.
        % - Topic: Inspect projection recipes
        variable
        % Scientific labels in basis-column order.
        % - Topic: Inspect projection recipes
        columnLabels
        % Snapshotted normalization convention.
        % - Topic: Inspect projection recipes
        normalization
        % Continuous signed target Gram matrix.
        % - Topic: Inspect projection recipes
        targetGramMatrix
        % Continuous positive-majorant Gram matrix.
        % - Topic: Inspect projection recipes
        majorantGramMatrix
        % Whether the family supports leakage measurements.
        % - Topic: Inspect projection recipes
        supportsLeakage
        % Whether the family supports quadratic measurements.
        % - Topic: Inspect projection recipes
        supportsQuadratic
        % Representation and reference-integration provenance.
        % - Topic: Inspect projection recipes
        provenance
    end
    properties (Access = private)
        zDomain
        evaluateFunction
        weightFunction
        spec
    end
    methods
        function self = IMProjectionRecipe(options)
            % Store a basis-owned recipe without changing its scientific state.
            % - Topic: Create projection recipes
            % - Parameter options.available: continuous pairing availability
            % - Parameter options.reason: unavailable pairing explanation
            % - Parameter options.variable: sampled scalar variable
            % - Parameter options.columnLabels: scientific labels in column order
            % - Parameter options.normalization: frozen normalization convention
            % - Parameter options.targetGramMatrix: continuous signed Gram target
            % - Parameter options.majorantGramMatrix: positive coefficient metric
            % - Parameter options.supportsLeakage: family leakage capability
            % - Parameter options.supportsQuadratic: family quadratic capability
            % - Parameter options.provenance: representation and reference quality
            % - Parameter options.zDomain: physical bounds
            % - Parameter options.evaluateFunction: normalized basis evaluator
            % - Parameter options.weightFunction: bound interior-weight evaluator
            % - Parameter options.spec: basis-owned endpoint pairing specification
            % - Returns self: immutable recipe
            arguments (Input)
                options.available (1,1) logical = true
                options.reason (1,1) string = ""
                options.variable (1,1) string = "u"
                options.columnLabels (1,:) string = strings(1,0)
                options.normalization (1,1) string = ""
                options.targetGramMatrix (:,:) double {mustBeFinite} = zeros(0)
                options.majorantGramMatrix (:,:) double {mustBeFinite} = zeros(0)
                options.supportsLeakage (1,1) logical = false
                options.supportsQuadratic (1,1) logical = false
                options.provenance (1,1) struct = struct()
                options.zDomain (1,2) double {mustBeReal,mustBeFinite} = [0 1]
                options.evaluateFunction (1,1) function_handle = @(z) zeros(numel(z),0)
                options.weightFunction (1,1) function_handle = @(z) zeros(size(z))
                options.spec (1,1) struct = struct()
            end
            arguments (Output)
                self (1,1) IMProjectionRecipe
            end
            fields = fieldnames(options);
            for i = 1:numel(fields)
                self.(fields{i}) = options.(fields{i});
            end
            n = numel(self.columnLabels);
            if self.available && (~isequal(size(self.targetGramMatrix),[n n]) || ~isequal(size(self.majorantGramMatrix),[n n]))
                error("IMProjectionRecipe:InvalidTargetShape","Both target matrices must have one row and column per scientific label.");
            end
            if self.zDomain(1) >= self.zDomain(2)
                error("IMProjectionRecipe:InvalidDomain","zDomain must be increasing.");
            end
        end
        function values = evaluate(self,z,options)
            % Evaluate normalized continuous columns at physical points.
            % - Topic: Evaluate projection recipes
            % - Parameter z: column vector of physical coordinates
            % - Parameter options.columns: selected column indices, all by default
            % - Returns values: nZ by nSelectedColumns values
            arguments (Input)
                self IMProjectionRecipe
                z (:,1) double {mustBeReal,mustBeFinite}
                options.columns (1,:) double {mustBeInteger,mustBePositive} = []
            end
            arguments (Output)
                values (:,:) double
            end
            self.assertAvailable();
            self.validatePoints(z);
            columns = self.selectedColumns(options.columns);
            values = self.evaluateFunction(z);
            values = values(:,columns);
        end
        function matrix = metric(self,z,weights,options)
            % Construct the signed or majorant pairing on fixed samples.
            %
            % Endpoint traces must be represented by supplied variable
            % samples. Derivative and companion-variable traces are explicit
            % unsupported capabilities, never inferred from the basis span.
            % Quadrature weights are unchanged; negative quadrature weights
            % can make the sampled majorant indefinite.
            % - Topic: Evaluate projection recipes
            % - Parameter z: strictly increasing physical sample coordinates
            % - Parameter weights: fixed quadrature weights
            % - Parameter options.majorant: use absolute endpoint coefficients
            % - Returns matrix: nZ by nZ sample metric
            arguments (Input)
                self IMProjectionRecipe
                z (:,1) double {mustBeReal,mustBeFinite}
                weights (:,1) double {mustBeReal,mustBeFinite}
                options.majorant (1,1) logical = false
            end
            arguments (Output)
                matrix (:,:) double
            end
            self.assertAvailable();
            self.validatePoints(z);
            if numel(z) < 2 || any(diff(z) <= 0)
                error("IMProjectionRecipe:InvalidGrid","z must contain at least two strictly increasing points.");
            end
            if numel(weights) ~= numel(z) || ~any(weights ~= 0)
                error("IMProjectionRecipe:InvalidWeights","Supply one weight per sample and at least one nonzero weight.");
            end
            weight = self.weightFunction(z);
            if isscalar(weight)
                weight = repmat(weight,size(z));
            end
            if numel(weight) ~= numel(z) || ~isreal(weight) || any(~isfinite(weight(:)))
                error("IMProjectionRecipe:InvalidInteriorWeight","The interior weight must return one finite real value per sample.");
            end
            if options.majorant
                weight = abs(weight);
            end
            matrix = diag(weights.*weight(:));
            endpoints = [self.spec.surfaceWeights;self.spec.bottomWeights];
            for i = 1:numel(endpoints)
                term = endpoints(i);
                if term.d ~= 0
                    error("IMProjectionRecipe:UnsupportedDerivativeTrace","The %s endpoint requires a derivative trace; value-only samples cannot represent this projection.",term.location);
                end
                if term.c ~= 0
                    matrix = self.addEndpoint(matrix,z,term.location,term.coefficient*term.c^2,options.majorant);
                end
            end
            if isfield(self.spec,"endpointInnerProductTerms")
                for i = 1:numel(self.spec.endpointInnerProductTerms)
                    term = self.spec.endpointInnerProductTerms(i);
                    if string(term.variable) ~= self.variable
                        error("IMProjectionRecipe:UnsupportedCompanionTrace","The %s endpoint requires %s values in the %s pairing; paired-state projection is unsupported.",term.location,term.variable,self.variable);
                    end
                    matrix = self.addEndpoint(matrix,z,term.location,term.coefficient,options.majorant);
                end
            end
        end
        function projection = projection(self,z,weights,options)
            % Construct a projection for explicit columns on an unchanged grid.
            % - Topic: Evaluate projection recipes
            % - Parameter z: physical sample coordinates
            % - Parameter weights: fixed quadrature weights
            % - Parameter options.columns: explicit column selection, all by default
            % - Returns projection: numerical projection and measurements
            arguments (Input)
                self IMProjectionRecipe
                z (:,1) double {mustBeReal,mustBeFinite}
                weights (:,1) double {mustBeReal,mustBeFinite}
                options.columns (1,:) double {mustBeInteger,mustBePositive} = []
            end
            arguments (Output)
                projection (1,1) IMProjection
            end
            columns = self.selectedColumns(options.columns);
            sampled = self.evaluate(z,columns=columns);
            matrix = self.metric(z,weights);
            target = self.targetGramMatrix(columns,columns);
            majorant = self.majorantGramMatrix(columns,columns);
            zeroTarget = abs(diag(target).') <= 1e3*eps(max(1,diag(majorant).'));
            columnNorms = vecnorm(sampled,2,1);
            zeroSample = columnNorms <= 1e3*eps(max(1,columnNorms));
            if any(zeroTarget & ~zeroSample)
                labels = self.columnLabels(columns);
                error("IMProjectionRecipe:UnsupportedZeroNormColumn","Nonzero sampled %s columns have zero continuous norm for label(s) %s.",self.variable,join(labels(zeroTarget & ~zeroSample),", "));
            end
            projection = IMProjection(sampled,matrix,target,majorantGramMatrix=majorant,activeColumnMask=~(zeroTarget & zeroSample),columnLabels=self.columnLabels(columns));
        end
    end
    methods (Access = private)
        function assertAvailable(self)
            if ~self.available
                error("IMProjectionRecipe:UnavailableProjection","%s",self.reason);
            end
        end
        function columns = selectedColumns(self,columns)
            if isempty(columns)
                columns = 1:numel(self.columnLabels);
            elseif any(columns > numel(self.columnLabels)) || numel(unique(columns)) ~= numel(columns)
                error("IMProjectionRecipe:InvalidColumns","Select distinct existing column indices; scientific labels are not array indices.");
            end
        end
        function validatePoints(self,z)
            tolerance = 100*eps(max(1,max(abs(self.zDomain))));
            if any(z < self.zDomain(1)-tolerance | z > self.zDomain(2)+tolerance)
                error("IMProjectionRecipe:InvalidGrid","All points must lie inside the basis domain.");
            end
        end
        function matrix = addEndpoint(self,matrix,z,location,coefficient,majorant)
            if ~isfinite(coefficient)
                error("IMProjectionRecipe:InvalidEndpointCoefficient","The %s endpoint coefficient must be finite.",location);
            end
            if string(location) == "surface"
                endpoint = self.zDomain(2);
            else
                endpoint = self.zDomain(1);
            end
            index = find(abs(z-endpoint) <= 100*eps(max(1,max(abs(self.zDomain)))),1);
            if isempty(index)
                error("IMProjectionRecipe:MissingEndpointSample","Include the %s endpoint at z=%g to represent its value-only pairing.",location,endpoint);
            end
            if majorant
                coefficient = abs(coefficient);
            end
            matrix(index,index) = matrix(index,index)+coefficient;
        end
    end
end

classdef IMVerticalSampling
    % Describe vertical point, cell-average and cell-integral observations.
    %
    % Each immutable table row is one declared observation of an existing
    % continuous vertical profile. Locations may differ between variables
    % and need not include the physical endpoints. No grid, physical field
    % composition, column activity or fitting/noise metric is inferred.
    %
    % ```matlab
    % points = IMVerticalSampling.pointValues("F",z,observationLabels=labels);
    % cells = IMVerticalSampling.cellAverages("G",bounds);
    % sampling = IMVerticalSampling([points.observations;cells.observations]);
    % [A,record] = sampling.evaluate(basis);
    % ```
    %
    % - Topic: Create vertical sampling
    % - Topic: Inspect vertical sampling
    % - Topic: Evaluate vertical sampling
    % - Declaration: classdef IMVerticalSampling

    properties (SetAccess = private)
        % Ordered observation definitions and caller-declared units.
        %
        % Required columns are observationLabel, variable, kind, zLower and
        % zUpper. derivativeOrder defaults to zero and units to unspecified.
        % Point rows have equal bounds; cell rows have increasing bounds.
        % Additional table columns are preserved as caller metadata.
        % - Topic: Inspect vertical sampling
        observations table
    end

    properties (Dependent, SetAccess = private)
        % Unique row identities in the exact declared observation order.
        % - Topic: Inspect vertical sampling
        observationLabels
        % Number of declared observations.
        % - Topic: Inspect vertical sampling
        observationCount
    end

    methods
        function self = IMVerticalSampling(observations)
            % Validate an ordered table of vertical observation functionals.
            %
            % Profile variables are F, G, APV and their supported aliases.
            % Availability and physical-domain containment are checked when
            % evaluating a particular basis. Geometry may repeat and rows
            % may appear in any order; their identities must remain unique.
            %
            % - Topic: Create vertical sampling
            % - Parameter observations: table of labelled point/cell observations
            % - Returns self: immutable vertical sampling definition
            arguments (Input)
                observations table
            end
            arguments (Output)
                self (1,1) IMVerticalSampling
            end
            required = ["observationLabel","variable","kind","zLower","zUpper"];
            if ~all(ismember(required,string(observations.Properties.VariableNames)))
                error("IMVerticalSampling:InvalidObservations","Supply observationLabel, variable, kind, zLower and zUpper table columns.");
            end
            count = height(observations);
            if ~ismember("derivativeOrder",string(observations.Properties.VariableNames))
                observations.derivativeOrder = zeros(count,1);
            end
            if ~ismember("units",string(observations.Properties.VariableNames))
                observations.units = repmat("unspecified",count,1);
            end
            for name = ["observationLabel","variable","kind","units"]
                value = string(observations.(name));
                if ~isequal(size(value),[count 1]) || any(ismissing(value))
                    error("IMVerticalSampling:InvalidObservations","%s must contain one nonmissing string per observation.",name);
                end
                observations.(name) = value;
            end
            if any(strlength(observations.observationLabel) == 0) || numel(unique(observations.observationLabel)) ~= count
                error("IMVerticalSampling:InvalidObservationLabels","Observation labels must be nonempty and distinct, including for repeated geometry.");
            end
            if any(~ismember(observations.variable,["F","Fz","G","Gz","Gzz","APV"]))
                error("IMVerticalSampling:InvalidVariable","Use vertical F, Fz, G, Gz, Gzz or APV profiles; physical field composition belongs to the caller.");
            end
            if any(~ismember(observations.kind,["point","cellAverage","cellIntegral"]))
                error("IMVerticalSampling:InvalidKind","Observation kind must be point, cellAverage or cellIntegral.");
            end
            for name = ["zLower","zUpper","derivativeOrder"]
                value = observations.(name);
                if ~isa(value,"double") || ~isequal(size(value),[count 1]) || ~isreal(value) || any(~isfinite(value))
                    error("IMVerticalSampling:InvalidObservations","%s must contain one finite real double per observation.",name);
                end
            end
            if any(observations.derivativeOrder < 0 | observations.derivativeOrder ~= fix(observations.derivativeOrder))
                error("IMVerticalSampling:InvalidDerivativeOrder","derivativeOrder must contain nonnegative integers.");
            end
            points = observations.kind == "point";
            widths = observations.zUpper-observations.zLower;
            midpoints = observations.zLower+widths/2;
            if any(widths(points) ~= 0) || any(widths(~points) <= 0 | ~isfinite(widths(~points))) || any(midpoints(~points) <= observations.zLower(~points) | midpoints(~points) >= observations.zUpper(~points))
                error("IMVerticalSampling:InvalidSupport","Point bounds must coincide; cells need a finite positive width and a representable interior midpoint.");
            end
            self.observations = observations;
        end

        function labels = get.observationLabels(self)
            labels = self.observations.observationLabel.';
        end

        function count = get.observationCount(self)
            count = height(self.observations);
        end

        [synthesis,record] = evaluate(self,basis,options)
    end

    methods (Static)
        function sampling = pointValues(variable,z,options)
            % Describe values or supported derivatives at supplied depths.
            % - Topic: Create vertical sampling
            % - Parameter variable: vertical profile or derivative alias
            % - Parameter z: physical point column in caller order
            % - Parameter options.observationLabels: optional distinct row labels
            % - Parameter options.derivativeOrder: additional derivative order
            % - Parameter options.units: caller-declared output units
            % - Returns sampling: point observations without required endpoints
            arguments (Input)
                variable (1,1) string
                z (:,1) double {mustBeReal,mustBeFinite}
                options.observationLabels (1,:) string = strings(1,0)
                options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 0
                options.units (1,1) string = "unspecified"
            end
            arguments (Output)
                sampling (1,1) IMVerticalSampling
            end
            sampling = IMVerticalSampling.build(variable,[z z],"point",options);
        end

        function sampling = cellAverages(variable,bounds,options)
            % Describe physical-depth cell means, not centre-point values.
            %
            % Each row means the integral over dz divided by cell width.
            % - Topic: Create vertical sampling
            % - Parameter variable: vertical profile or derivative alias
            % - Parameter bounds: one increasing physical-depth pair per cell
            % - Parameter options.observationLabels: optional distinct row labels
            % - Parameter options.derivativeOrder: additional derivative order
            % - Parameter options.units: caller-declared output units
            % - Returns sampling: physical-depth cell-average observations
            arguments (Input)
                variable (1,1) string
                bounds (:,2) double {mustBeReal,mustBeFinite}
                options.observationLabels (1,:) string = strings(1,0)
                options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 0
                options.units (1,1) string = "unspecified"
            end
            arguments (Output)
                sampling (1,1) IMVerticalSampling
            end
            sampling = IMVerticalSampling.build(variable,bounds,"cellAverage",options);
        end

        function sampling = cellIntegrals(variable,bounds,options)
            % Describe vertical cell integrals with the physical dz measure.
            % - Topic: Create vertical sampling
            % - Parameter variable: vertical profile or derivative alias
            % - Parameter bounds: one increasing physical-depth pair per cell
            % - Parameter options.observationLabels: optional distinct row labels
            % - Parameter options.derivativeOrder: additional derivative order
            % - Parameter options.units: caller-declared output units, including the depth integral
            % - Returns sampling: physical-depth cell-integral observations
            arguments (Input)
                variable (1,1) string
                bounds (:,2) double {mustBeReal,mustBeFinite}
                options.observationLabels (1,:) string = strings(1,0)
                options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 0
                options.units (1,1) string = "unspecified"
            end
            arguments (Output)
                sampling (1,1) IMVerticalSampling
            end
            sampling = IMVerticalSampling.build(variable,bounds,"cellIntegral",options);
        end
    end

    methods (Static, Access = private)
        function sampling = build(variable,bounds,kind,options)
            count = size(bounds,1);
            labels = options.observationLabels;
            if isempty(labels)
                labels = variable+":"+kind+":"+string(1:count);
            elseif numel(labels) ~= count
                error("IMVerticalSampling:InvalidObservationLabels","Supply one observation label per point or cell.");
            end
            observations = table(labels(:),repmat(variable,count,1),repmat(options.derivativeOrder,count,1),repmat(kind,count,1),bounds(:,1),bounds(:,2),repmat(options.units,count,1),VariableNames=["observationLabel","variable","derivativeOrder","kind","zLower","zUpper","units"]);
            sampling = IMVerticalSampling(observations);
        end
    end
end

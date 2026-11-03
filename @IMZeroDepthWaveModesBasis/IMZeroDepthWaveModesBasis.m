classdef IMZeroDepthWaveModesBasis
    % Evaluate a finite auxiliary basis of exact zero-depth waves.
    %
    % The caller declares maximal critical intervals where N2=f0^2 exactly.
    % DD sine profiles G vanish outside these intervals. The manuscript
    % convention F=D Gz is imposed directly, including one-sided interval
    % endpoint traces. Every column has h=0 and frequency abs(f0); there is
    % no finite inverse-depth eigenvalue. Counts truncate each interval's
    % infinite-dimensional kernel explicitly, with no extra constant modes.
    %
    % ```matlab
    % basis = IMZeroDepthWaveModesBasis(background,k,criticalIntervals=[-200 0],nModes=3);
    % F = basis.F(z,normalization="depthRMSF");
    % ```
    %
    % - Topic: Create zero-depth waves
    % - Topic: Inspect zero-depth waves
    % - Topic: Evaluate zero-depth waves
    % - Declaration: classdef IMZeroDepthWaveModesBasis

    properties (SetAccess = private)
        % Authoritative physical column.
        % - Topic: Inspect zero-depth waves
        background
        % Positive horizontal wavenumber, in inverse meters.
        % - Topic: Inspect zero-depth waves
        k (1,1) double
        % Declared exact critical intervals.
        % - Topic: Inspect zero-depth waves
        criticalIntervals (:,2) double
    end

    properties
        % Active profile normalization, applied lazily to aligned F and G.
        % - Topic: Evaluate zero-depth waves
        normalization (1,1) string = "waveEnergy"
    end

    properties (Access = private)
        auxiliary
        normalizationRules = configureDictionary("string","cell")
    end

    properties (Dependent, SetAccess = private)
        % Retained column count.
        % - Topic: Inspect zero-depth waves
        nColumns
        % Physical bottom and surface coordinates.
        % - Topic: Inspect zero-depth waves
        zDomain
        % Scientific family identifier.
        % - Topic: Inspect zero-depth waves
        family
        % Distinct interval and local-mode column identities.
        % - Topic: Inspect zero-depth waves
        columnLabels
        % Branch and interval metadata in retained order.
        % - Topic: Inspect zero-depth waves
        modeInfo
        % Exact zero equivalent depths, in meters.
        % - Topic: Inspect zero-depth waves
        h
        % Physical inertial frequencies, in inverse seconds.
        % - Topic: Inspect zero-depth waves
        frequency
        % Supported profile and derivative channels.
        % - Topic: Inspect zero-depth waves
        capabilities
    end

    methods
        function self = IMZeroDepthWaveModesBasis(background,k,options)
            % Construct a truncated exact kernel from declared critical intervals.
            %
            % nModes is scalar or contains one exact count per interval.
            % A zero count retains no columns on that interval.
            %
            % - Topic: Create zero-depth waves
            % - Parameter background: authoritative physical column
            % - Parameter k: positive horizontal wavenumber
            % - Parameter options.criticalIntervals: exact N2=f0^2 intervals
            % - Parameter options.nModes: nonnegative count per interval
            % - Parameter options.normalization: active normalization
            % - Returns self: capable zero-depth wave basis
            arguments
                background (1,1) IMBackground
                k (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
                options.criticalIntervals (:,2) double {mustBeReal,mustBeFinite}
                options.nModes (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 8
                options.normalization (1,1) string = "waveEnergy"
            end
            IMDeclaredIntervals.validate(background,options.criticalIntervals,background.f0^2,"critical","IMZeroDepthWaveModesBasis");
            self.background = background;
            self.k = k;
            self.criticalIntervals = options.criticalIntervals;
            self.auxiliary = IMNeutralLabelModesBasis(options.criticalIntervals,zDomain=background.zDomain,nModes=options.nModes,rotationRate=background.rotationRate,latitude=background.latitude);
            self.normalization = options.normalization;
        end

        function value = get.nColumns(self)
            value = self.auxiliary.nColumns;
        end
        function value = get.zDomain(self)
            value = self.background.zDomain;
        end
        function value = get.family(~)
            value = "zeroDepthWave";
        end
        function value = get.columnLabels(self)
            value = "zeroDepth:"+self.auxiliary.columnLabels;
        end
        function info = get.modeInfo(self)
            info = self.auxiliary.modeInfo;
            info.columnLabel = self.columnLabels(:);
            info.branch(:) = "zeroDepth";
        end
        function value = get.h(self)
            value = zeros(1,self.nColumns);
        end
        function value = get.frequency(self)
            value = repmat(abs(self.background.f0),1,self.nColumns);
        end
        function value = get.capabilities(~)
            value = ["F","Fz","G","Gz","Gzz"];
        end

        function self = selectColumns(self,columns)
            % Select an aligned subset without solving or resampling.
            % - Topic: Inspect zero-depth waves
            source = self;
            for name = reshape(string(keys(self.normalizationRules)),1,[])
                self.normalizationRules{name} = @(~,j) IMInternalModesNormalizationTools.selectedNormFactor(source,name,columns(j));
            end
            self.auxiliary = self.auxiliary.selectColumns(columns);
        end

        function self = addNormalization(self,name,rule)
            % Install factor = rule(basis,j), dividing both F and G.
            % - Topic: Evaluate zero-depth waves
            arguments
                self IMZeroDepthWaveModesBasis
                name (1,1) string
                rule (1,1) function_handle
            end
            self.normalizationRules{name} = rule;
        end

        function names = normalizationNames(self)
            % List exact and custom profile normalization rules.
            % - Topic: Evaluate zero-depth waves
            names = unique(["unity";"depthRMSF";"depthRMSG";"maxAbsF";"maxAbsG";"waveEnergy";reshape(string(keys(self.normalizationRules)),[],1)]);
        end

        function factors = normalizationFactors(self,normalization)
            % Return exact analytical factors without profile quadrature.
            % - Topic: Evaluate zero-depth waves
            arguments
                self IMZeroDepthWaveModesBasis
                normalization (1,1) string = self.normalization
            end
            D = self.background.D;
            m = self.auxiliary.m_n;
            switch normalization
                case {"unity","depthRMSG"}
                    factors = ones(1,self.nColumns);
                case "depthRMSF"
                    factors = D*m;
                case "waveEnergy"
                    factors = D*sqrt(m.^2+self.k^2);
                case "maxAbsG"
                    factors = self.auxiliary.normalizationFactors("maxAbsG");
                case "maxAbsF"
                    factors = D*m.*self.auxiliary.normalizationFactors("maxAbsG");
                otherwise
                    if ~isKey(self.normalizationRules,normalization)
                        error("IMZeroDepthWaveModesBasis:UnknownNormalization","Choose an available normalization or install a custom rule.");
                    end
                    rule = self.normalizationRules{normalization};
                    factors = arrayfun(@(j) rule(self,j),1:self.nColumns);
            end
            if any(~isfinite(factors) | factors <= 0)
                error("IMZeroDepthWaveModesBasis:NullNormalization","Normalization factors must be positive and finite.");
            end
        end

        function values = evaluate(self,z,options)
            % Evaluate manuscript profiles and classical interval derivatives.
            %
            % G supports orders zero through two; F supports zero and one.
            % Endpoint derivatives use the inside trace on each interval.
            %
            % - Topic: Evaluate zero-depth waves
            % - Declaration: values = evaluate(self,z,options)
            % - Parameter z: physical coordinate column within zDomain
            % - Parameter options.variable: required channel F or G
            % - Parameter options.derivativeOrder: supported nonnegative vertical derivative order, zero by default
            % - Parameter options.normalization: normalization rule for this call
            % - Returns values: sample-by-mode profiles or derivatives
            arguments (Input)
                self IMZeroDepthWaveModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.variable (1,1) string = ""
                options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 0
                options.normalization (1,1) string = self.normalization
            end
            arguments (Output)
                values double
            end
            variable = options.variable;
            derivativeOrder = options.derivativeOrder;
            if variable == ""
                error("IMZeroDepthWaveModesBasis:MissingVariable","Specify the channel with variable= followed by a supported channel name.");
            end
            if ~ismember(variable,["F","G"])
                error("IMZeroDepthWaveModesBasis:UnsupportedVariable","Unsupported channel %s. Inspect capabilities for supported channels and derivatives.",variable);
            end
            if variable == "G" && derivativeOrder <= 2
                values = self.auxiliary.evaluate(z,variable="G",derivativeOrder=derivativeOrder,normalization="unity");
            elseif variable == "F" && derivativeOrder <= 1
                values = self.background.D*self.auxiliary.evaluate(z,variable="G",derivativeOrder=derivativeOrder+1,normalization="unity");
            else
                error("IMZeroDepthWaveModesBasis:UnsupportedDerivative","G supports orders zero through two; F supports orders zero and one.");
            end
            values = values./self.normalizationFactors(options.normalization);
        end

        function values = F(self,z,options)
            % Evaluate F=D Gz.
            % - Topic: Evaluate zero-depth waves
            arguments
                self IMZeroDepthWaveModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.normalization (1,1) string = self.normalization
            end
            values = self.evaluate(z,variable="F",derivativeOrder=0,normalization=options.normalization);
        end

        function values = Fz(self,z,options)
            % Evaluate the classical interval derivative Fz=D Gzz.
            % - Topic: Evaluate zero-depth waves
            arguments
                self IMZeroDepthWaveModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.normalization (1,1) string = self.normalization
            end
            values = self.evaluate(z,variable="F",derivativeOrder=1,normalization=options.normalization);
        end

        function values = G(self,z,options)
            % Evaluate DD displacement profiles.
            % - Topic: Evaluate zero-depth waves
            arguments
                self IMZeroDepthWaveModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.normalization (1,1) string = self.normalization
            end
            values = self.evaluate(z,variable="G",derivativeOrder=0,normalization=options.normalization);
        end

        function values = Gz(self,z,options)
            % Evaluate first derivatives using interval-interior endpoint traces.
            % - Topic: Evaluate zero-depth waves
            arguments
                self IMZeroDepthWaveModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.normalization (1,1) string = self.normalization
            end
            values = self.evaluate(z,variable="G",derivativeOrder=1,normalization=options.normalization);
        end

        function values = Gzz(self,z,options)
            % Evaluate classical interval-interior second derivatives.
            % - Topic: Evaluate zero-depth waves
            arguments
                self IMZeroDepthWaveModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
                options.normalization (1,1) string = self.normalization
            end
            values = self.evaluate(z,variable="G",derivativeOrder=2,normalization=options.normalization);
        end

        function gram = gramMatrix(self,options)
            % Return exact depth-mean F, G or wave-energy pairings.
            % - Topic: Inspect zero-depth waves
            arguments
                self IMZeroDepthWaveModesBasis
                options.variable (1,1) string {mustBeMember(options.variable,["F","G","waveEnergy"])} = "waveEnergy"
                options.normalization (1,1) string = self.normalization
            end
            switch options.variable
                case "F", squares = (self.background.D*self.auxiliary.m_n).^2;
                case "G", squares = ones(1,self.nColumns);
                case "waveEnergy", squares = self.background.D^2*(self.auxiliary.m_n.^2+self.k^2);
            end
            factors = self.normalizationFactors(options.normalization);
            gram = diag(squares./factors.^2);
        end
    end
end

classdef IMNeutralLabelModesBasis
    % Evaluate independent zero-energy labels on exact neutral intervals.
    %
    % The caller specifies maximal intervals where N2 is exactly zero; no
    % inference from sampled density or stratification is performed. With
    % D = diff(zDomain), dimensionless profiles satisfy
    % $$G_{zz} = -m_{\mathrm n}^2 G,\quad D^{-1}\int G_i G_j\,dz = \delta_{ij}.$$
    % Each profile is zero outside its interval. Joins to stratified water
    % are Dirichlet. A physical surface or bottom is Neumann only when its
    % g0 or gd is zero; every nonzero value, including infinity, selects
    % Dirichlet. These modes have no equivalent depth or energy norm.
    %
    % A fully neutral Neumann-Neumann column has an APV-invisible constant
    % label. Construction requires constantLabelPolicy="exclude" to fix
    % that independent label to zero and retain only positive wavenumbers.
    % Amplitudes have units of length and require particle-label information
    % in addition to velocity and density. At finite amplitude, deleting a
    % component in this fixed linear basis need not preserve physical density.
    %
    % ```matlab
    % basis = IMNeutralLabelModesBasis([-200 -100],zDomain=[-1000 0],nModes=4);
    % G = basis.G(linspace(-1000,0,1001).');
    % projection = basis.projectionOnGrid(z,weights);
    % amplitudes = projection.project(APV);
    % ```
    %
    % - Topic: Create neutral labels
    % - Topic: Inspect neutral labels
    % - Topic: Evaluate neutral labels
    % - Topic: Project APV
    % - Declaration: classdef IMNeutralLabelModesBasis

    properties (SetAccess = private)
        % Increasing physical bottom and surface coordinates, in meters.
        % - Topic: Inspect neutral labels
        zDomain (1,2) double
        % Exact maximal neutral intervals, ordered from bottom to surface.
        % - Topic: Inspect neutral labels
        neutralIntervals (:,2) double
        % Positive retained mode count for each neutral interval.
        % - Topic: Inspect neutral labels
        nModes (1,:) double
        % Interval index for each retained column.
        % - Topic: Inspect neutral labels
        intervalIndex (1,:) double
        % Positive local mode number within each interval.
        % - Topic: Inspect neutral labels
        localModeNumber (1,:) double
        % Distinct interval/local-mode identities in retained-column order.
        % - Topic: Inspect neutral labels
        columnLabels (1,:) string
        % Positive vertical wavenumbers, in inverse meters.
        % - Topic: Inspect neutral labels
        m_n (1,:) double
        % Lower/upper boundary types for each interval: D or N.
        % - Topic: Inspect neutral labels
        boundaryTypes (:,2) string
        % Planetary rotation rate in radians per second.
        % - Topic: Inspect neutral labels
        rotationRate (1,1) double
        % Central latitude in degrees north.
        % - Topic: Inspect neutral labels
        latitude (1,1) double
        % Surface generalized-energy acceleration, in meters per second squared.
        % - Topic: Inspect neutral labels
        g0 (1,1) double
        % Bottom generalized-energy acceleration, in meters per second squared.
        % - Topic: Inspect neutral labels
        gd (1,1) double
        % Explicit policy for the fully neutral constant-label kernel.
        % - Topic: Inspect neutral labels
        constantLabelPolicy (1,1) string
        % Fixed depth normalization: integral G_i G_j / D equals delta_ij.
        % - Topic: Inspect neutral labels
        normalization (1,1) string = "depth"
    end

    properties (Dependent, SetAccess = private)
        % Derived Coriolis frequency, f0 = 2 rotationRate sin(latitude).
        % - Topic: Inspect neutral labels
        f0
    end

    methods
        function self = IMNeutralLabelModesBasis(neutralIntervals,options)
            % Create analytical neutral labels on explicitly specified intervals.
            %
            % Intervals must be nonempty, strictly separated, increasing and
            % contained in zDomain. Merge touching intervals before calling.
            % Counts retain modes 1:nModes per interval, in interval order.
            %
            % - Topic: Create neutral labels
            % - Parameter neutralIntervals: exact neutral intervals, one [lower upper] row each
            % - Parameter options.zDomain: physical [bottom surface] domain in meters
            % - Parameter options.nModes: positive scalar count or one count per interval
            % - Parameter options.rotationRate: nonnegative planetary rotation rate
            % - Parameter options.latitude: central latitude in degrees north
            % - Parameter options.g0: surface generalized-energy acceleration, zero by default
            % - Parameter options.gd: bottom generalized-energy acceleration, zero by default
            % - Parameter options.constantLabelPolicy: reject (default) or explicitly exclude the NN constant label
            % - Returns self: continuous neutral-label value basis
            arguments (Input)
                neutralIntervals (:,2) double {mustBeReal,mustBeFinite}
                options.zDomain (1,2) double {mustBeReal,mustBeFinite}
                options.nModes (1,:) double {mustBeInteger,mustBePositive,mustBeFinite} = 8
                options.rotationRate (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative} = 7.2921e-5
                options.latitude (1,1) double {mustBeReal,mustBeFinite,mustBeGreaterThanOrEqual(options.latitude,-90),mustBeLessThanOrEqual(options.latitude,90)} = 33
                options.g0 (1,1) double {mustBeReal} = 0
                options.gd (1,1) double {mustBeReal} = 0
                options.constantLabelPolicy (1,1) string {mustBeMember(options.constantLabelPolicy,["reject","exclude"])} = "reject"
            end
            arguments (Output)
                self (1,1) IMNeutralLabelModesBasis
            end
            if options.zDomain(1) >= options.zDomain(2)
                error("IMNeutralLabelModesBasis:InvalidDomain","zDomain must be strictly increasing.");
            end
            if isempty(neutralIntervals) || any(neutralIntervals(:,1) >= neutralIntervals(:,2)) || any(neutralIntervals(:,1) < options.zDomain(1)) || any(neutralIntervals(:,2) > options.zDomain(2))
                error("IMNeutralLabelModesBasis:InvalidIntervals","Supply nonempty increasing neutral intervals contained in zDomain.");
            end
            if any(neutralIntervals(2:end,1) <= neutralIntervals(1:end-1,2))
                error("IMNeutralLabelModesBasis:InvalidIntervals","Neutral intervals must be ordered, disjoint and strictly separated. Merge touching intervals.");
            end
            if isnan(options.g0) || isnan(options.gd)
                error("IMNeutralLabelModesBasis:InvalidEndpointAcceleration","g0 and gd must be real numbers or infinities, not NaN.");
            end
            counts = options.nModes;
            if isscalar(counts)
                counts = repmat(counts,1,size(neutralIntervals,1));
            elseif numel(counts) ~= size(neutralIntervals,1)
                error("IMNeutralLabelModesBasis:InvalidModeCounts","nModes must be scalar or contain one positive count per neutral interval.");
            end
            if isempty(counts)
                error("IMNeutralLabelModesBasis:InvalidModeCounts","Supply at least one positive mode count.");
            end
            self.zDomain = options.zDomain;
            self.neutralIntervals = neutralIntervals;
            self.nModes = counts;
            self.rotationRate = options.rotationRate;
            self.latitude = options.latitude;
            self.g0 = options.g0;
            self.gd = options.gd;
            self.constantLabelPolicy = options.constantLabelPolicy;
            self.boundaryTypes = repmat("D",size(neutralIntervals));
            self.boundaryTypes(neutralIntervals(:,1) == self.zDomain(1) & self.gd == 0,1) = "N";
            self.boundaryTypes(neutralIntervals(:,2) == self.zDomain(2) & self.g0 == 0,2) = "N";
            if any(all(self.boundaryTypes == "N",2)) && self.constantLabelPolicy ~= "exclude"
                error("IMNeutralLabelModesBasis:ConstantLabelKernel","A fully neutral NN column has an APV-invisible constant label. Set constantLabelPolicy=""exclude"" to fix that label to zero.");
            end
            self.intervalIndex = repelem(1:numel(counts),counts);
            self.localModeNumber = zeros(1,sum(counts));
            self.m_n = zeros(1,sum(counts));
            for iInterval = 1:numel(counts)
                columns = self.intervalIndex == iInterval;
                local = 1:counts(iInterval);
                self.localModeNumber(columns) = local;
                if sum(self.boundaryTypes(iInterval,:) == "N") == 1
                    local = local-0.5;
                end
                self.m_n(columns) = pi*local/diff(neutralIntervals(iInterval,:));
            end
            self.columnLabels = "interval" + string(self.intervalIndex) + ":mode" + string(self.localModeNumber);
        end

        function value = get.f0(self)
            value = 2*self.rotationRate*sind(self.latitude);
        end

        function values = G(self,z)
            % Evaluate dimensionless label profiles, zero outside their interval.
            % - Topic: Evaluate neutral labels
            % - Parameter z: physical coordinate column within zDomain
            % - Returns values: sample-by-mode label profiles
            arguments (Input)
                self IMNeutralLabelModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
            end
            arguments (Output)
                values (:,:) double
            end
            values = self.evaluateDerivative(z,0);
        end

        function values = Gz(self,z)
            % Evaluate interval-interior first derivatives, in inverse meters.
            %
            % At interval endpoints use the one-sided trace from inside the
            % neutral interval; derivatives may jump to zero outside it.
            % Quadrature must resolve intervals and use these traces consistently.
            %
            % - Topic: Evaluate neutral labels
            % - Parameter z: physical coordinate column within zDomain
            % - Returns values: sample-by-mode first derivatives
            arguments (Input)
                self IMNeutralLabelModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
            end
            arguments (Output)
                values (:,:) double
            end
            values = self.evaluateDerivative(z,1);
        end

        function values = Gzz(self,z)
            % Evaluate classical interval-interior second derivatives.
            %
            % Values use the inside trace at interval endpoints. This is not
            % the distributional derivative of the zero-extended Gz.
            %
            % - Topic: Evaluate neutral labels
            % - Parameter z: physical coordinate column within zDomain
            % - Returns values: sample-by-mode second derivatives in inverse meters squared
            arguments (Input)
                self IMNeutralLabelModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
            end
            arguments (Output)
                values (:,:) double
            end
            values = self.evaluateDerivative(z,2);
        end

        function values = APV(self,z)
            % Evaluate APV per unit displacement amplitude, q = -f0 Gz.
            % - Topic: Evaluate neutral labels
            % - Parameter z: physical coordinate column within zDomain
            % - Returns values: sample-by-mode APV profiles in inverse meter seconds
            arguments (Input)
                self IMNeutralLabelModesBasis
                z (:,1) double {mustBeReal,mustBeFinite}
            end
            arguments (Output)
                values (:,:) double
            end
            values = -self.f0*self.Gz(z);
        end

        function gram = gramMatrix(self,options)
            % Return exact continuous label or potential-enstrophy pairings.
            %
            % G uses integral G_i G_j / D = delta_ij. APV uses the
            % unscaled integral q_i q_j = f0^2 D m_n^2 delta_ij, with no
            % factor of one half. APV is degenerate when f0 is zero.
            %
            % - Topic: Inspect neutral labels
            % - Parameter options.variable: G (default) or APV
            % - Returns gram: exact continuous Gram matrix
            arguments (Input)
                self IMNeutralLabelModesBasis
                options.variable (1,1) string {mustBeMember(options.variable,["G","APV"])} = "G"
            end
            arguments (Output)
                gram (:,:) double
            end
            if options.variable == "G"
                gram = eye(numel(self.m_n));
            else
                gram = diag(self.f0^2*diff(self.zDomain)*self.m_n.^2);
            end
        end

        projection = projectionOnGrid(self,z,weights,options)
    end

    methods (Hidden)
        recipe = projectionRecipe(self,options)
    end

    methods (Access = private)
        function values = evaluateDerivative(self,z,order)
            if any(z < self.zDomain(1) | z > self.zDomain(2))
                error("IMNeutralLabelModesBasis:OutsideDomain","Evaluation coordinates must lie inside zDomain.");
            end
            values = zeros(numel(z),numel(self.m_n));
            for iInterval = 1:size(self.neutralIntervals,1)
                bounds = self.neutralIntervals(iInterval,:);
                rows = z >= bounds(1) & z <= bounds(2);
                if ~any(rows)
                    continue
                end
                columns = self.intervalIndex == iInterval;
                m = self.m_n(columns);
                phase = (z(rows)-bounds(1))*m;
                amplitude = sqrt(2*diff(self.zDomain)/diff(bounds));
                if self.boundaryTypes(iInterval,1) == "D"
                    shape = sin(phase + order*pi/2);
                else
                    shape = cos(phase + order*pi/2);
                end
                values(rows,columns) = amplitude*(m.^order).*shape;
            end
        end
    end
end

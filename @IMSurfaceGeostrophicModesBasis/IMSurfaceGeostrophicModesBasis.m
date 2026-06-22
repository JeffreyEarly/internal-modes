classdef IMSurfaceGeostrophicModesBasis
    % Store projected surface-geostrophic boundary modes.
    %
    % `IMSurfaceGeostrophicModesBasis` evaluates the projected zero-APV
    % boundary modes computed by a numerical solver from an
    % `IMSurfaceGeostrophicModes` problem. The solved mode amplitude is
    % `F(z)`. The diagnostic displacement variable is
    %
    % $$
    % G(z)=
    % -\frac{g}{N^2(z)}
    % \frac{\partial F}{\partial z}(z).
    % $$
    %
    % The equivalent boundary depths satisfy
    %
    % $$
    % h_0^a=2k^2\gamma_a,
    % $$
    %
    % where $$\gamma_a$$ are the eigenvalues of the boundary-energy
    % matrix used to project the raw endpoint modes.
    %
    % ```matlab
    % basisSet = solver.solveSurfaceGeostrophicModes(problem);
    % F = basisSet.F(z);
    % G = basisSet.G(z);
    % h = basisSet.h;
    % ```
    %
    % - Topic: Evaluate surface-geostrophic modes
    % - Topic: Inspect surface-geostrophic modes
    % - Declaration: classdef IMSurfaceGeostrophicModesBasis

    properties (SetAccess = private)
        % Surface-geostrophic problem descriptor.
        %
        % - Topic: Inspect surface-geostrophic modes
        problem

        % Numerical solver used to compute the native mode columns.
        %
        % - Topic: Inspect surface-geostrophic modes
        solver

        % Native projected mode columns before interpolation.
        %
        % - Topic: Inspect surface-geostrophic modes
        nativeModes

        % Mode-aligned horizontal wavenumbers.
        %
        % `k(j)` is the horizontal wavenumber for column `j` of `F(z)`
        % and `G(z)`.
        %
        % - Topic: Inspect surface-geostrophic modes
        k

        % Equivalent boundary depths.
        %
        % `h(j)` is the projected zero-APV boundary depth for column `j`
        % of `F(z)` and `G(z)`.
        %
        % - Topic: Inspect surface-geostrophic modes
        h

        % Projected mode labels within each wavenumber.
        %
        % - Topic: Inspect surface-geostrophic modes
        modeNumber

        % Raw endpoint-mode mixing coefficients.
        %
        % `mixingCoefficients` is a `2 x nModes` matrix. Row 1 is the
        % raw surface mode coefficient and row 2 is the raw bottom mode
        % coefficient; omitted endpoint modes have zero coefficients.
        %
        % - Topic: Inspect surface-geostrophic modes
        mixingCoefficients

        % Boundary-energy eigenvalues.
        %
        % `energyEigenvalues(j)` is the value $$\gamma_a$$ used in
        % $$h_j=2k_j^2\gamma_a$$.
        %
        % - Topic: Inspect surface-geostrophic modes
        energyEigenvalues

        % Surface buoyancy-anomaly weight.
        %
        % - Topic: Inspect surface-geostrophic modes
        g0

        % Bottom buoyancy-anomaly weight.
        %
        % - Topic: Inspect surface-geostrophic modes
        gd

        % Gravitational acceleration.
        %
        % - Topic: Inspect surface-geostrophic modes
        g

        % Surface anomaly convention.
        %
        % - Topic: Inspect surface-geostrophic modes
        surfaceAnomaly

        % Physical vertical domain.
        %
        % - Topic: Inspect surface-geostrophic modes
        zDomain

        % Buoyancy frequency squared function.
        %
        % - Topic: Inspect surface-geostrophic modes
        N2

        % Additional creation metadata.
        %
        % - Topic: Inspect surface-geostrophic modes
        metadata
    end

    methods
        function self = IMSurfaceGeostrophicModesBasis(options)
            % Create a solved surface-geostrophic basis.
            %
            % - Topic: Evaluate surface-geostrophic modes
            % - Declaration: basisSet = IMSurfaceGeostrophicModesBasis(options)
            % - Parameter options.problem: surface-geostrophic problem descriptor
            % - Parameter options.solver: configured numerical solver
            % - Parameter options.nativeModes: native projected mode columns
            % - Parameter options.k: mode-aligned horizontal wavenumbers
            % - Parameter options.h: equivalent boundary depths
            % - Parameter options.modeNumber: projected mode labels
            % - Parameter options.mixingCoefficients: raw endpoint-mode coefficients
            % - Parameter options.energyEigenvalues: boundary-energy eigenvalues
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: solved surface-geostrophic basis
            arguments
                options.problem IMSurfaceGeostrophicModes
                options.solver IMSolver
                options.nativeModes (:,:) double {mustBeReal, mustBeFinite}
                options.k double {mustBeReal, mustBeFinite, mustBePositive}
                options.h double {mustBeReal, mustBeFinite}
                options.modeNumber double {mustBeInteger, mustBePositive}
                options.mixingCoefficients (2,:) double {mustBeReal, mustBeFinite}
                options.energyEigenvalues double {mustBeReal, mustBeFinite}
                options.metadata struct = struct()
            end

            nModes = size(options.nativeModes,2);
            if numel(options.k) ~= nModes || numel(options.h) ~= nModes || numel(options.modeNumber) ~= nModes || numel(options.energyEigenvalues) ~= nModes
                error("IMSurfaceGeostrophicModesBasis:InvalidModeCount", "k, h, modeNumber, and energyEigenvalues must have one entry per native mode column.");
            end
            if size(options.mixingCoefficients,2) ~= nModes
                error("IMSurfaceGeostrophicModesBasis:InvalidModeCount", "mixingCoefficients must have one column per native mode column.");
            end

            self.problem = options.problem;
            self.solver = options.solver;
            self.nativeModes = options.nativeModes;
            self.k = reshape(options.k, 1, []);
            self.h = reshape(options.h, 1, []);
            self.modeNumber = reshape(options.modeNumber, 1, []);
            self.mixingCoefficients = options.mixingCoefficients;
            self.energyEigenvalues = reshape(options.energyEigenvalues, 1, []);
            self.g0 = options.problem.g0;
            self.gd = options.problem.gd;
            self.g = options.problem.g;
            self.surfaceAnomaly = options.problem.surfaceAnomaly;
            self.zDomain = options.problem.zDomain;
            self.N2 = options.problem.N2;
            self.metadata = options.metadata;
        end

        function values = F(self, z)
            % Evaluate projected SQG streamfunction modes.
            %
            % The returned array has one row per `z` value and one column
            % per projected boundary mode.
            %
            % - Topic: Evaluate surface-geostrophic modes
            % - Declaration: values = F(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: projected `F` values
            arguments
                self IMSurfaceGeostrophicModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            values = self.solver.evaluateNativeModes(self.nativeModes, z);
        end

        function values = G(self, z)
            % Evaluate diagnostic SQG displacement modes.
            %
            % `G(z)` is recovered from the projected `F` modes using
            % $$G=-gN^{-2}\partial F/\partial z$$.
            %
            % - Topic: Evaluate surface-geostrophic modes
            % - Declaration: values = G(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: diagnostic `G` values
            arguments
                self IMSurfaceGeostrophicModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            dFdz = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
            N2Values = self.N2(z(:));
            values = -(self.g./N2Values(:)).*dFdz;
        end

        function summarize(self)
            % Print a readable surface-geostrophic basis summary.
            %
            % - Topic: Inspect surface-geostrophic modes
            % - Declaration: summarize(basisSet)
            arguments
                self IMSurfaceGeostrophicModesBasis
            end

            fprintf("%s\n", class(self));
            fprintf("  zDomain: [%g, %g]\n", self.zDomain(1), self.zDomain(2));
            fprintf("  surfaceAnomaly: %s\n", self.surfaceAnomaly);
            fprintf("  g0: %s\n", self.formatEndpointWeight(self.g0));
            fprintf("  gd: %s\n", self.formatEndpointWeight(self.gd));
            fprintf("  nModes: %d\n", length(self.h));
            fprintf("  solver: %s\n", class(self.solver));
        end
    end

    methods (Access = private)
        function text = formatEndpointWeight(~, value)
            if isfinite(value)
                text = sprintf("%g", value);
            else
                text = "not represented";
            end
        end
    end
end

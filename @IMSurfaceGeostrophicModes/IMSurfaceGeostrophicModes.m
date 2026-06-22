classdef IMSurfaceGeostrophicModes
    % Describe projected surface-geostrophic boundary modes at fixed wavenumber.
    %
    % `IMSurfaceGeostrophicModes` stores the zero-APV boundary-mode
    % problem for horizontal wavenumber magnitude $$k$$. The raw endpoint
    % modes satisfy
    %
    % $$
    % \frac{\partial}{\partial z}
    % \left(
    % \frac{f_0^2}{N^2(z)}
    % \frac{\partial F}{\partial z}
    % \right)
    % -k^2F=0.
    % $$
    %
    % A finite `g0` includes the surface buoyancy-anomaly mode; a finite
    % `gd` includes the bottom buoyancy-anomaly mode. Omitted or infinite
    % endpoint weights mean that boundary anomaly is not represented. The
    % solver forms the raw endpoint modes, projects them with the
    % boundary-energy matrix, and returns an
    % `IMSurfaceGeostrophicModesBasis` with `F(z)`, `G(z)`, and `h`.
    %
    % The default surface anomaly includes the free-surface stretching term:
    %
    % $$
    % \eta_0[F]=
    % -\frac{f_0}{N^2(0)}
    % \frac{\partial F}{\partial z}(0)
    % -\frac{f_0}{g}F(0).
    % $$
    %
    % Use `surfaceAnomaly="noFreeSurface"` to omit the second term:
    %
    % $$
    % \eta_0[F]=
    % -\frac{f_0}{N^2(0)}
    % \frac{\partial F}{\partial z}(0).
    % $$
    %
    % The bottom anomaly is
    %
    % $$
    % \eta_d[F]=
    % -\frac{f_0}{N^2(z_b)}
    % \frac{\partial F}{\partial z}(z_b).
    % $$
    %
    % ```matlab
    % problem = IMSurfaceGeostrophicModes.atWavenumber(N2=N2,zDomain=[-4000 0],f0=1e-4,k=1e-4,g0=-0.035);
    % solver = IMSolverSpectral(nEVP=128);
    % basisSet = solver.solveSurfaceGeostrophicModes(problem);
    % F = basisSet.F(z);
    % G = basisSet.G(z);
    % h = basisSet.h;
    % ```
    %
    % - Topic: Create surface-geostrophic problems
    % - Topic: Summarize surface-geostrophic problems
    % - Topic: Inspect surface-geostrophic problems
    % - Declaration: classdef IMSurfaceGeostrophicModes

    properties (SetAccess = private)
        % Buoyancy frequency squared function.
        %
        % `N2` has signature `values = N2(z)` and defines the continuous
        % stratification used by numerical solvers.
        %
        % - Topic: Inspect surface-geostrophic problems
        N2

        % Physical vertical domain.
        %
        % - Topic: Inspect surface-geostrophic problems
        zDomain

        % Coriolis parameter.
        %
        % - Topic: Inspect surface-geostrophic problems
        f0

        % Gravitational acceleration.
        %
        % - Topic: Inspect surface-geostrophic problems
        g

        % Horizontal wavenumbers.
        %
        % - Topic: Inspect surface-geostrophic problems
        k

        % Surface buoyancy-anomaly weight.
        %
        % A finite nonzero `g0` includes the surface endpoint mode.
        % `Inf` means the surface anomaly is not represented.
        %
        % - Topic: Inspect surface-geostrophic problems
        g0

        % Bottom buoyancy-anomaly weight.
        %
        % A finite nonzero `gd` includes the bottom endpoint mode. `Inf`
        % means the bottom anomaly is not represented.
        %
        % - Topic: Inspect surface-geostrophic problems
        gd

        % Surface anomaly convention.
        %
        % `surfaceAnomaly` is `"freeSurface"` or `"noFreeSurface"`.
        % The default includes the free-surface stretching term
        % $$-f_0F(0)/g$$ in the surface anomaly.
        %
        % - Topic: Inspect surface-geostrophic problems
        surfaceAnomaly

        % Additional creation metadata.
        %
        % - Topic: Inspect surface-geostrophic problems
        metadata
    end

    methods
        function self = IMSurfaceGeostrophicModes(options)
            % Create a projected surface-geostrophic boundary-mode problem.
            %
            % - Topic: Create surface-geostrophic problems
            % - Declaration: problem = IMSurfaceGeostrophicModes(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.k: horizontal wavenumbers
            % - Parameter options.g0: surface buoyancy-anomaly weight
            % - Parameter options.gd: bottom buoyancy-anomaly weight
            % - Parameter options.surfaceAnomaly: surface anomaly convention
            % - Parameter options.metadata: additional metadata
            % - Returns problem: surface-geostrophic boundary-mode problem
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.f0 (1,1) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.k double {mustBeReal, mustBeFinite, mustBePositive}
                options.g0 (1,1) double {mustBeReal} = Inf
                options.gd (1,1) double {mustBeReal} = Inf
                options.surfaceAnomaly {mustBeTextScalar, mustBeMember(options.surfaceAnomaly, ["freeSurface", "noFreeSurface"])} = "freeSurface"
                options.metadata struct = struct()
            end

            if options.f0 == 0
                error("IMSurfaceGeostrophicModes:InvalidCoriolis", "f0 must be nonzero for surface-geostrophic modes.");
            end
            if any(isnan([options.g0 options.gd]))
                error("IMSurfaceGeostrophicModes:InvalidEndpointWeight", "g0 and gd must be finite nonzero values or Inf.");
            end
            if any([isfinite(options.g0) && options.g0 == 0, isfinite(options.gd) && options.gd == 0])
                error("IMSurfaceGeostrophicModes:InvalidEndpointWeight", "Finite g0 and gd values must be nonzero.");
            end
            if ~isfinite(options.g0) && ~isfinite(options.gd)
                error("IMSurfaceGeostrophicModes:NoBoundaryAnomaly", "At least one of g0 or gd must be finite to define surface-geostrophic modes.");
            end

            self.N2 = options.N2;
            self.zDomain = sort(options.zDomain);
            self.f0 = options.f0;
            self.g = options.g;
            self.k = reshape(options.k, 1, []);
            self.g0 = options.g0;
            self.gd = options.gd;
            self.surfaceAnomaly = string(options.surfaceAnomaly);
            self.metadata = options.metadata;
        end

        function summarize(self)
            % Print a readable problem summary.
            %
            % - Topic: Summarize surface-geostrophic problems
            % - Declaration: summarize(problem)
            arguments
                self IMSurfaceGeostrophicModes
            end

            fprintf("%s\n", class(self));
            fprintf("  zDomain: [%g, %g]\n", self.zDomain(1), self.zDomain(2));
            fprintf("  f0: %g\n", self.f0);
            fprintf("  g: %g\n", self.g);
            fprintf("  surfaceAnomaly: %s\n", self.surfaceAnomaly);
            fprintf("  g0: %s\n", self.formatEndpointWeight(self.g0));
            fprintf("  gd: %s\n", self.formatEndpointWeight(self.gd));
            fprintf("  nWavenumbers: %d\n", numel(self.k));
            fprintf("  modesPerWavenumber: %d\n", self.modesPerWavenumber());
        end

        function count = modesPerWavenumber(self)
            % Return the number of projected modes for each wavenumber.
            %
            % - Topic: Inspect surface-geostrophic problems
            % - Declaration: count = modesPerWavenumber(problem)
            % - Returns count: number of included endpoint anomaly modes
            arguments
                self IMSurfaceGeostrophicModes
            end

            count = double(isfinite(self.g0)) + double(isfinite(self.gd));
        end
    end

    methods (Static)
        function problem = atWavenumber(options)
            % Create projected surface-geostrophic modes at fixed wavenumber.
            %
            % - Topic: Create surface-geostrophic problems
            % - Declaration: problem = IMSurfaceGeostrophicModes.atWavenumber(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.k: horizontal wavenumbers
            % - Parameter options.g0: surface buoyancy-anomaly weight
            % - Parameter options.gd: bottom buoyancy-anomaly weight
            % - Parameter options.surfaceAnomaly: surface anomaly convention
            % - Parameter options.metadata: additional metadata
            % - Returns problem: projected SQG problem
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.f0 (1,1) double {mustBeReal, mustBeFinite}
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
                options.k double {mustBeReal, mustBeFinite, mustBePositive}
                options.g0 (1,1) double {mustBeReal} = Inf
                options.gd (1,1) double {mustBeReal} = Inf
                options.surfaceAnomaly {mustBeTextScalar, mustBeMember(options.surfaceAnomaly, ["freeSurface", "noFreeSurface"])} = "freeSurface"
                options.metadata struct = struct()
            end

            problem = IMSurfaceGeostrophicModes(N2=options.N2, zDomain=options.zDomain, f0=options.f0, g=options.g, k=options.k, g0=options.g0, gd=options.gd, surfaceAnomaly=options.surfaceAnomaly, metadata=options.metadata);
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

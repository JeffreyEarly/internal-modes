classdef IMSurfaceGeostrophicModes
    % Describe surface-geostrophic boundary modes at fixed wavenumber.
    %
    % `IMSurfaceGeostrophicModes` stores the physical boundary-value
    % problem for SQG boundary modes. For each horizontal wavenumber
    % $$k$$, the streamfunction $$\psi(z)$$ satisfies
    %
    % $$
    % \frac{\partial}{\partial z}
    % \left(
    % \frac{f_0^2}{N^2(z)}
    % \frac{\partial \psi}{\partial z}
    % \right)
    % -k^2\psi=0.
    % $$
    %
    % The active boundary is either the surface or the bottom. A surface
    % mode uses
    %
    % $$
    % f_0\frac{\partial\psi}{\partial z}(z_s)=1,\qquad
    % f_0\frac{\partial\psi}{\partial z}(z_b)=0,
    % $$
    %
    % while a bottom mode swaps the two endpoint conditions.
    %
    % ```matlab
    % problem = IMSurfaceGeostrophicModes.surfaceModesAtWavenumber(N2=N2,zDomain=[-4000 0],f0=1e-4,k=1e-4);
    % solver = IMSolverSpectral(nEVP=128);
    % basisSet = solver.solveSurfaceGeostrophicModes(problem);
    % psi = basisSet.psi(z);
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

        % Horizontal wavenumbers.
        %
        % - Topic: Inspect surface-geostrophic problems
        k

        % Active SQG boundary, `"surface"` or `"bottom"`.
        %
        % - Topic: Inspect surface-geostrophic problems
        boundary

        % Additional creation metadata.
        %
        % - Topic: Inspect surface-geostrophic problems
        metadata
    end

    methods
        function self = IMSurfaceGeostrophicModes(options)
            % Create a surface-geostrophic boundary-mode problem.
            %
            % - Topic: Create surface-geostrophic problems
            % - Declaration: problem = IMSurfaceGeostrophicModes(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.k: horizontal wavenumbers
            % - Parameter options.boundary: active boundary, `"surface"` or `"bottom"`
            % - Parameter options.metadata: additional metadata
            % - Returns problem: surface-geostrophic boundary-mode problem
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.f0 (1,1) double {mustBeReal, mustBeFinite}
                options.k double {mustBeReal, mustBeFinite, mustBePositive}
                options.boundary {mustBeTextScalar, mustBeMember(options.boundary, ["surface", "bottom"])} = "surface"
                options.metadata struct = struct()
            end

            if options.f0 == 0
                error("IMSurfaceGeostrophicModes:InvalidCoriolis", "f0 must be nonzero for surface-geostrophic modes.");
            end
            self.N2 = options.N2;
            self.zDomain = sort(options.zDomain);
            self.f0 = options.f0;
            self.k = reshape(options.k, 1, []);
            self.boundary = string(options.boundary);
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
            fprintf("  boundary: %s\n", self.boundary);
            fprintf("  zDomain: [%g, %g]\n", self.zDomain(1), self.zDomain(2));
            fprintf("  f0: %g\n", self.f0);
            fprintf("  nWavenumbers: %d\n", numel(self.k));
        end
    end

    methods (Static)
        function problem = surfaceModesAtWavenumber(options)
            % Create surface SQG modes at fixed wavenumber.
            %
            % - Topic: Create surface-geostrophic problems
            % - Declaration: problem = IMSurfaceGeostrophicModes.surfaceModesAtWavenumber(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.k: horizontal wavenumbers
            % - Parameter options.metadata: additional metadata
            % - Returns problem: surface SQG problem
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.f0 (1,1) double {mustBeReal, mustBeFinite}
                options.k double {mustBeReal, mustBeFinite, mustBePositive}
                options.metadata struct = struct()
            end

            problem = IMSurfaceGeostrophicModes(N2=options.N2, zDomain=options.zDomain, f0=options.f0, k=options.k, boundary="surface", metadata=options.metadata);
        end

        function problem = bottomModesAtWavenumber(options)
            % Create bottom SQG modes at fixed wavenumber.
            %
            % - Topic: Create surface-geostrophic problems
            % - Declaration: problem = IMSurfaceGeostrophicModes.bottomModesAtWavenumber(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.k: horizontal wavenumbers
            % - Parameter options.metadata: additional metadata
            % - Returns problem: bottom SQG problem
            arguments
                options.N2 function_handle
                options.zDomain (1,2) double {mustBeReal, mustBeFinite}
                options.f0 (1,1) double {mustBeReal, mustBeFinite}
                options.k double {mustBeReal, mustBeFinite, mustBePositive}
                options.metadata struct = struct()
            end

            problem = IMSurfaceGeostrophicModes(N2=options.N2, zDomain=options.zDomain, f0=options.f0, k=options.k, boundary="bottom", metadata=options.metadata);
        end
    end
end

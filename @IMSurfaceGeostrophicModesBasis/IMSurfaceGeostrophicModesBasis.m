classdef IMSurfaceGeostrophicModesBasis
    % Store solved surface-geostrophic boundary modes.
    %
    % `IMSurfaceGeostrophicModesBasis` evaluates the SQG streamfunction
    % modes computed by a numerical solver from an
    % `IMSurfaceGeostrophicModes` problem.
    %
    % ```matlab
    % basisSet = solver.solveSurfaceGeostrophicModes(problem);
    % psi = basisSet.psi(z);
    % psiz = basisSet.psiz(z);
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

        % Native mode columns before interpolation.
        %
        % - Topic: Inspect surface-geostrophic modes
        nativeModes

        % Horizontal wavenumbers.
        %
        % - Topic: Inspect surface-geostrophic modes
        k

        % Active SQG boundary, `"surface"` or `"bottom"`.
        %
        % - Topic: Inspect surface-geostrophic modes
        boundary

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
            % - Parameter options.nativeModes: native mode columns
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: solved surface-geostrophic basis
            arguments
                options.problem IMSurfaceGeostrophicModes
                options.solver IMSolver
                options.nativeModes (:,:) double {mustBeReal, mustBeFinite}
                options.metadata struct = struct()
            end

            if size(options.nativeModes,2) ~= numel(options.problem.k)
                error("IMSurfaceGeostrophicModesBasis:InvalidModeCount", "nativeModes must contain one column for each wavenumber.");
            end
            self.problem = options.problem;
            self.solver = options.solver;
            self.nativeModes = options.nativeModes;
            self.k = options.problem.k;
            self.boundary = options.problem.boundary;
            self.zDomain = options.problem.zDomain;
            self.N2 = options.problem.N2;
            self.metadata = options.metadata;
        end

        function values = psi(self, z)
            % Evaluate SQG streamfunction modes.
            %
            % The returned array has one row per `z` value and one column
            % per retained wavenumber in `k`.
            %
            % - Topic: Evaluate surface-geostrophic modes
            % - Declaration: values = psi(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: streamfunction values
            arguments
                self IMSurfaceGeostrophicModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            values = self.solver.evaluateNativeModes(self.nativeModes, z);
        end

        function values = psiz(self, z)
            % Evaluate vertical derivatives of SQG streamfunction modes.
            %
            % The returned array has one row per `z` value and one column
            % per retained wavenumber in `k`.
            %
            % - Topic: Evaluate surface-geostrophic modes
            % - Declaration: values = psiz(basisSet,z)
            % - Parameter z: physical coordinate
            % - Returns values: streamfunction derivative values
            arguments
                self IMSurfaceGeostrophicModesBasis
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            values = self.solver.evaluatePhysicalDerivative(self.nativeModes, z, 1);
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
            fprintf("  boundary: %s\n", self.boundary);
            fprintf("  zDomain: [%g, %g]\n", self.zDomain(1), self.zDomain(2));
            fprintf("  nWavenumbers: %d\n", numel(self.k));
            fprintf("  solver: %s\n", class(self.solver));
        end
    end
end

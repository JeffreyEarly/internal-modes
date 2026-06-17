classdef IMAnalyticalSQGBasis
    % Store exact SQG boundary modes from an analytical solution family.
    %
    % `IMAnalyticalSQGBasis` evaluates the boundary-trapped SQG streamfunction
    % modes available for a closed-form stratification family.
    %
    % ```matlab
    % solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-5000 0]);
    % sqg = solution.sqgModesAtWavenumber(1e-4,boundary="surface");
    % psi = sqg.psi(z);
    % ```
    %
    % - Topic: Evaluate SQG modes
    % - Topic: Inspect SQG modes
    % - Declaration: classdef IMAnalyticalSQGBasis

    properties (SetAccess = private)
        % Analytical solution family that created this SQG basis.
        %
        % - Topic: Inspect SQG modes
        solution

        % Horizontal wavenumbers.
        %
        % - Topic: Inspect SQG modes
        k

        % Active SQG boundary, `"surface"` or `"bottom"`.
        %
        % - Topic: Inspect SQG modes
        boundary

        % Physical vertical domain.
        %
        % - Topic: Inspect SQG modes
        zDomain

        % Buoyancy frequency squared function.
        %
        % - Topic: Inspect SQG modes
        N2

        % Additional creation metadata.
        %
        % - Topic: Inspect SQG modes
        metadata
    end

    properties (Access = private)
        psiFunction
    end

    methods
        function self = IMAnalyticalSQGBasis(options)
            % Create an exact SQG basis.
            %
            % - Topic: Evaluate SQG modes
            % - Declaration: sqg = IMAnalyticalSQGBasis(options)
            % - Parameter options.solution: analytical solution family
            % - Parameter options.k: horizontal wavenumbers
            % - Parameter options.boundary: active boundary
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.psiFunction: exact streamfunction evaluator
            % - Parameter options.metadata: creation metadata
            % - Returns sqg: exact SQG basis
            arguments
                options.solution IMAnalyticalSolution
                options.k double {mustBeReal, mustBeFinite, mustBePositive}
                options.boundary {mustBeTextScalar, mustBeMember(options.boundary, ["surface", "bottom"])}
                options.N2 (1,1) function_handle
                options.psiFunction (1,1) function_handle
                options.metadata struct = struct()
            end

            self.solution = options.solution;
            self.k = reshape(options.k,1,[]);
            self.boundary = string(options.boundary);
            self.zDomain = options.solution.zDomain;
            self.N2 = options.N2;
            self.psiFunction = options.psiFunction;
            self.metadata = options.metadata;
        end

        function values = psi(self, z)
            % Evaluate exact SQG streamfunction modes.
            %
            % The returned array has one row per `z` value and one column per
            % retained wavenumber in `k`.
            %
            % - Topic: Evaluate SQG modes
            % - Declaration: values = psi(sqg,z)
            % - Parameter z: physical coordinate
            % - Returns values: SQG streamfunction values
            arguments
                self IMAnalyticalSQGBasis
                z (:,1) double {mustBeReal, mustBeFinite}
            end

            values = self.psiFunction(z(:));
        end

        function summarize(self)
            % Print a readable SQG basis summary.
            %
            % - Topic: Inspect SQG modes
            % - Declaration: summarize(sqg)
            arguments
                self IMAnalyticalSQGBasis
            end

            fprintf("%s\n", class(self));
            fprintf("  solution: %s\n", class(self.solution));
            fprintf("  boundary: %s\n", self.boundary);
            fprintf("  zDomain: [%g, %g]\n", self.zDomain(1), self.zDomain(2));
            fprintf("  nWavenumbers: %d\n", numel(self.k));
        end
    end
end

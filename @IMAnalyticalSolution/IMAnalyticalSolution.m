classdef IMAnalyticalSolution
    % Describe a stratification family with closed-form mode solutions.
    %
    % `IMAnalyticalSolution` is the base class for exact stratification
    % families. A concrete solution object owns the continuous stratification
    % parameters and advertises which internal-mode and SQG branches are
    % available.
    %
    % ```matlab
    % solution = IMConstantStratificationSolution(N0=5.2e-3,zDomain=[-5000 0]);
    % availability = solution.internalModeAvailability(evp);
    % ```
    %
    % - Topic: Create analytical solutions
    % - Topic: Inspect analytical solutions
    % - Topic: Compute internal modes
    % - Topic: Compute SQG modes
    % - Declaration: classdef IMAnalyticalSolution

    properties (SetAccess = protected)
        % Physical vertical domain.
        %
        % - Topic: Inspect analytical solutions
        zDomain

        % Coriolis parameter in radians per second.
        %
        % - Topic: Inspect analytical solutions
        f0

        % Gravitational acceleration in meters per second squared.
        %
        % - Topic: Inspect analytical solutions
        g
    end

    methods
        function self = IMAnalyticalSolution(options)
            % Create an analytical solution family.
            %
            % - Topic: Create analytical solutions
            % - Declaration: solution = IMAnalyticalSolution(options)
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Returns solution: analytical solution family
            arguments
                options.zDomain (1,2) double {mustBeReal, mustBeFinite} = [-1 0]
                options.f0 (1,1) double {mustBeReal, mustBeFinite} = 0
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
            end

            self.zDomain = sort(options.zDomain);
            self.f0 = options.f0;
            self.g = options.g;
        end

        function summarize(self, evp)
            % Print a readable summary of the analytical solution.
            %
            % With an EVP input, the summary also reports internal-mode
            % availability for that problem.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: summarize(solution,evp)
            % - Parameter evp: optional internal-mode EVP
            arguments
                self IMAnalyticalSolution
                evp = []
            end

            fprintf("%s\n", class(self));
            fprintf("  zDomain: [%g, %g]\n", self.zDomain(1), self.zDomain(2));
            fprintf("  f0: %g\n", self.f0);
            fprintf("  g: %g\n", self.g);
            if ~isempty(evp)
                availability = self.internalModeAvailability(evp);
                if availability.isAvailable
                    fprintf("  internal modes: available\n");
                else
                    fprintf("  internal modes: unavailable\n");
                end
                fprintf("  reason: %s\n", availability.reason);
            end
        end
    end

    methods
        function values = N2(self, z)
            % Evaluate the stratification profile.
            %
            % Concrete solution families override this method with their
            % closed-form $$N^2(z)$$.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: values = N2(solution,z)
            % - Parameter z: physical coordinate
            % - Returns values: buoyancy frequency squared
            arguments
                self IMAnalyticalSolution
                z double {mustBeReal, mustBeFinite}
            end

            error("IMAnalyticalSolution:AbstractMethod", "%s must implement N2.", class(self));
        end

        function availability = internalModeAvailability(self, evp)
            % Report whether exact internal modes are available.
            %
            % Concrete solution families return a struct with
            % `isAvailable`, `reason`, `solutionKind`, `stratification`,
            % `supportedVariables`, `supportedInnerProducts`, and
            % `supportedNormalizations`.
            %
            % - Topic: Compute internal modes
            % - Declaration: availability = internalModeAvailability(solution,evp)
            % - Parameter evp: internal-mode EVP
            % - Returns availability: availability report struct
            arguments
                self IMAnalyticalSolution
                evp = []
            end

            error("IMAnalyticalSolution:AbstractMethod", "%s must implement internalModeAvailability.", class(self));
        end

        function basisSet = internalModes(self, evp, options)
            % Create an exact internal-mode basis when available.
            %
            % Concrete solution families return an
            % `IMAnalyticalInternalModesBasis` for supported EVPs and throw a
            % class-specific unavailable or unsupported error otherwise.
            %
            % - Topic: Compute internal modes
            % - Declaration: basisSet = internalModes(solution,evp,options)
            % - Parameter evp: internal-mode EVP
            % - Parameter options.nModes: number of retained modes
            % - Parameter options.normalization: active normalization
            % - Parameter options.metadata: additional metadata
            % - Returns basisSet: exact analytical internal-mode basis
            arguments
                self IMAnalyticalSolution
                evp = []
                options.nModes (1,1) double {mustBeInteger, mustBePositive} = 64
                options.normalization = []
                options.metadata struct = struct()
            end

            error("IMAnalyticalSolution:AbstractMethod", "%s must implement internalModes.", class(self));
        end

        function availability = sqgAvailability(self, options)
            % Report whether exact SQG boundary modes are available.
            %
            % - Topic: Compute SQG modes
            % - Declaration: availability = sqgAvailability(solution,options)
            % - Parameter options.boundary: `"surface"` or `"bottom"`
            % - Returns availability: availability report struct
            arguments
                self IMAnalyticalSolution
                options.boundary {mustBeTextScalar, mustBeMember(options.boundary, ["surface", "bottom"])} = "surface"
            end

            error("IMAnalyticalSolution:AbstractMethod", "%s must implement sqgAvailability.", class(self));
        end

        function basisSet = sqgModesAtWavenumber(self, k, options)
            % Create exact SQG boundary modes at fixed wavenumber.
            %
            % Concrete solution families return an `IMAnalyticalSQGBasis`
            % when the requested boundary mode is available.
            %
            % - Topic: Compute SQG modes
            % - Declaration: sqg = sqgModesAtWavenumber(solution,k,options)
            % - Parameter k: horizontal wavenumbers
            % - Parameter options.boundary: `"surface"` or `"bottom"`
            % - Parameter options.metadata: additional metadata
            % - Returns sqg: exact SQG basis
            arguments
                self IMAnalyticalSolution
                k double {mustBeReal, mustBeFinite, mustBePositive}
                options.boundary {mustBeTextScalar, mustBeMember(options.boundary, ["surface", "bottom"])} = "surface"
                options.metadata struct = struct()
            end

            error("IMAnalyticalSolution:AbstractMethod", "%s must implement sqgModesAtWavenumber.", class(self));
        end
    end

    methods (Static, Access = protected)
        function availability = availabilityStruct(isAvailable, solutionKind, stratification, reason)
            availability = struct();
            availability.isAvailable = isAvailable;
            availability.solutionKind = string(solutionKind);
            availability.stratification = string(stratification);
            availability.reason = string(reason);
            availability.supportedVariables = ["F", "G"];
            availability.supportedInnerProducts = ["F", "G"];
            availability.supportedNormalizations = ["unity", "uMax", "wMax", "surfacePressure"];
        end
    end
end

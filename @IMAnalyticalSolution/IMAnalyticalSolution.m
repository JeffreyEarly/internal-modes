classdef (Abstract) IMAnalyticalSolution
    % Describe common state for a stratification family with exact solutions.
    %
    % `IMAnalyticalSolution` owns only physical state shared by analytical
    % stratification families. Concrete classes advertise exact solution
    % operations by implementing those construction methods directly.
    %
    % - Topic: Create analytical solutions
    % - Topic: Inspect analytical solutions
    % - Declaration: classdef (Abstract) IMAnalyticalSolution

    properties (SetAccess = protected)
        % Physical vertical domain.
        %
        % - Topic: Inspect analytical solutions
        zDomain

        % Gravitational acceleration in meters per second squared.
        %
        % - Topic: Inspect analytical solutions
        g
    end

    properties (SetAccess = private)
        % Planetary rotation rate in radians per second.
        %
        % - Topic: Inspect analytical solutions
        rotationRate

        % Latitude in degrees north.
        %
        % - Topic: Inspect analytical solutions
        latitude
    end

    properties (Dependent, SetAccess = private)
        % Coriolis parameter in radians per second.
        %
        % $$f_0=2\Omega\sin\varphi$$ is derived from `rotationRate` and
        % `latitude`.
        %
        % - Topic: Inspect analytical solutions
        f0
    end

    methods
        function self = IMAnalyticalSolution(options)
            % Initialize common analytical-stratification state.
            %
            % - Topic: Create analytical solutions
            % - Declaration: solution = IMAnalyticalSolution(options)
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.rotationRate: nonnegative planetary rotation rate in radians per second
            % - Parameter options.latitude: latitude in degrees north, between -90 and 90
            % - Parameter options.g: gravitational acceleration
            % - Returns solution: analytical solution family
            arguments
                options.zDomain (1,2) double {mustBeReal, mustBeFinite} = [-1 0]
                options.rotationRate (1,1) double {mustBeReal, mustBeFinite, mustBeNonnegative} = 7.2921e-5
                options.latitude (1,1) double {mustBeReal, mustBeFinite, mustBeGreaterThanOrEqual(options.latitude,-90), mustBeLessThanOrEqual(options.latitude,90)} = 33
                options.g (1,1) double {mustBeReal, mustBeFinite, mustBePositive} = 9.81
            end

            self.zDomain = sort(options.zDomain);
            self.rotationRate = options.rotationRate;
            self.latitude = options.latitude;
            self.g = options.g;
        end

        function value = get.f0(self)
            value = 2*self.rotationRate*sind(self.latitude);
        end

        function summarize(self)
            % Print common analytical-stratification state.
            %
            % Concrete families extend this summary with their parameters and
            % exact construction methods.
            %
            % - Topic: Inspect analytical solutions
            % - Declaration: summarize(solution)
            arguments
                self IMAnalyticalSolution
            end

            fprintf("%s\n",class(self));
            fprintf("  zDomain: [%g, %g]\n",self.zDomain(1),self.zDomain(2));
            fprintf("  rotationRate: %g rad/s\n",self.rotationRate);
            fprintf("  latitude: %g degrees north\n",self.latitude);
            fprintf("  f0: %g\n",self.f0);
            fprintf("  g: %g\n",self.g);
        end
    end

    methods (Abstract)
        % Evaluate buoyancy frequency squared $$N^2(z)$$.
        %
        % - Topic: Inspect analytical solutions
        % - Declaration: values = N2(solution,z)
        % - Parameter z: physical coordinate
        % - Returns values: buoyancy frequency squared
        values = N2(self,z)
    end
end

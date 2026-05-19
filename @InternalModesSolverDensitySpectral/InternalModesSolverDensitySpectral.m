classdef InternalModesSolverDensitySpectral < InternalModesSolverSpectral
    % Solve physical-coordinate EVPs in the density-like coordinate.
    %
    % The native coordinate satisfies $$dx/dz=N^2(z)$$. The inherited solver
    % applies this coordinate pullback automatically.
    %
    % ```matlab
    % solver = InternalModesSolverDensitySpectral(N2=@(z) 1e-5*exp(z/1000), zDomain=[-1000 0], nEVP=64);
    % ```
    %
    % - Topic: Create solvers
    % - Declaration: classdef InternalModesSolverDensitySpectral < InternalModesSolverSpectral

    methods
        function self = InternalModesSolverDensitySpectral(options)
            % Create a density-coordinate spectral solver.
            %
            % - Topic: Create solvers
            % - Declaration: solver = InternalModesSolverDensitySpectral(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.nEVP: number of EVP coefficients
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Returns solver: initialized density spectral solver
            arguments
                options.N2 function_handle = @(z) 1e-5*ones(size(z))
                options.zDomain (1,2) double = [-1 0]
                options.nEVP (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.nEVP, 4)} = 64
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
            end

            self@InternalModesSolverSpectral(N2=options.N2, zDomain=options.zDomain, ...
                nEVP=options.nEVP, f0=options.f0, g=options.g, coordinateKind="density");
        end
    end
end

classdef IMSolverWKBSpectral < IMSolverSpectral
    % Solve physical-coordinate EVPs in the WKB stretched coordinate.
    %
    % The native coordinate satisfies $$dx/dz=N(z)$$. The inherited solver
    % applies this coordinate pullback automatically.
    %
    % ```matlab
    % solver = IMSolverWKBSpectral(N2=@(z) 1e-5*exp(z/1000), zDomain=[-1000 0], nEVP=64);
    % ```
    %
    % - Topic: Create solvers
    % - Declaration: classdef IMSolverWKBSpectral < IMSolverSpectral

    methods
        function self = IMSolverWKBSpectral(options)
            % Create a WKB-coordinate spectral solver.
            %
            % - Topic: Create solvers
            % - Declaration: solver = IMSolverWKBSpectral(options)
            % - Parameter options.N2: buoyancy frequency squared function
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.nEVP: number of EVP coefficients
            % - Returns solver: initialized WKB spectral solver
            arguments
                options.N2 function_handle = @(z) 1e-5*ones(size(z))
                options.zDomain (1,2) double = [-1 0]
                options.nEVP (1,1) double {mustBeInteger, mustBeGreaterThanOrEqual(options.nEVP, 4)} = 64
            end

            self@IMSolverSpectral(N2=options.N2, zDomain=options.zDomain, nEVP=options.nEVP, coordinateKind="wkb");
        end
    end
end

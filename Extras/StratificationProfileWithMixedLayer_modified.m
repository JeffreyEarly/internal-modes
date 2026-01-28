function [N2,zDomain] = StratificationProfileWithMixedLayer_modified(options)
%StratificationProfileWithMixedLayer
%   mixedLayerDepth is the depth of the mixed layer 
% (should be a negative number, i.e., below the surface) (m)
%   dRhoTransitionLayer is the density jump across the transition layer, 
% 3 would be the high-end (kg/m^3)
%   width of the transition layer (m)
%   buoyancy frequency at the surface (m/s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
arguments
    options.mixedLayerDepth (1,1) double {mustBeNegative} = -100
    options.dRhoTransitionLayer (1,1) double {mustBeNonnegative} = 0.5
    options.dZTransitionLayer (1,1) double {mustBePositive} = 15
    options.N0 (1,1) double {mustBePositive} = 5e-4
end

z_p = options.mixedLayerDepth;
dRhoML = options.dRhoTransitionLayer;
N0 = options.N0; 
delta_p = options.dZTransitionLayer; 

Lgm = 1300; % e-fold scale of stratification below the mixed layer
Nb = (6.5e-3)*exp(-4000/Lgm); % buoyancy freqency at the bottom. 
% This is tuned so that the equivalent depth is approx 80 cm.
D = 4000; % depth of the ocean

N2_ml = dRhoML*(9.81)/(2*1025*delta_p);

% Stratification profile---deep exponential, then rapid transition to mixed
% layer at the surface.
N2 = @(z) N0*N0*exp( ((z+D)/Lgm+log(Nb)-log(N0)) .* (1-tanh((z-z_p)/delta_p)) ) ...
    + N2_ml*sech( (z-z_p)/delta_p ).^2;

zDomain = [-D, 0];
end


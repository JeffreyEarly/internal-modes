function [collection,diagnostics] = solveWaveModesAtWavenumbers(self,varargin)
% Solve continuous wave bases for an ordered collection of wavenumbers.
%
% Each distinct requested wavenumber is solved exactly once. Coordinates,
% stratification samples, derivative matrices and the shared pencil terms
% are prepared once for this call. The wave pencil is
% $$A(\kappa)=A_0+\kappa^2 A_2,\qquad B(\kappa)=B_0,$$
% with zero endpoint rows in $$A_2$$ because boundary equations replace the
% interior operator there. No approximate sharing or global cache is used.
% Shared preparation supports the exact built-in IMSolverSpectral and
% IMSolverFiniteDifference classes. Custom solvers must use independent
% solveEVP calls until they expose a validated preparation-reuse capability.
%
% Positive wavenumbers request `nModes` columns, either a uniform scalar or
% one count per positive entry of `kappa`, preserving its order and repeats.
% Repeated wavenumbers must request identical counts. A zero-wavenumber request
% requires an explicit independent `nInertialModes` count. Both counts are
% exact requests: an insufficient returned family raises an error. Stored
% bases retain the same scientific labels, normalization and boundary
% conditions as independent `solveEVP` calls. Frequency signs and physical
% wave polarizations belong to the caller, not this vertical-basis solve.
%
% ```matlab
% [bases,costs] = solver.solveWaveModesAtWavenumbers([2e-4 0 1e-4 2e-4],N2=N2,zDomain=[-1000 0],nModes=8,nInertialModes=5);
% G = bases.evaluate(z,variable="G",pages=[1 3]);
% ```
%
% - Topic: Solve EVPs
% - Declaration: [collection,diagnostics] = solveWaveModesAtWavenumbers(solver,kappa,options)
% - Parameter kappa: nonempty row of requested nonnegative wavenumbers, in radians per meter
% - Parameter options.N2: buoyancy frequency squared evaluator
% - Parameter options.zDomain: increasing physical vertical bounds
% - Parameter options.rotationRate: planetary rotation rate in radians per second; default Earth rotation
% - Parameter options.latitude: central latitude in degrees north; default 33
% - Parameter options.g: gravitational acceleration
% - Parameter options.surfaceBoundary: surface condition of the G EVP
% - Parameter options.bottomBoundary: bottom condition of the G EVP
% - Parameter options.nModes: required scalar or row of nonnegative counts aligned with kappa(kappa>0); may be omitted only for all-zero kappa
% - Parameter options.nInertialModes: independent count, required when zero is requested
% - Returns collection: continuous bases and exact requested-page mapping
% - Returns diagnostics: construction costs and solve provenance, not accuracy qualification
arguments
    self (1,1) IMSolver
end
arguments (Repeating)
    varargin
end
args = varargin;
if ~isempty(args) && isa(args{1},"IMBackground")
    background = args{1}; args(1) = [];
    if isempty(args), error("IMSolver:MissingWavenumbers","Supply background and kappa."); end
    kappa = args{1}; args(1) = [];
    if mod(numel(args),2) ~= 0, error("IMSolver:InvalidOptions","Use name-value options after kappa."); end
    names = string(args(1:2:end));
    physicalArgs = args;
    for name = ["nModes","nInertialModes","normalization"]
        index = find(string(physicalArgs(1:2:end)) == name);
        physicalArgs(reshape([2*index-1;2*index],1,[])) = [];
    end
    problem = IMInternalModes.waveModesAtWavenumber(background,physicalArgs{:},k=0);
    defaults = struct(N2=problem.N2,zDomain=problem.zDomain,rotationRate=problem.rotationRate,latitude=problem.latitude,g=problem.g,surfaceBoundary=problem.surfaceBoundary,bottomBoundary=problem.bottomBoundary);
    for name = reshape(string(fieldnames(defaults)),1,[])
        if ~ismember(name,names), args(end+(1:2)) = {name,defaults.(name)}; end
    end
    args = [{kappa} args {"background",background}];
end
[collection,diagnostics] = self.solveWaveModesAtWavenumbersFromOptions(args{:});
end

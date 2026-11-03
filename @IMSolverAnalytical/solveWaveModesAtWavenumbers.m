function [collection,diagnostics] = solveWaveModesAtWavenumbers(self,background,kappa,options)
% Solve each distinct exact analytical wave page once.
%
% Counts include every retained column. Positive pages use nModes, scalar
% or aligned with kappa(kappa>0). Zero-wavenumber pages require the explicit
% independent nInertialModes count. Repeated pages must request equal counts.
% No approximate sharing or global cache is used.
%
% - Topic: Solve analytical problems
% - Parameter background: authoritative source-specific physical column
% - Parameter kappa: nonempty requested nonnegative horizontal wavenumbers
% - Parameter options.nModes: required exact total count per positive page; may be omitted only for all-zero kappa
% - Parameter options.nInertialModes: independent exact total count for zero pages
% - Parameter options.normalization: optional active profile normalization
% - Returns collection: capable analytical bases and requested-page mapping
% - Returns diagnostics: exact repeated-page reuse and retained counts
arguments
    self (1,1) IMSolverAnalytical
    background (1,1) IMBackground
    kappa (1,:) double {mustBeReal,mustBeFinite,mustBeNonnegative}
    options.nModes (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = []
    options.nInertialModes double {mustBeInteger,mustBeNonnegative,mustBeFinite} = []
    options.normalization = []
end
if isempty(kappa)
    error("IMSolver:EmptyWavenumberCollection","Supply at least one requested wavenumber.");
end
positive = kappa > 0;
if any(positive) && isempty(options.nModes)
    error("IMSolver:MissingModeCount","Positive-wavenumber pages require an explicit nModes count. Supply a scalar or one count per positive entry of kappa; zero is valid.");
end
if ~isempty(options.nInertialModes) && ~isscalar(options.nInertialModes)
    error("IMSolver:InvalidInertialModeCount","nInertialModes must be a scalar count.");
end
if any(kappa == 0) && isempty(options.nInertialModes)
    error("IMSolver:MissingInertialModeCount","A zero-wavenumber page requires an explicit nInertialModes count independent of nModes.");
end
if isscalar(options.nModes)
    counts = repmat(options.nModes,size(kappa));
elseif numel(options.nModes) == nnz(positive)
    counts = zeros(size(kappa));
    counts(positive) = options.nModes;
else
    error("IMSolver:InvalidWaveModeCounts","nModes must be scalar or contain one count per positive entry of kappa.");
end
if any(~positive), counts(~positive) = options.nInertialModes; end
[distinct,first,mapping] = unique(kappa,"stable");
bases = cell(1,numel(distinct));
for j = 1:numel(distinct)
    if numel(unique(counts(kappa == distinct(j)))) ~= 1
        error("IMSolver:ConflictingWaveModeCounts","Repeated kappa values must request identical mode counts.");
    end
    problem = IMInternalModes.waveModesAtWavenumber(background,k=distinct(j));
    bases{j} = self.solveEVP(problem,nModes=counts(first(j)),normalization=options.normalization);
end
collection = IMBasisCollection(bases,kappa=kappa,basisIndex=reshape(mapping,1,[]));
diagnostics = struct(requestedKappa=kappa,solvedKappa=distinct,basisIndex=reshape(mapping,1,[]),requestedPageCount=numel(kappa),solveCount=numel(distinct),reusedPageCount=numel(kappa)-numel(distinct),reuseKind="exactEqualKappa",solverClass="IMSolverAnalytical",nModes=options.nModes,nInertialModes=options.nInertialModes,requestedCounts=counts(first));
end

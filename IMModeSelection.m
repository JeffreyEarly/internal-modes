classdef (Hidden, Sealed) IMModeSelection
    % Select scientific ladders separately from finite exceptional sectors.
    methods (Static)
        function selection = selectCounted(evp,lambda,nModes,branch,options)
            % Select an exact total inventory, or explicit counted branches.
            arguments (Input)
                evp IMEigenvalueProblem
                lambda (:,1) double {mustBeReal,mustBeFinite}
                nModes (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite}
                branch (1,:) string = ""
                options.solver = []
                options.A double = []
                options.diagnostics (1,1) struct = struct()
                options.zeroIndices (1,:) double = []
            end
            explicit = ~(isscalar(branch) && ismember(branch,["","all"]));
            if (~explicit && ~isscalar(nModes)) || (explicit && numel(branch) ~= numel(nModes)) || numel(unique(branch)) ~= numel(branch)
                error("IMSolver:InvalidModeRequest","Use one total nModes count, or one count for each distinct explicitly requested branch.");
            end
            diagnostics = options.diagnostics;
            wave = isa(evp,"IMInternalModes") && ismember(string(evp.name),["waveModesAtWavenumber","waveModesAtFrequency"]);
            if wave
                positive = find(lambda > 0);
                [~,order] = sort(lambda(positive)); positive = positive(order);
                negative = find(lambda < 0);
                [~,order] = sort(lambda(negative),"descend"); negative = negative(order);
                surface = zeros(0,1);
                if evp.surfaceBoundary.isEigenvalueDependent() && ~isempty(positive)
                    surface = positive(1); positive(1) = [];
                end
                if evp.name == "waveModesAtWavenumber"
                    k = evp.parameters.k;
                    if k > 0
                        negative = negative(evp.f0^2+evp.g*k^2./lambda(negative) > 0);
                    end
                    index = [surface;positive;negative];
                    branches = [repmat("surface",numel(surface),1);repmat("internal",numel(positive),1);repmat("subinertial",numel(negative),1)];
                    numbers = [ones(numel(surface),1);(1:numel(positive))';(1:numel(negative))'];
                else
                    k2 = (evp.parameters.omega^2-evp.f0^2)*lambda/evp.g;
                    propagating = find(k2 > 0 & ~ismember((1:numel(lambda))',surface));
                    evanescent = find(k2 < 0 & ~ismember((1:numel(lambda))',surface));
                    [~,order] = sort(abs(lambda(propagating))); propagating = propagating(order);
                    [~,order] = sort(abs(lambda(evanescent))); evanescent = evanescent(order);
                    index = [surface;propagating;evanescent];
                    branches = [repmat("surface",numel(surface),1);repmat("propagating",numel(propagating),1);repmat("evanescent",numel(evanescent),1)];
                    numbers = [ones(numel(surface),1);(1:numel(propagating))';(1:numel(evanescent))'];
                end
                labels = 1:numel(index);
                labels(lambda(index)<0) = -1:-1:-nnz(lambda(index)<0);
            else
                if ~isempty(options.solver) && ~isempty(lambda)
                    base = evp.selectModes(lambda,max(1,numel(lambda)),options.solver,options.A,diagnostics=diagnostics);
                    index = base.sortIndex;
                    labels = base.modeNumber;
                    diagnostics = base.modeSelectionDiagnostics;
                else
                    negative = find(lambda < 0); [~,order] = sort(lambda(negative)); negative = negative(order);
                    zero = options.zeroIndices(:);
                    positive = find(lambda > 0 & ~ismember((1:numel(lambda))',zero));
                    [~,order] = sort(lambda(positive)); positive = positive(order);
                    index = [negative;zero;positive];
                    labels = [-1:-1:-numel(negative),zeros(1,numel(zero)),1:numel(positive)];
                end
                branches = repmat("positive",numel(index),1);
                branches(labels < 0) = "negative";
                branches(labels == 0) = "zeroEigenvalue";
                if isa(evp,"IMInternalModes") && evp.modeFamily == "meanDensityAnomaly"
                    branches(labels == 0) = "constant";
                end
                if isa(evp,"IMInternalModes") && evp.name == "hydrostaticFModes"
                    branches(labels == 0) = "barotropic";
                end
                numbers = abs(labels(:));
            end
            inventory = table(branches+":"+string(numbers),branches,numbers,lambda(index),VariableNames=["columnLabel","branch","branchIndex","eigenvalue"]);
            if explicit
                positions = zeros(0,1);
                for j = 1:numel(branch)
                    available = find(inventory.branch == branch(j));
                    if numel(available) < nModes(j)
                        error("IMSolver:InsufficientModes","Requested %d %s columns, but only %d are available. Increase resolution or revise the request.",nModes(j),branch(j),numel(available));
                    end
                    positions = [positions;available(1:nModes(j))]; %#ok<AGROW>
                end
            else
                if numel(index) < nModes
                    error("IMSolver:InsufficientModes","Requested %d total columns, but only %d are available. Increase resolution or revise the request.",nModes,numel(index));
                end
                positions = (1:nModes)';
            end
            selection = struct("sortIndex",index(positions),"modeNumber",reshape(labels(positions),1,[]),"modeSelectionDiagnostics",diagnostics,"modeInfo",inventory(positions,:));
        end

    end
end

classdef (Hidden, Sealed) IMModeSelection
    % Select scientific ladders separately from finite exceptional sectors.
    methods (Static)
        function selection = select(evp,lambda,request,options)
            arguments (Input)
                evp IMInternalModes
                lambda (:,1) double {mustBeReal,mustBeFinite}
                request (1,1) struct
                options.solver = []
                options.A double = []
                options.diagnostics (1,1) struct = struct()
                options.zeroIndices (1,:) double = []
            end
            diagnostics = options.diagnostics;
            kind = request.kind;
            if ismember(kind,["wave","inertial","frequency"])
                positive = find(lambda > 0);
                [~,order] = sort(lambda(positive));
                positive = positive(order);
                negative = find(lambda < 0);
                [~,order] = sort(lambda(negative),"descend");
                negative = negative(order);
                if kind == "frequency"
                    surface = zeros(0,1);
                    if request.hasSurface
                        if isempty(positive), error("IMModes:MissingSpecialMode","The fixed-frequency surface mode was not resolved."); end
                        surface = positive(1);
                    end
                    numerator = evp.parameters.omega^2-evp.f0^2;
                    k2 = numerator*lambda/evp.g;
                    propagating = find(k2 > 0 & ~ismember((1:numel(lambda))',surface));
                    evanescent = find(k2 < 0 & ~ismember((1:numel(lambda))',surface));
                    [~,order] = sort(abs(lambda(propagating))); propagating = propagating(order);
                    [~,order] = sort(abs(lambda(evanescent))); evanescent = evanescent(order);
                    p = IMModeSelection.take(propagating,request.nPropagatingModes,"propagating");
                    e = IMModeSelection.take(evanescent,request.nEvanescentModes,"evanescent");
                    index = [surface(:);p(:);e(:)];
                    branch = [repmat("surface",numel(surface),1);repmat("propagating",numel(p),1);repmat("evanescent",numel(e),1)];
                    number = [ones(numel(surface),1);(1:numel(p))';(1:numel(e))'];
                else
                    surface = zeros(0,1);
                    if request.hasSurface
                        if isempty(positive)
                            error("IMModes:MissingSpecialMode","The free-surface mode was not resolved. Increase the solver resolution.");
                        end
                        surface = positive(1); positive(1) = [];
                    end
                    if kind == "wave"
                        k = evp.parameters.k;
                        frequencySquared = evp.f0^2+evp.g*k^2./lambda(negative);
                        negative = negative(frequencySquared > 0);
                    end
                    p = IMModeSelection.take(positive,request.nPositiveModes,"internal/positive-depth");
                    n = IMModeSelection.take(negative,request.nNegativeModes,"subinertial/negative-depth");
                    index = [surface(:);p(:);n(:)];
                    names = ["surface","internal","subinertial"];
                    if kind == "inertial", names = ["external","positiveDepth","negativeDepth"]; end
                    branch = [repmat(names(1),numel(surface),1);repmat(names(2),numel(p),1);repmat(names(3),numel(n),1)];
                    number = [ones(numel(surface),1);(1:numel(p))';(1:numel(n))'];
                end
                labels = 1:numel(index);
                labels(lambda(index)<0) = -1:-1:-nnz(lambda(index)<0);
            else
                if ~isempty(options.solver)
                    base = evp.selectModes(lambda,max(1,numel(lambda)),options.solver,options.A,diagnostics=diagnostics);
                    candidates = base.sortIndex;
                    diagnostics = base.modeSelectionDiagnostics;
                    zero = candidates(base.modeNumber == 0);
                else
                    candidates = (1:numel(lambda))';
                    zero = options.zeroIndices(:);
                end
                negative = candidates(lambda(candidates)<0 & ~ismember(candidates,zero));
                [~,order] = sort(lambda(negative)); negative = negative(order);
                positive = candidates(lambda(candidates)>0 & ~ismember(candidates,zero));
                [~,order] = sort(lambda(positive)); positive = positive(order);
                surface = zeros(0,1);
                if kind == "hydrostaticG" && request.hasSurface
                    if isempty(positive), error("IMModes:MissingSpecialMode","The hydrostatic surface mode was not resolved."); end
                    surface = positive(1); positive(1) = [];
                end
                p = IMModeSelection.take(positive,request.nPositiveModes,"positive/internal");
                index = [negative(:);zero(:);surface(:);p(:)];
                zeroName = "zeroEigenvalue";
                if kind == "mda", zeroName = "constant"; end
                if kind == "hydrostaticF", zeroName = "barotropic"; end
                branch = [repmat("negative",numel(negative),1);repmat(zeroName,numel(zero),1);repmat("surface",numel(surface),1);repmat("positive",numel(p),1)];
                number = [(1:numel(negative))';zeros(numel(zero),1);ones(numel(surface),1);(1:numel(p))'];
                labels = [-1:-1:-numel(negative),zeros(1,numel(zero)),1:(numel(surface)+numel(p))];
            end
            columnLabel = branch+":"+string(number);
            info = table(columnLabel,branch,number,lambda(index),VariableNames=["columnLabel","branch","branchIndex","eigenvalue"]);
            selection = struct("sortIndex",index(:),"modeNumber",labels,"modeSelectionDiagnostics",diagnostics,"modeInfo",info);
        end
    end
    methods (Static, Access = private)
        function selected = take(candidates,count,branch)
            if numel(candidates) < count
                error("IMModes:InsufficientModes","Requested %d %s modes, but only %d were available. Increase solver resolution or revise the request.",count,branch,numel(candidates));
            end
            selected = candidates(1:count);
        end
    end
end

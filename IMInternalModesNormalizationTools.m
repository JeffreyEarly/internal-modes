classdef (Hidden, Sealed) IMInternalModesNormalizationTools
    % Define normalization availability and defaults for internal-mode EVPs.

    methods (Static)
        function names = available(evp)
            names = ["unity","uMax","wMax","depthRMSF","depthRMSG","maxAbsF","maxAbsG","generalizedEnergy"];
            if IMInternalModesNormalizationTools.isWave(evp) && isfield(evp.parameters,"k")
                names(end+1) = "waveEnergy";
            end
            if evp.modeFamily ~= "meanDensityAnomaly"
                names(end+1) = "surfacePressure";
            end
            if evp.modeFamily == "hydrostatic"
                names(end+1) = "geostrophic";
                if evp.innerProduct("F").hasInnerProduct
                    names(end+1) = "depth";
                end
            end
            if string(evp.name) == "waveModesAtWavenumber"
                names(end+1) = "kConstant";
            end
            if string(evp.name) == "geostrophicGeneralizedPotentialEnstrophyModes"
                names(end+1) = "generalizedPotentialEnstrophy";
            end
        end

        function name = default(evp)
            if string(evp.name) == "geostrophicGeneralizedPotentialEnstrophyModes"
                name = "generalizedPotentialEnstrophy";
            elseif string(evp.name) == "geostrophicAPVModes"
                name = "depthRMSF";
            elseif evp.modeFamily == "hydrostatic"
                name = "depthRMSF";
            elseif string(evp.name) == "waveModesAtWavenumber"
                name = "waveEnergy";
            elseif string(evp.name) == "waveModesAtFrequency"
                name = "depthRMSF";
            else
                name = "unity";
            end
        end

        function tf = isWave(evp)
            tf = ismember(string(evp.name),["waveModesAtWavenumber","waveModesAtFrequency"]);
        end

        function tf = hasHydrostaticDerivativeRelations(evp)
            tf = ismember(string(evp.name),["hydrostaticGModes","hydrostaticFModes", ...
                "geostrophicAPVModes","geostrophicGeneralizedPotentialEnstrophyModes"]);
        end

        function info = modeInfo(evp,modeNumber,eigenvalues,metadata,diagnostics)
            if isfield(diagnostics,"modeInfo")
                info = diagnostics.modeInfo;
            elseif isfield(metadata,"modeInfo")
                info = metadata.modeInfo;
            else
                branchIndex = abs(modeNumber(:));
                branch = repmat("positive",numel(modeNumber),1);
                branch(modeNumber < 0) = "negative";
                branch(modeNumber == 0) = "zeroEigenvalue";
                if startsWith(string(evp.name),"hydrostatic") || isfield(metadata,"includeBarotropic")
                    branch(modeNumber == 0) = "barotropic";
                elseif evp.modeFamily == "meanDensityAnomaly"
                    branch(modeNumber == 0) = "constant";
                end
                info = table(string(modeNumber(:)),branch,branchIndex,eigenvalues(:),VariableNames=["columnLabel","branch","branchIndex","eigenvalue"]);
            end
        end

        function frequency = frequency(evp,h)
            if isfield(evp.parameters,"omega")
                frequency = repmat(evp.parameters.omega,size(h));
            elseif IMInternalModesNormalizationTools.isWave(evp) && isfield(evp.parameters,"k")
                frequency = sqrt(evp.f0^2+evp.g*evp.parameters.k^2*h);
            else
                frequency = zeros(size(h));
            end
        end

        function access = normalizationCache()
            % Immutable basis copies share only values derived from identical state.
            % Mutating rules or selecting columns installs a new closure.
            values = configureDictionary("string","cell");
            access = @cached;
            function value = cached(name,varargin)
                if ~isempty(varargin)
                    values{name} = varargin{1};
                end
                if isKey(values,name)
                    value = values{name};
                else
                    value = [];
                end
            end
        end

        function factors = maximumAmplitudes(evaluator,z,values)
            [z,index] = sort(z);
            values = values(index,:);
            factors = max(abs(values),[],1);
            settings = optimset("TolX",max(1e-12,eps(diff([z(1),z(end)]))),"Display","off");
            for j = 1:size(values,2)
                magnitude = abs(values(:,j));
                peaks = find(magnitude(2:end-1) >= magnitude(1:end-2) & magnitude(2:end-1) >= magnitude(3:end) & (magnitude(2:end-1) > magnitude(1:end-2) | magnitude(2:end-1) > magnitude(3:end)))+1;
                for index = peaks(:).'
                    [~,negativePeak] = fminbnd(@(zz) -IMInternalModesNormalizationTools.columnMagnitude(evaluator,zz,j),z(index-1),z(index+1),settings);
                    factors(j) = max(factors(j),-negativePeak);
                end
            end
        end

        function value = columnMagnitude(evaluator,z,j)
            values = evaluator(z);
            value = abs(values(:,j));
        end

        function selected = columns(values,columns)
            selected = values(:,columns);
        end

        function factor = selectedNormFactor(basis,name,column)
            factors = basis.normalizationFactors(name);
            factor = factors(column);
        end
    end
end

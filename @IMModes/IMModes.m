classdef (Sealed) IMModes
    % Create continuous mode families from one authoritative background.
    %
    % Count options retain regular ladders. Every admissible finite special
    % sector is returned automatically and does not consume those counts.
    % Numerical resolution is independent of mode counts and physical sampling.
    %
    % ```matlab
    % background = IMModeBackground(N2=N2,zDomain=[-1000 0]);
    % solver = IMSolverSpectral(nEVP=128,coordinateKind="bandwidth");
    % waves = IMModes.waveModesAtWavenumber(background,1e-3,solver=solver,nInternalModes=8,nSubinertialModes=3);
    % ```
    %
    % - Topic: Create modes
    % - Declaration: classdef IMModes
    methods (Static)
        function collection = waveModesAtWavenumbers(background,k,options)
            % Preserve requested pages and solve each distinct positive k once.
            % Counts may be scalar or aligned with the requested k row. Repeated
            % k values must request identical ladders. Inertial modes are separate.
            % - Topic: Create modes
            arguments (Input)
                background (1,1) IMModeBackground
                k (1,:) double {mustBeReal,mustBeFinite,mustBePositive}
                options.solver = IMSolverSpectral()
                options.nInternalModes (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 8
                options.nSubinertialModes (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 0
                options.normalization (1,1) string = "waveEnergy"
            end
            if isempty(k), error("IMModes:EmptyWavenumbers","Supply at least one positive wavenumber."); end
            counts = {options.nInternalModes,options.nSubinertialModes};
            for j = 1:2
                if isscalar(counts{j})
                    counts{j} = repmat(counts{j},size(k));
                elseif numel(counts{j}) ~= numel(k)
                    error("IMModes:InvalidModeCounts","Each count must be scalar or contain one value per requested wavenumber.");
                end
            end
            [distinct,first,mapping] = unique(k,"stable");
            bases = cell(1,numel(distinct));
            for j = 1:numel(distinct)
                pages = k == distinct(j);
                if numel(unique(counts{1}(pages)))~=1 || numel(unique(counts{2}(pages)))~=1
                    error("IMModes:ConflictingModeCounts","Repeated wavenumbers must request identical ladder counts.");
                end
                bases{j} = IMModes.waveModesAtWavenumber(background,distinct(j),solver=options.solver,nInternalModes=counts{1}(first(j)),nSubinertialModes=counts{2}(first(j)),normalization=options.normalization);
            end
            collection = IMBasisCollection(bases,kappa=k,basisIndex=reshape(mapping,1,[]));
        end

        function basis = waveModesAtWavenumber(background,k,options)
            % Return the surface mode and independently requested wave ladders.
            % Internal frequencies decrease; subinertial frequencies increase.
            % Both profiles use the fixed definition F = D Gz.
            % - Topic: Create modes
            arguments (Input)
                background (1,1) IMModeBackground
                k (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
                options.solver = IMSolverSpectral()
                options.nInternalModes (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 8
                options.nSubinertialModes (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 0
                options.normalization (1,1) string = "waveEnergy"
            end
            args = IMModes.physicalArguments(background);
            boundary = IMModes.waveBoundary(background);
            evp = IMInternalModes.waveModesAtWavenumber(args{:},k=k,surfaceBoundary=boundary);
            request = struct("kind","wave","hasSurface",background.surfaceCondition=="freeSurface","nPositiveModes",options.nInternalModes,"nNegativeModes",options.nSubinertialModes);
            [source,columns,info,h,lambda] = IMModes.solve(background,evp,options.solver,request);
            frequency = sqrt(background.f0^2+background.g*k^2*h);
            basis = IMModes.wrap(source,columns,info,background,"wave",options.normalization,h,lambda,true,k,frequency);
        end

        function basis = inertialModes(background,options)
            % Return an auxiliary vertical basis with physical frequency |f0|.
            % Positive/negative counts refer to auxiliary equivalent depths.
            % - Topic: Create modes
            arguments (Input)
                background (1,1) IMModeBackground
                options.solver = IMSolverSpectral()
                options.nPositiveModes (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 8
                options.nNegativeModes (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 0
                options.normalization (1,1) string = "depthRMSF"
            end
            args = IMModes.physicalArguments(background);
            evp = IMInternalModes.waveModesAtWavenumber(args{:},k=0,surfaceBoundary=IMModes.waveBoundary(background));
            request = struct("kind","inertial","hasSurface",background.surfaceCondition=="freeSurface","nPositiveModes",options.nPositiveModes,"nNegativeModes",options.nNegativeModes);
            [source,columns,info,h,lambda] = IMModes.solve(background,evp,options.solver,request);
            basis = IMModes.wrap(source,columns,info,background,"inertial",options.normalization,h,lambda,true,0,repmat(abs(background.f0),size(h)));
        end

        function basis = waveModesAtFrequency(background,omega,options)
            % Return propagating and evanescent roots at one fixed frequency.
            % The inertial-frequency degeneracy has its own inertial factory.
            % - Topic: Create modes
            arguments (Input)
                background (1,1) IMModeBackground
                omega (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
                options.solver = IMSolverSpectral()
                options.nPropagatingModes (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 8
                options.nEvanescentModes (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 0
                options.normalization (1,1) string = "depthRMSF"
            end
            if omega == abs(background.f0)
                error("IMModes:InertialFrequency","Use inertialModes for omega=abs(f0); the fixed-frequency horizontal dispersion is degenerate.");
            end
            if options.normalization == "waveEnergy"
                error("IMModes:UnavailableNormalization","Fixed-frequency roots have distinct or evanescent horizontal wavenumbers. Use depthRMSF or another supported profile normalization.");
            end
            args = IMModes.physicalArguments(background);
            evp = IMInternalModes.waveModesAtFrequency(args{:},omega=omega,surfaceBoundary=IMModes.waveBoundary(background));
            request = struct("kind","frequency","hasSurface",background.surfaceCondition=="freeSurface","nPropagatingModes",options.nPropagatingModes,"nEvanescentModes",options.nEvanescentModes);
            [source,columns,info,h,lambda] = IMModes.solve(background,evp,options.solver,request);
            basis = IMModes.wrap(source,columns,info,background,"fixedFrequencyWave",options.normalization,h,lambda,true,0,repmat(omega,size(h)));
        end

        function basis = geostrophicAPVModes(background,options)
            % Return all negative/zero eigenvalues plus nPositiveModes.
            % A zero eigenvalue means h=Inf, not a zero-depth wave.
            % - Topic: Create modes
            arguments (Input)
                background (1,1) IMModeBackground
                options.solver = IMSolverSpectral()
                options.energyBoundary (1,1) IMGeneralizedEnergyBoundary
                options.nPositiveModes (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 8
                options.normalization (1,1) string = "depthRMSF"
            end
            args = IMModes.physicalArguments(background);
            energy = options.energyBoundary;
            evp = IMInternalModes.geostrophicAPVModes(args{:},g0=energy.g0,gd=energy.gd,surfaceBoundary=background.surfaceCondition);
            request = struct("kind","apv","nPositiveModes",options.nPositiveModes);
            [source,columns,info,h,lambda] = IMModes.solve(background,evp,options.solver,request);
            basis = IMModes.wrap(source,columns,info,background,"geostrophicAPV",options.normalization,h,lambda,false,0,zeros(size(h)));
        end

        function basis = meanDensityAnomalyModes(background,options)
            % Return signed boundary modes, the constant, and a positive ladder.
            % - Topic: Create modes
            arguments (Input)
                background (1,1) IMModeBackground
                options.solver = IMSolverSpectral()
                options.energyBoundary (1,1) IMGeneralizedEnergyBoundary
                options.nPositiveModes (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 8
                options.normalization (1,1) string = "generalizedEnergy"
            end
            args = IMModes.physicalArguments(background);
            energy = options.energyBoundary;
            evp = IMInternalModes.meanDensityAnomalyModes(args{:},g0=energy.g0,gd=energy.gd);
            request = struct("kind","mda","nPositiveModes",options.nPositiveModes);
            [source,columns,info,h,lambda] = IMModes.solve(background,evp,options.solver,request);
            basis = IMModes.wrap(source,columns,info,background,"meanDensityAnomaly",options.normalization,h,lambda,false,0,zeros(size(h)));
        end

        function basis = geostrophicGeneralizedPotentialEnstrophyModes(background,k,options)
            % Create the generalized-potential-enstrophy eigenvalue ladder.
            % - Topic: Create modes
            arguments (Input)
                background (1,1) IMModeBackground
                k (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
                options.solver = IMSolverSpectral()
                options.nModes (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 8
                options.alpha0 double {mustBeReal} = []
                options.alphaD (1,1) double {mustBeReal} = Inf
                options.normalization (1,1) string = "generalizedPotentialEnstrophy"
            end
            if background.surfaceCondition ~= "freeSurface"
                error("IMModes:UnsupportedSurfaceCondition","Generalized-potential-enstrophy modes currently require a free surface.");
            end
            args = IMModes.physicalArguments(background);
            evp = IMInternalModes.geostrophicGeneralizedPotentialEnstrophyModes(args{:},k=k,alpha0=options.alpha0,alphaD=options.alphaD);
            request = struct("kind","generalizedEnstrophy","nPositiveModes",options.nModes);
            [source,columns,info,h,lambda] = IMModes.solve(background,evp,options.solver,request);
            basis = IMModes.wrap(source,columns,info,background,"geostrophicGeneralizedPotentialEnstrophy",options.normalization,h,lambda,false,k,zeros(size(h)));
        end

        function basis = hydrostaticGModes(background,options)
            % Return the external surface mode and nInternalModes G modes.
            % - Topic: Create modes
            arguments (Input)
                background (1,1) IMModeBackground
                options.solver = IMSolverSpectral()
                options.nInternalModes (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 8
                options.normalization (1,1) string = "depthRMSF"
            end
            args = IMModes.physicalArguments(background);
            evp = IMInternalModes.hydrostaticGModes(args{:},surfaceBoundary=IMModes.waveBoundary(background));
            request = struct("kind","hydrostaticG","hasSurface",background.surfaceCondition=="freeSurface","nPositiveModes",options.nInternalModes);
            [source,columns,info,h,lambda] = IMModes.solve(background,evp,options.solver,request);
            basis = IMModes.wrap(source,columns,info,background,"hydrostaticG",options.normalization,h,lambda,false,0,zeros(size(h)));
        end

        function basis = hydrostaticFModes(background,options)
            % Return the barotropic/external sector and nInternalModes F modes.
            % - Topic: Create modes
            arguments (Input)
                background (1,1) IMModeBackground
                options.solver = IMSolverSpectral()
                options.nInternalModes (1,1) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 8
                options.normalization (1,1) string = "depthRMSF"
            end
            args = IMModes.physicalArguments(background);
            % The regular pair avoids inverse N2 at neutral depths under
            % either physical surface convention.
            evp = IMInternalModes.geostrophicAPVModes(args{:},g0=Inf,gd=Inf,surfaceBoundary=background.surfaceCondition);
            if background.surfaceCondition == "freeSurface"
                request = struct("kind","hydrostaticG","hasSurface",true,"nPositiveModes",options.nInternalModes);
            else
                request = struct("kind","hydrostaticF","hasSurface",false,"nPositiveModes",options.nInternalModes);
            end
            [source,columns,info,h,lambda] = IMModes.solve(background,evp,options.solver,request);
            basis = IMModes.wrap(source,columns,info,background,"hydrostaticF",options.normalization,h,lambda,false,0,zeros(size(h)));
        end

        function basis = geostrophicZeroAPVModesAtWavenumber(background,k,options)
            % Return every active endpoint response in canonical endpoint order.
            % Energy coefficients select active endpoints, not response shapes.
            % - Topic: Create modes
            arguments (Input)
                background (1,1) IMModeBackground
                k (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
                options.solver = IMSolverSpectral()
                options.energyBoundary (1,1) IMGeneralizedEnergyBoundary
                options.normalization (1,1) string = "endpointResponse"
            end
            endpoints = options.energyBoundary.activeEndpoints;
            info = IMModes.info(repmat("endpointResponse",numel(endpoints),1),ones(numel(endpoints),1));
            info.columnLabel = endpoints(:);
            info.endpoint = endpoints(:);
            if isempty(endpoints)
                source = [];
            elseif isa(options.solver,"IMSolverAnalytical")
                solution = background.analyticalSolution();
                source = solution.geostrophicZeroAPVModesAtWavenumber(k,endpoints=endpoints,surfaceBoundary=background.surfaceCondition);
            elseif isa(options.solver,"IMSolver")
                args = IMModes.physicalArguments(background);
                problem = IMGeostrophicZeroAPVModes.atWavenumber(args{:},k=k,endpoints=endpoints,surfaceBoundary=background.surfaceCondition);
                source = options.solver.solveGeostrophicZeroAPVModes(problem);
            else
                error("IMModes:InvalidSolver","Use an IMSolver or IMSolverAnalytical backend.");
            end
            basis = IMModes.wrap(source,1:numel(endpoints),info,background,"geostrophicZeroAPV",options.normalization,[],[],false,k,zeros(1,numel(endpoints)));
        end

        function basis = neutralLabelModes(background,options)
            % Return analytical labels and any fully neutral constant kernel.
            % Intervals are caller declarations of exact N2=0, never inferred.
            % - Topic: Create modes
            arguments (Input)
                background (1,1) IMModeBackground
                options.neutralIntervals (:,2) double {mustBeReal,mustBeFinite}
                options.energyBoundary (1,1) IMGeneralizedEnergyBoundary
                options.nModesPerInterval (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 8
                options.normalization (1,1) string = "depthRMSG"
            end
            energy = options.energyBoundary;
            IMModes.validateIntervals(background,options.neutralIntervals,0,"neutral");
            source = IMNeutralLabelModesBasis(options.neutralIntervals,zDomain=background.zDomain,nModes=options.nModesPerInterval,rotationRate=background.rotationRate,latitude=background.latitude,g0=energy.g0,gd=energy.gd);
            branch = repmat("neutralLabel",numel(source.m_n),1);
            branch(source.m_n==0) = "constantLabel";
            info = IMModes.info(branch,source.localModeNumber(:));
            info.columnLabel = source.columnLabels(:);
            info.intervalIndex = source.intervalIndex(:);
            info.localModeNumber = source.localModeNumber(:);
            basis = IMModes.wrap(source,1:numel(source.m_n),info,background,"neutralLabel",options.normalization,[],[],false,0,zeros(size(source.m_n)));
        end

        function basis = zeroDepthWaveModesAtWavenumber(background,k,options)
            % Create a finite DD auxiliary basis on exact N2=f0^2 intervals.
            % This infinite-dimensional wave kernel needs explicit truncation.
            % - Topic: Create modes
            arguments (Input)
                background (1,1) IMModeBackground
                k (1,1) double {mustBeReal,mustBeFinite,mustBePositive}
                options.criticalIntervals (:,2) double {mustBeReal,mustBeFinite}
                options.nModesPerInterval (1,:) double {mustBeInteger,mustBeNonnegative,mustBeFinite} = 8
                options.normalization (1,1) string = "waveEnergy"
            end
            IMModes.validateIntervals(background,options.criticalIntervals,background.f0^2,"critical");
            % Nonzero endpoint weights enforce DD at physical boundaries too.
            auxiliary = IMNeutralLabelModesBasis(options.criticalIntervals,zDomain=background.zDomain,nModes=options.nModesPerInterval,rotationRate=background.rotationRate,latitude=background.latitude,g0=1,gd=1);
            info = IMModes.info(repmat("zeroDepth",numel(auxiliary.m_n),1),auxiliary.localModeNumber(:));
            info.intervalIndex = auxiliary.intervalIndex(:);
            info.columnLabel = "zeroDepth:"+auxiliary.columnLabels(:);
            evaluator = @(variable,z,order) IMModes.zeroDepthValues(auxiliary,background.D,variable,z,order);
            basis = IMModeBasis(background=background,family="zeroDepthWave",modeInfo=info,normalization=options.normalization,waveConvention=true,k=k,h=zeros(size(auxiliary.m_n)),hasEquivalentDepths=true,frequency=repmat(abs(background.f0),size(auxiliary.m_n)),intrinsicEvaluator=evaluator,variables=["F","G"],integrationIntervals=options.criticalIntervals);
        end
    end

    methods (Static, Access = private)
        function argumentsCell = physicalArguments(background)
            options = struct("N2",background.N2,"zDomain",background.zDomain,"rotationRate",background.rotationRate,"latitude",background.latitude,"g",background.g);
            argumentsCell = namedargs2cell(options);
        end

        function boundary = waveBoundary(background)
            if background.surfaceCondition == "freeSurface"
                boundary = IMBoundaryCondition(a=0,b=1,c=1,d=0);
            else
                boundary = IMBoundaryCondition.dirichlet();
            end
        end

        function [source,columns,info,h,lambda] = solve(background,evp,solver,request)
            if isa(solver,"IMSolverAnalytical")
                solution = background.analyticalSolution();
                if ismember(request.kind,["mda","generalizedEnstrophy"])
                    error("IMModes:UnsupportedAnalyticalFamily","The analytical backend does not implement %s. Choose a numerical solver.",request.kind);
                end
                if request.kind == "frequency"
                    count = request.nPropagatingModes+request.nEvanescentModes+3;
                else
                    count = request.nPositiveModes+3;
                    if isfield(request,"nNegativeModes"), count = count+request.nNegativeModes; end
                end
                source = solution.internalModes(evp,nModes=max(1,count));
                zero = find(source.modeNumber==0);
                selection = IMModeSelection.select(evp,source.eigenvalues(:),request,zeroIndices=zero);
                columns = selection.sortIndex(:).';
            elseif isa(solver,"IMSolver")
                [source,selection] = solver.solveRequestedEVP(evp,request);
                columns = 1:numel(selection.sortIndex);
            else
                error("IMModes:InvalidSolver","Use an IMSolver or IMSolverAnalytical backend.");
            end
            info = selection.modeInfo;
            if isempty(columns)
                source = [];
                h = zeros(1,0); lambda = zeros(1,0);
            else
                h = source.h(columns); lambda = source.eigenvalues(columns);
            end
            info.eigenvalue = lambda(:);
        end

        function basis = wrap(source,columns,info,background,family,normalization,h,lambda,wave,k,frequency)
            if isempty(source)
                evaluator = @(~,z,~) zeros(numel(z),height(info));
                hasSpectrum = ~ismember(family,["neutralLabel","geostrophicZeroAPV"]);
                basis = IMModeBasis(background=background,family=family,modeInfo=info,normalization=normalization,h=h,eigenvalues=lambda,hasEquivalentDepths=hasSpectrum,hasEigenvalues=hasSpectrum,waveConvention=wave,k=k,frequency=frequency,intrinsicEvaluator=evaluator,variables=["F","G"]);
            else
                basis = IMModeBasis(source=source,background=background,family=family,columns=columns,modeInfo=info,normalization=normalization,h=h,eigenvalues=lambda,waveConvention=wave,k=k,frequency=frequency);
            end
        end

        function info = info(branch,number)
            columnLabel = branch+":"+string(number);
            info = table(columnLabel,branch,number,VariableNames=["columnLabel","branch","branchIndex"]);
        end

        function validateIntervals(background,intervals,value,kind)
            if isempty(intervals) || any(intervals(:,1)>=intervals(:,2)) || any(intervals(:,1)<background.zDomain(1)) || any(intervals(:,2)>background.zDomain(2)) || any(intervals(2:end,1)<=intervals(1:end-1,2))
                error("IMModes:InvalidIntervals","Declare ordered, separated nonempty %s intervals inside zDomain.",kind);
            end
            % Reject contradictions to the supplied profile. Sampling cannot
            % prove an arbitrary function is exactly constant; the declaration
            % remains the caller's scientific assertion, never an inference.
            for j = 1:size(intervals,1)
                z = linspace(intervals(j,1),intervals(j,2),33).';
                values = background.N2(z);
                if ~isequal(size(values),size(z)) || ~isreal(values) || any(~isfinite(values)) || any(values ~= value)
                    error("IMModes:InconsistentInterval","The declared %s interval contradicts the background N2 profile. Exact N2=%g is required throughout each interval.",kind,value);
                end
            end
        end

        function values = zeroDepthValues(auxiliary,D,variable,z,order)
            if variable == "G"
                switch order
                    case 0, values = auxiliary.G(z);
                    case 1, values = auxiliary.Gz(z);
                    case 2, values = auxiliary.Gzz(z);
                    otherwise, error("IMModeBasis:UnsupportedDerivative","Only G derivatives through order two are available.");
                end
            elseif variable == "F" && order <= 1
                if order == 0, values = D*auxiliary.Gz(z); else, values = D*auxiliary.Gzz(z); end
            else
                error("IMModeBasis:UnsupportedVariable","Zero-depth waves provide F and G profiles.");
            end
        end
    end
end

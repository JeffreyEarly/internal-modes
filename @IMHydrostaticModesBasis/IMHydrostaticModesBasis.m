classdef IMHydrostaticModesBasis < IMInternalModesBasis
    % Retain one rigid-lid G solve with its analytical barotropic F mode.
    %
    % The native columns are $$G_0=0,G_1,\ldots,G_n$$. Their aligned
    % horizontal profiles obey $$F_j=h_jG_j'$$ for positive modes and
    % $$F_0=1$$. The zero native column has $$h_0=\infty$$ and is inactive
    % in the G pairing. Native coefficients, solver diagnostics, custom
    % normalizations and Gram matrices remain available through the parent.
    %
    % ```matlab
    % source = solver.solveEVP(IMInternalModes.hydrostaticGModes(background),nModes=12);
    % basis = IMHydrostaticModesBasis(source,normalization="depthRMSF");
    % ```
    %
    % - Topic: Create internal-mode bases
    % - Topic: Evaluate modes
    % - Topic: Build discrete transforms
    % - Declaration: classdef IMHydrostaticModesBasis < IMInternalModesBasis

    methods
        function self = IMHydrostaticModesBasis(source,options)
            % Append the barotropic sector without solving another EVP.
            %
            % - Topic: Create internal-mode bases
            % - Parameter source: rigid-lid hydrostatic G basis of positive modes
            % - Parameter options.normalization: shared F/G normalization
            % - Returns self: aligned native hydrostatic basis
            arguments
                source (1,1) IMInternalModesBasis
                options.normalization = "depthRMSF"
            end
            surface = source.evp.surfaceBoundary;
            bottom = source.evp.bottomBoundary;
            if source.evp.name ~= "hydrostaticGModes" || source.evp.formulation ~= "G" || any(source.h <= 0 | ~isfinite(source.h)) || surface.a == 0 || surface.b ~= 0 || surface.c ~= 0 || surface.d ~= 0 || bottom.a == 0 || bottom.b ~= 0 || bottom.c ~= 0 || bottom.d ~= 0
                error("IMHydrostaticModesBasis:InvalidSource","Supply positive rigid-lid hydrostatic G modes with Dirichlet conditions at both endpoints.");
            end
            metadata = source.metadata;
            metadata.family = "hydrostaticFModes";
            metadata.includeBarotropic = true;
            metadata.profileConvention = "F=hGz with F0=1 and G0=0";
            info = source.modeInfo;
            metadata.modeInfo = table(["barotropic:0";info.columnLabel],["barotropic";info.branch],[0;info.branchIndex],[0;info.eigenvalue],VariableNames=["columnLabel","branch","branchIndex","eigenvalue"]);
            diagnostics = source.modeSelectionDiagnostics;
            diagnostics.modeInfo = metadata.modeInfo;
            diagnostics.zeroModeStatus = "analyticalBarotropic";
            diagnostics.zeroModeCount = 1;
            self@IMInternalModesBasis(solver=source.solver,evp=source.evp,nativeModes=[zeros(size(source.nativeModes,1),1) source.nativeModes],eigenvalues=[0 source.eigenvalues],modeNumber=[0 source.modeNumber],modeSelectionDiagnostics=diagnostics,normalization=options.normalization,metadata=metadata);
        end

        function values = rawVariable(self,variable,z)
            % Evaluate the native G family and its aligned analytical F sector.
            %
            % - Topic: Evaluate modes
            % - Parameter variable: F or G
            % - Parameter z: physical coordinate column
            % - Returns values: unnormalized sample-by-column profiles
            arguments
                self IMHydrostaticModesBasis
                variable (1,1) string {mustBeMember(variable,["F","G"])}
                z (:,1) double {mustBeReal,mustBeFinite}
            end
            if variable == "G"
                values = self.solver.evaluateNativeModes(self.nativeModes,z);
                values(z == self.zDomain(1) | z == self.zDomain(2),:) = 0;
            else
                active = self.modeNumber ~= 0;
                values = ones(numel(z),numel(self.modeNumber));
                if any(active)
                    values(:,active) = self.solver.evaluatePhysicalDerivative(self.nativeModes(:,active),z,1).*self.h(active);
                end
            end
        end

        function [z,design] = modeRootPoints(self,generatingMode)
            % Find endpoints and G roots from an already retained spectral mode.
            %
            % - Topic: Build discrete transforms
            % - Parameter generatingMode: positive physical mode number
            % - Returns z: increasing endpoints and interior G roots
            % - Returns design: point-construction provenance
            arguments
                self IMHydrostaticModesBasis
                generatingMode (1,1) double {mustBeInteger,mustBePositive}
            end
            if ~isa(self.solver,"IMSolverSpectral")
                error("IMHydrostaticModesBasis:UnavailableModeRoots","Mode-root points require a spectral hydrostatic basis.");
            end
            column = find(self.modeNumber == generatingMode,1);
            if isempty(column)
                error("IMHydrostaticModesBasis:MissingGeneratingMode","Retain internal mode %d before requesting its roots.",generatingMode);
            end
            roots = sort(self.solver.rootsOfNativeMode(self.nativeModes(:,column)));
            tolerance = max(100*eps(max(1,max(abs(self.zDomain)))),1e-10*max(1,diff(self.zDomain)));
            roots = roots(roots > self.zDomain(1)+tolerance & roots < self.zDomain(2)-tolerance);
            roots = roots(:);
            if ~isempty(roots), roots = roots([true;diff(roots) > tolerance]); end
            z = [self.zDomain(1);roots;self.zDomain(2)];
            design = IMDiscreteTransformTools.modeRootGridDesign(self,z,column-1,generatingMode);
        end

        function [z,design] = pointsFromModeRoots(self,options)
            % Include the analytical barotropic column in auxiliary root solves.
            %
            % - Topic: Build discrete transforms
            % - Parameter options.nModes: represented leading column count
            % - Returns z: endpoints and roots of the next native G column
            % - Returns design: point-construction provenance
            arguments
                self IMHydrostaticModesBasis
                options.nModes (1,1) double {mustBeInteger,mustBePositive} = self.nColumns
            end
            nModes = options.nModes;
            if nModes ~= self.nColumns
                [z,design] = pointsFromModeRoots@IMBasisSet(self,nModes=nModes);
                return
            end
            try
                source = self.solver.solveEVP(self.evp,nModes=nModes);
                auxiliary = IMHydrostaticModesBasis(source);
            catch exception
                wrapped = MException("IMBasisSet:AuxiliaryModeUnavailable", ...
                    "The solver could not obtain the auxiliary mode required for a %d-mode quadrature grid.",nModes);
                throwAsCaller(addCause(wrapped,exception))
            end
            candidateEigenvalues = auxiliary.eigenvalues(1:nModes);
            scale = max(1,max(abs([self.eigenvalues candidateEigenvalues]),[],"all"));
            if ~isequal(auxiliary.modeNumber(1:nModes),self.modeNumber) || any(abs(candidateEigenvalues-self.eigenvalues) > 1e-8*scale)
                error("IMBasisSet:AuxiliaryModeMismatch","The auxiliary solve did not reproduce the retained mode prefix used by this basis set.");
            end
            [z,design] = auxiliary.pointsFromModeRoots(nModes=nModes);
        end

        function [transform,assessment] = discreteTransform(self,varargin)
            % Build a fixed-grid transform or use the inherited qualification API.
            %
            % A positional IMVerticalGrid uses its points and weights exactly,
            % without a solve or weight fit. Name-value requests retain the
            % parent's full transform qualification interface.
            %
            % - Topic: Build discrete transforms
            % - Parameter varargin: one IMVerticalGrid or inherited name-value options
            % - Returns transform: aligned F/G transform with inactive G0
            % - Returns assessment: inherited qualification result, empty for a fixed grid
            if isempty(varargin) || ~isa(varargin{1},"IMVerticalGrid")
                [transform,assessment] = discreteTransform@IMInternalModesBasis(self,varargin{:});
                return
            end
            if numel(varargin) ~= 1 || ~isscalar(varargin{1})
                error("IMHydrostaticModesBasis:InvalidGrid","Supply one scalar IMVerticalGrid.");
            end
            grid = varargin{1};
            if ~isequal(grid.zDomain,self.zDomain)
                error("IMHydrostaticModesBasis:IncompatibleDomain","The basis and sampling grid must have the same zDomain.");
            end
            F = self.F(grid.z); G = self.G(grid.z);
            if ismember(string(self.normalization),["depth","depthRMSF"])
                HF = diff(self.zDomain)*ones(1,numel(self.h));
            else
                HF = diag(self.gramMatrix(variable="F")).';
            end
            HG = HF./self.h;
            HG(~isfinite(self.h)) = 0;
            activeG = isfinite(self.h);
            N2 = self.N2(grid.z);
            channel = @(interior,target,active) struct(available=true,reason="",activeModeMask=active,interiorWeight=interior,endpointMetricMatrix=zeros(grid.nPoints),metricMatrix=diag(interior.*grid.weights),targetGramMatrix=diag(target),targetMajorantGramMatrix=diag(target));
            channels = struct(F=channel(ones(size(grid.z)),HF,true(1,numel(self.h))),G=channel(N2/self.evp.g,HG,activeG));
            metadata = struct(construction="IMHydrostaticModesBasis.discreteTransform",surfaceCondition="rigidLid",columnLabels=self.columnLabels,relation="F = h dG/dz with F0=1 and G0=0");
            transform = IMInternalModesDiscreteTransform(z=grid.z,weights=grid.weights,modeNumber=self.modeNumber,h=self.h,normalization=string(self.normalization),inverseF=F,inverseG=G,endpointF=self.F(self.zDomain([2 1]).'),endpointG=self.G(self.zDomain([2 1]).'),channelData=channels,primaryVariable="G",zDomain=self.zDomain,g=self.evp.g,modeFamily="hydrostatic",N2Values=N2,problemMetadata=metadata,basis=self);
            assessment = [];
        end
    end
end

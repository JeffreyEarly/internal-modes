classdef IMInternalModes < IMEigenvalueProblem
    % Describe canonical EVPs with internal-mode interpretation.
    %
    % `IMInternalModes` translates standard `F` and `G` internal-mode
    % problems into the canonical scalar EVP. The solved scalar `u` is
    % either `F` or `G`; the other variable is recovered diagnostically by
    % the resulting `IMInternalModesBasis`.
    %
    % - Topic: Create internal-mode EVPs
    % - Topic: Inspect internal-mode EVPs
    % - Declaration: classdef IMInternalModes < IMEigenvalueProblem

    properties (SetAccess = private)
        % Solved physical variable, `"F"` or `"G"`.
        %
        % - Topic: Inspect internal-mode EVPs
        formulation = "G"

        % Coriolis parameter.
        %
        % - Topic: Inspect internal-mode EVPs
        f0 = 0

        % Gravitational acceleration.
        %
        % - Topic: Inspect internal-mode EVPs
        g = 9.81
    end

    methods
        function self = IMInternalModes(options)
            % Create an internal-mode canonical EVP.
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes(options)
            % - Parameter options.name: short EVP name
            % - Parameter options.formulation: solved variable, `"F"` or `"G"`
            % - Parameter options.p: canonical derivative-flux coefficient
            % - Parameter options.q: canonical left-side value coefficient
            % - Parameter options.r: canonical metric coefficient
            % - Parameter options.surfaceBoundary: surface endpoint condition
            % - Parameter options.bottomBoundary: bottom endpoint condition
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.metadata: additional scalar parameters
            % - Returns evp: internal-mode EVP descriptor
            arguments
                options.name {mustBeTextScalar} = "internalModes"
                options.formulation {mustBeTextScalar} = "G"
                options.p = @(z,~) ones(size(z))
                options.q = @(z,~) zeros(size(z))
                options.r = @(z,~) ones(size(z))
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.hFromEigenvalue = @(lambda) 1 ./ lambda
                options.hasZeroMode (1,1) logical = false
                options.defaultNormalization = []
                options.normalizations struct = struct()
                options.metadata struct = struct()
            end

            formulation = IMInternalModes.validateFormulation(options.formulation);
            metadata = options.metadata;
            metadata.f0 = options.f0;
            metadata.g = options.g;
            metadata.formulation = formulation;

            self@IMEigenvalueProblem(name=options.name, p=options.p, q=options.q, r=options.r, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                hFromEigenvalue=options.hFromEigenvalue, hasZeroMode=options.hasZeroMode, ...
                defaultNormalization=options.defaultNormalization, ...
                normalizations=options.normalizations, metadata=metadata);
            self.formulation = formulation;
            self.f0 = options.f0;
            self.g = options.g;
        end

        function context = contextForSolver(self, solver)
            % Return the internal-mode coefficient context.
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: context = contextForSolver(evp,solver)
            % - Parameter solver: canonical solver
            % - Returns context: coefficient context
            context = contextForSolver@IMEigenvalueProblem(self, solver);
            context.f0 = self.f0;
            context.g = self.g;
            context.formulation = self.formulation;
        end

        function spec = innerProduct(self, variable)
            % Return the `F` or `G` inner-product recipe.
            %
            % - Topic: Inspect internal-mode EVPs
            % - Declaration: spec = innerProduct(evp,variable)
            % - Parameter variable: `"F"` or `"G"`
            % - Returns spec: struct with interior and endpoint metric terms
            variable = IMInternalModes.validateFormulation(variable);
            spec.variable = variable;
            if variable == "G"
                spec.interiorWeight = @(z,ctx) ctx.N2(z)/ctx.g;
            else
                spec.interiorWeight = @(z,~) ones(size(z));
            end
            if variable == self.formulation
                spec.surfaceWeights = self.endpointWeights("surface");
                spec.bottomWeights = self.endpointWeights("bottom");
            else
                spec.surfaceWeights = struct("location", {}, "coefficient", {}, "c", {}, "d", {});
                spec.bottomWeights = struct("location", {}, "coefficient", {}, "c", {}, "d", {});
            end
            spec.hasKnownBoundaryWeights = true;
        end

        function basisSet = makeBasisSet(self, solver, nativeModes, eigenvalues, h, modeNumber, index)
            % Create an internal-mode basis set.
            %
            % - Topic: Developer topics
            % - Declaration: basisSet = makeBasisSet(evp,solver,nativeModes,eigenvalues,h,modeNumber,index)
            % - Returns basisSet: internal-mode basis set
            % - Developer: true
            basisSet = IMInternalModesBasis(solver=solver, evp=self, nativeModes=nativeModes, ...
                eigenvalues=eigenvalues, h=h, modeNumber=modeNumber, index=index, ...
                zDomain=solver.zDomain, N2Function=@(z) solver.N2(z));
        end
    end

    methods (Static)
        function evp = hydrostaticGModes(options)
            % Create the hydrostatic `G` internal-mode EVP.
            %
            % The canonical scalar form is
            % $$-G''=\lambda N^2G/g.$$
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.hydrostaticGModes(options)
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceBoundary: surface endpoint condition
            % - Parameter options.bottomBoundary: bottom endpoint condition
            % - Returns evp: hydrostatic `G` EVP
            arguments
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            normalizations = IMInternalModes.standardNormalizations();
            evp = IMInternalModes(name="hydrostaticGModes", formulation="G", ...
                p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ctx.N2(z)/ctx.g, f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                defaultNormalization=Normalization.geostrophic, normalizations=normalizations);
        end

        function evp = hydrostaticFModes(options)
            % Create the hydrostatic `F` internal-mode EVP.
            %
            % The canonical scalar form is
            % $$-\partial_z(N^{-2}F_z)=\lambda F/g.$$
            % It declares the barotropic zero mode.
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.hydrostaticFModes(options)
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceBoundary: surface endpoint condition
            % - Parameter options.bottomBoundary: bottom endpoint condition
            % - Returns evp: hydrostatic `F` EVP
            arguments
                options.g (1,1) double {mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.neumann()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.neumann()
            end

            normalizations = IMInternalModes.standardNormalizations();
            evp = IMInternalModes(name="hydrostaticFModes", formulation="F", ...
                p=@(z,ctx) 1./ctx.N2(z), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) ones(size(z))/ctx.g, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                defaultNormalization=Normalization.geostrophic, normalizations=normalizations, ...
                hasZeroMode=true);
        end

        function evp = waveModesAtWavenumber(options)
            % Create the fixed-wavenumber wave-mode EVP.
            %
            % The canonical scalar form is
            % $$-G''+K^2G=\lambda(N^2-f_0^2)G/g.$$
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.waveModesAtWavenumber(options)
            % - Parameter options.k: horizontal wavenumber
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceBoundary: surface endpoint condition
            % - Parameter options.bottomBoundary: bottom endpoint condition
            % - Returns evp: fixed-wavenumber `G` EVP
            arguments
                options.k (1,1) double {mustBeNonnegative}
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            k = options.k;
            metadata = struct("k", k);
            normalizations = IMInternalModes.standardNormalizations();
            evp = IMInternalModes(name="waveModesAtWavenumber", formulation="G", ...
                p=@(z,~) ones(size(z)), q=@(z,~) k*k*ones(size(z)), ...
                r=@(z,ctx) (ctx.N2(z) - ctx.f0*ctx.f0)/ctx.g, ...
                f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                defaultNormalization=Normalization.kConstant, normalizations=normalizations, ...
                metadata=metadata);
        end

        function evp = waveModesAtFrequency(options)
            % Create the fixed-frequency wave-mode EVP.
            %
            % The canonical scalar form is
            % $$-G''=\lambda(N^2-\omega^2)G/g.$$
            %
            % - Topic: Create internal-mode EVPs
            % - Declaration: evp = IMInternalModes.waveModesAtFrequency(options)
            % - Parameter options.omega: wave frequency
            % - Parameter options.f0: Coriolis parameter
            % - Parameter options.g: gravitational acceleration
            % - Parameter options.surfaceBoundary: surface endpoint condition
            % - Parameter options.bottomBoundary: bottom endpoint condition
            % - Returns evp: fixed-frequency `G` EVP
            arguments
                options.omega (1,1) double {mustBeNonnegative}
                options.f0 (1,1) double = 0
                options.g (1,1) double {mustBePositive} = 9.81
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
            end

            omega = options.omega;
            metadata = struct("omega", omega);
            normalizations = IMInternalModes.standardNormalizations();
            evp = IMInternalModes(name="waveModesAtFrequency", formulation="G", ...
                p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
                r=@(z,ctx) (ctx.N2(z) - omega*omega)/ctx.g, ...
                f0=options.f0, g=options.g, ...
                surfaceBoundary=options.surfaceBoundary, bottomBoundary=options.bottomBoundary, ...
                defaultNormalization=Normalization.omegaConstant, normalizations=normalizations, ...
                metadata=metadata);
        end
    end

    methods (Static, Access = private)
        function formulation = validateFormulation(formulation)
            formulation = string(formulation);
            if formulation ~= "F" && formulation ~= "G"
                error("IMInternalModes:InvalidFormulation", ...
                    "formulation must be ""F"" or ""G"".");
            end
        end

        function normalizations = standardNormalizations()
            normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor(basisSet.evp.formulation, iMode);
            normalizations.geostrophic = @(basisSet,iMode) basisSet.geostrophicNormFactor(iMode);
            normalizations.kConstant = @(basisSet,iMode) basisSet.innerProductNormFactor("G", iMode);
            normalizations.omegaConstant = @(basisSet,iMode) basisSet.innerProductNormFactor("F", iMode);
            normalizations.wMax = @(basisSet,iMode) basisSet.maxAbsFactor("G", iMode);
            normalizations.uMax = @(basisSet,iMode) basisSet.maxAbsFactor("F", iMode);
            normalizations.surfacePressure = @(basisSet,iMode) basisSet.surfacePressureNormFactor(iMode);
        end
    end
end

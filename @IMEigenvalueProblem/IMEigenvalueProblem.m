classdef IMEigenvalueProblem
    % Describe a canonical scalar eigenvalue problem.
    %
    % `IMEigenvalueProblem` is the solver-independent description of the
    % scalar problem
    % $$-(p u')' + q u = \lambda r u,$$
    % with endpoint conditions
    % $$-[a_i u-b_i(pu')]=\lambda[c_i u-d_i(pu')].$$
    % The EVP owns the continuous problem: the physical interval, coefficient
    % functions, endpoint conditions, normalization rules, and diagnostic
    % definiteness checks. A solver owns only the numerical choices used to
    % discretize this problem.
    %
    % ```matlab
    % evp = IMEigenvalueProblem(zDomain=[-1 0], ...
    %     p=@(z,~) ones(size(z)), q=@(z,~) zeros(size(z)), ...
    %     r=@(z,~) ones(size(z)), ...
    %     surfaceBoundary=IMBoundaryCondition.dirichlet(), ...
    %     bottomBoundary=IMBoundaryCondition.dirichlet());
    % solver = IMSolverSpectral(nEVP=64);
    % basisSet = solver.solveEVP(evp,nModes=4);
    % ```
    %
    % - Topic: Create EVPs
    % - Topic: Define canonical coefficients
    % - Topic: Inspect EVP metadata
    % - Topic: Inspect inner products
    % - Topic: Inspect definiteness diagnostics
    % - Topic: Select modes
    % - Topic: Developer topics
    % - Declaration: classdef IMEigenvalueProblem

    properties (SetAccess = private)
        % Short EVP name.
        %
        % `name` appears in solver diagnostics and generated error messages.
        %
        % - Topic: Inspect EVP metadata
        name = "canonical"

        % Physical vertical domain.
        %
        % `zDomain` is sorted as `[bottom surface]` and defines the
        % interval on which the canonical EVP is posed.
        %
        % - Topic: Define canonical coefficients
        zDomain = [-1 0]

        % Coefficient multiplying the derivative flux.
        %
        % `p` defines the flux term in $$-(p u')'$$. It may be a scalar,
        % a vector on the solver grid, or a function handle with signature
        % `values = p(z,ctx)`.
        %
        % - Topic: Define canonical coefficients
        p = @(z,~) ones(size(z))

        % Coefficient multiplying the solved variable on the left side.
        %
        % `q` defines the multiplication term in $$q u$$. It may be a
        % scalar, a vector on the solver grid, or a function handle with
        % signature `values = q(z,ctx)`.
        %
        % - Topic: Define canonical coefficients
        q = @(z,~) zeros(size(z))

        % Metric coefficient multiplying the eigenvalue side.
        %
        % `r` defines the interior metric in $$\lambda r u$$ and in the
        % default scalar inner product. It may be a scalar, a vector on the
        % solver grid, or a function handle with signature
        % `values = r(z,ctx)`.
        %
        % - Topic: Define canonical coefficients
        r = @(z,~) ones(size(z))

        % Surface endpoint condition.
        %
        % `surfaceBoundary` stores the scalar endpoint condition at
        % `zDomain(2)` in canonical `IMBoundaryCondition` form.
        %
        % - Topic: Define canonical coefficients
        surfaceBoundary = IMBoundaryCondition.dirichlet()

        % Bottom endpoint condition.
        %
        % `bottomBoundary` stores the scalar endpoint condition at
        % `zDomain(1)` in canonical `IMBoundaryCondition` form.
        %
        % - Topic: Define canonical coefficients
        bottomBoundary = IMBoundaryCondition.dirichlet()

        % Whether the scalar problem declares a zero mode.
        %
        % When true, mode selection retains one eigenvalue near zero and
        % labels it with mode number `0`.
        %
        % - Topic: Select modes
        hasZeroMode = false

        % Natural default normalization.
        %
        % If a basis set is created without an explicit normalization,
        % `defaultNormalization` becomes the active `basisSet.normalization`.
        % Evaluated modes are always raw modes divided by a per-mode scale,
        % $$u_j^{\mathrm{out}}(z)=u_j^{\mathrm{raw}}(z)/s_j.$$
        %
        % - Topic: Inspect EVP metadata
        defaultNormalization = []

        % Additional EVP metadata.
        %
        % Fields are copied into the coefficient context so custom
        % coefficient functions can read scalar parameters without new
        % public properties. Standard internal-mode factories add fields
        % `f0`, `g`, and `formulation`; wave-mode factories also add `k` or
        % `omega`. User-supplied metadata fields are preserved on
        % `evp.metadata` and are visible to coefficient handles through
        % `ctx.fieldName`.
        %
        % - Topic: Inspect EVP metadata
        metadata = struct()
    end

    properties
        % Named normalization rules.
        %
        % Each field stores a function handle with signature
        % `scale = rule(basisSet,iMode)`. The returned value is the raw
        % scale $$s_j$$ for one mode, and basis-set evaluation divides all
        % variables for that mode by $$s_j$$. The `unity` rule is supplied
        % automatically when omitted. Internal-mode factories add rules for
        % `geostrophic`, `kConstant`, `omegaConstant`, `wMax`, `uMax`, and
        % `surfacePressure`.
        %
        % ```matlab
        % normalizations.unity = @(basisSet,iMode) ...
        %     basisSet.innerProductNormFactor(iMode);
        % evp = IMEigenvalueProblem(normalizations=normalizations, ...
        %     defaultNormalization=Normalization.unity);
        % ```
        %
        % - Topic: Inspect EVP metadata
        normalizations = struct()
    end

    methods
        function self = IMEigenvalueProblem(options)
            % Create a canonical scalar EVP.
            %
            % - Topic: Create EVPs
            % - Declaration: evp = IMEigenvalueProblem(options)
            % - Parameter options.name: short EVP name
            % - Parameter options.zDomain: physical vertical domain
            % - Parameter options.p: derivative-flux coefficient
            % - Parameter options.q: left-side value coefficient
            % - Parameter options.r: eigenvalue-side metric coefficient
            % - Parameter options.surfaceBoundary: surface endpoint condition
            % - Parameter options.bottomBoundary: bottom endpoint condition
            % - Parameter options.hasZeroMode: whether one zero mode should be retained
            % - Parameter options.defaultNormalization: natural normalization
            % - Parameter options.normalizations: named normalization handles
            % - Parameter options.metadata: additional scalar parameters
            % - Returns evp: canonical EVP descriptor
            arguments
                options.name {mustBeTextScalar} = "canonical"
                options.zDomain (1,2) double = [-1 0]
                options.p = @(z,~) ones(size(z))
                options.q = @(z,~) zeros(size(z))
                options.r = @(z,~) ones(size(z))
                options.surfaceBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.bottomBoundary (1,1) IMBoundaryCondition = IMBoundaryCondition.dirichlet()
                options.hasZeroMode (1,1) logical = false
                options.defaultNormalization = []
                options.normalizations struct = struct()
                options.metadata struct = struct()
            end

            self.name = string(options.name);
            self.zDomain = sort(options.zDomain);
            self.p = options.p;
            self.q = options.q;
            self.r = options.r;
            self.surfaceBoundary = options.surfaceBoundary;
            self.bottomBoundary = options.bottomBoundary;
            self.hasZeroMode = options.hasZeroMode;
            self.defaultNormalization = options.defaultNormalization;
            self.normalizations = IMEigenvalueProblem.resolveNormalizations(options.normalizations);
            self.metadata = options.metadata;
        end

        function [A, B] = assemble(self, solver)
            % Build the canonical matrix pair on a solver grid.
            %
            % Interior rows discretize
            % $$-(p u')' + q u = \lambda r u.$$
            % The surface and bottom rows are replaced by the endpoint
            % conditions using endpoint values of `p`, producing the matrix
            % pencil $$A q = \lambda B q$$. This method is mainly for solver
            % implementations, diagnostics, and external eigensolver
            % experiments; ordinary workflows call `solver.solveEVP`.
            %
            % - Topic: Developer topics
            % - Declaration: [A,B] = assemble(evp,solver)
            % - Parameter solver: canonical EVP solver
            % - Returns A: left matrix
            % - Returns B: right matrix
            % - Developer: true
            solver = solver.configuredForEVP(self);
            context = self.contextForSolver(solver);
            z = solver.zNative(:);
            [pValues, qValues, rValues] = self.coefficientValues(z, context);
            pzValues = solver.differentiateGridValues(pValues, 1);
            D0 = solver.physicalDerivativeMatrix(0);
            D1 = solver.physicalDerivativeMatrix(1);
            D2 = solver.physicalDerivativeMatrix(2);

            A = -diag(pValues)*D2 - diag(pzValues)*D1 + diag(qValues)*D0;
            B = diag(rValues)*D0;
            [A, B] = self.applyBoundaryRow(A, B, solver, "surface", self.surfaceBoundary, pValues);
            [A, B] = self.applyBoundaryRow(A, B, solver, "bottom", self.bottomBoundary, pValues);
        end

        function context = contextForSolver(self, solver)
            % Return the coefficient context for this EVP and solver.
            %
            % The context begins with `solver.context()`, adds `zDomain`,
            % then copies each field of `metadata`. Coefficient handles such
            % as `p(z,ctx)`, `q(z,ctx)`, and `r(z,ctx)` receive this struct.
            %
            % - Topic: Developer topics
            % - Declaration: context = contextForSolver(evp,solver)
            % - Parameter solver: canonical solver
            % - Returns context: coefficient context
            % - Developer: true
            context = solver.context();
            context.zDomain = self.zDomain;
            metadataFields = fieldnames(self.metadata);
            for iField = 1:numel(metadataFields)
                fieldName = metadataFields{iField};
                context.(fieldName) = self.metadata.(fieldName);
            end
        end

        function profile = coordinateProfile(~, coordinateKind)
            % Return fields needed by a solver coordinate map.
            %
            % Generic canonical EVPs are independent of stratification and
            % only support the physical `z` coordinate.
            %
            % - Topic: Developer topics
            % - Declaration: profile = coordinateProfile(evp,coordinateKind)
            % - Parameter coordinateKind: solver coordinate kind
            % - Returns profile: struct with coordinate resources
            % - Developer: true
            coordinateKind = string(coordinateKind);
            if coordinateKind ~= "z"
                error("IMEigenvalueProblem:UnsupportedCoordinateKind", ...
                    "Generic canonical EVPs support coordinateKind=""z"" only.");
            end
            profile.N2 = [];
            profile.dzLogN2 = [];
        end

        function spec = innerProduct(self, variable)
            % Return the scalar inner-product recipe.
            %
            % The canonical basis set uses `r` in the interior and the
            % endpoint metric terms implied by active endpoint conditions:
            % $$M_{ij}=\int r u_i u_j\,dz+
            % \sum_\ell \gamma_\ell L_\ell[u_i]L_\ell[u_j].$$
            % The returned struct has fields `variable`, `interiorWeight`,
            % `surfaceWeights`, `bottomWeights`, and
            % `hasKnownBoundaryWeights`.
            %
            % - Topic: Inspect inner products
            % - Declaration: spec = innerProduct(evp,variable)
            % - Parameter variable: scalar variable name; only `"u"` is accepted
            % - Returns spec: struct with interior and endpoint metric terms
            arguments
                self IMEigenvalueProblem
                variable {mustBeTextScalar} = "u"
            end

            variable = string(variable);
            if variable ~= "u"
                error("IMEigenvalueProblem:InvalidVariable", ...
                    "Canonical scalar inner products use variable ""u"".");
            end
            spec.variable = "u";
            spec.interiorWeight = self.r;
            spec.surfaceWeights = self.endpointWeights("surface");
            spec.bottomWeights = self.endpointWeights("bottom");
            spec.hasKnownBoundaryWeights = true;
        end

        function weights = endpointWeights(self, location)
            % Return endpoint metric terms implied by active conditions.
            %
            % For an active endpoint condition
            % $$-[a_i u-b_i(pu')]=\lambda[c_i u-d_i(pu')],$$
            % Yassin's indexing uses `z_1` for the bottom and `z_2` for the
            % surface, so
            % $$D_i=(-1)^{i+1}(a_i d_i-b_i c_i).$$
            % Each returned struct has fields `location`, `coefficient`,
            % `c`, and `d`, representing the endpoint metric contribution
            % $$D_i^{-1}(c_i u-d_i p u_z)^2.$$
            %
            % - Topic: Inspect inner products
            % - Declaration: weights = endpointWeights(evp,location)
            % - Parameter location: `"surface"`, `"bottom"`, or omitted for both endpoints
            % - Returns weights: endpoint metric terms
            arguments
                self IMEigenvalueProblem
                location {mustBeTextScalar} = "all"
            end

            switch string(location)
                case "surface"
                    weights = self.weightForBoundary("surface", self.surfaceBoundary);
                case "bottom"
                    weights = self.weightForBoundary("bottom", self.bottomBoundary);
                case "all"
                    weights = [self.weightForBoundary("surface", self.surfaceBoundary); ...
                        self.weightForBoundary("bottom", self.bottomBoundary)];
                otherwise
                    error("IMEigenvalueProblem:InvalidBoundaryLocation", ...
                        "Boundary location must be ""surface"", ""bottom"", or ""all"".");
            end
        end

        function value = metricIndex(self, options)
            % Count negative endpoint directions in the metric.
            %
            % `metricIndex` is the number of active endpoint metric weights
            % with negative coefficient, after the supplied sign tolerance.
            % It is the finite endpoint part of the indefinite metric index
            % used by `negativeEigenvalueBounds`.
            %
            % - Topic: Inspect inner products
            % - Declaration: value = metricIndex(evp,options)
            % - Parameter options.tolerance: scalar sign tolerance
            % - Returns value: number of active endpoint terms with negative metric weight
            arguments
                self IMEigenvalueProblem
                options.tolerance (1,1) double {mustBeNonnegative} = 0
            end

            weights = self.endpointWeights();
            value = 0;
            for iWeight = 1:numel(weights)
                if weights(iWeight).coefficient < -options.tolerance
                    value = value + 1;
                end
            end
        end

        function info = definitenessInfo(self, solver)
            % Check grid-level signs for the canonical coefficients.
            %
            % This diagnostic certifies the assembled finite-dimensional
            % problem on the solver grid. It does not claim continuum signs
            % between grid points. The returned struct includes sampled
            % minima `pMin`, `qMin`, `rMin`; sign flags `pPositive`,
            % `qNonnegative`, and `rPositive`; metric fields
            % `endpointWeights`, `metricIndex`, `metricPositive`, and
            % `hasDegenerateEndpointMetric`; numerator fields
            % `endpointNumeratorNegativeDirections`,
            % `endpointNumeratorNonnegative`, `interiorNonnegative`, and
            % `qNonnegativeCertified`; and status fields
            % `certificationLevel` and `reason`.
            %
            % - Topic: Inspect definiteness diagnostics
            % - Declaration: info = definitenessInfo(evp,solver)
            % - Parameter solver: canonical EVP solver
            % - Returns info: struct with sign, metric, and endpoint checks
            solver = solver.configuredForEVP(self);
            context = self.contextForSolver(solver);
            z = solver.zNative(:);
            [pValues, qValues, rValues] = self.coefficientValues(z, context);
            tolerance = 100*eps;
            pTol = IMEigenvalueProblem.signTolerance(pValues, tolerance);
            qTol = IMEigenvalueProblem.signTolerance(qValues, tolerance);
            rTol = IMEigenvalueProblem.signTolerance(rValues, tolerance);

            info.pMin = min(pValues);
            info.qMin = min(qValues);
            info.rMin = min(rValues);
            info.pPositive = all(isfinite(pValues)) && info.pMin > pTol;
            info.qNonnegative = all(isfinite(qValues)) && info.qMin >= -qTol;
            info.rPositive = all(isfinite(rValues)) && info.rMin > rTol;
            info.endpointWeights = self.endpointWeights();
            info.metricIndex = self.metricIndex(tolerance=0);
            info.metricPositive = info.rPositive && info.metricIndex == 0;
            info.hasDegenerateEndpointMetric = self.hasDegenerateEndpointMetric();
            info.endpointNumeratorNegativeDirections = self.endpointNegativeDirections();
            info.endpointNumeratorNonnegative = info.endpointNumeratorNegativeDirections == 0;
            info.interiorNonnegative = info.pPositive && info.qNonnegative;
            info.qNonnegativeCertified = info.interiorNonnegative && info.endpointNumeratorNonnegative;
            info.certificationLevel = "grid";
            info.reason = "Grid signs and scalar endpoint conditions were checked.";
            if info.hasDegenerateEndpointMetric
                info.certificationLevel = "unknown";
                info.reason = "An active endpoint determinant is numerically degenerate.";
            elseif ~(info.pPositive && info.rPositive && all(isfinite(qValues)))
                info.certificationLevel = "unknown";
                info.reason = "One or more coefficient samples are nonfinite or fail the required signs.";
            end
        end

        function bounds = negativeEigenvalueBounds(self, solver, A)
            % Bound negative eigenvalues using grid-level certification.
            %
            % The returned counts describe how many negative eigenvalues are
            % certified by the discretized canonical problem, rather than by
            % raw negative finite-real eigenvalues alone. The returned struct
            % includes `certificationLevel`, `metricIndex`,
            % `zeroEigenvalueStatus`, `minNegativeEigenvalueCount`,
            % `maxNegativeEigenvalueCount`, and `reason`.
            % `maxNegativeEigenvalueCount` may be the string `"unknown"`
            % when coefficient signs or endpoint determinants cannot be
            % certified on the grid.
            %
            % - Topic: Inspect definiteness diagnostics
            % - Declaration: bounds = negativeEigenvalueBounds(evp,solver,A)
            % - Parameter solver: canonical EVP solver
            % - Parameter A: assembled left matrix, used for the zero-eigenvalue check
            % - Returns bounds: struct with min/max counts and a reason
            arguments
                self IMEigenvalueProblem
                solver IMSolver
                A double = []
            end

            info = self.definitenessInfo(solver);
            bounds.certificationLevel = info.certificationLevel;
            bounds.metricIndex = info.metricIndex;
            bounds.zeroEigenvalueStatus = "unchecked";
            bounds.minNegativeEigenvalueCount = 0;
            bounds.maxNegativeEigenvalueCount = "unknown";
            bounds.reason = info.reason;

            if info.certificationLevel == "unknown"
                return;
            end

            if info.metricPositive && info.qNonnegativeCertified
                bounds.maxNegativeEigenvalueCount = 0;
                bounds.reason = "The grid metric is positive and the quadratic form is nonnegative.";
                return;
            end

            if info.qNonnegativeCertified && info.metricIndex > 0
                bounds.zeroEigenvalueStatus = self.zeroEigenvalueStatus(A);
                if bounds.zeroEigenvalueStatus == "absent"
                    bounds.minNegativeEigenvalueCount = info.metricIndex;
                    bounds.maxNegativeEigenvalueCount = info.metricIndex;
                    bounds.reason = "The left-definite endpoint metric has a certified finite index and zero is absent.";
                else
                    bounds.maxNegativeEigenvalueCount = info.metricIndex;
                    bounds.reason = "The metric index is certified, but the zero-eigenvalue check is not conclusive.";
                end
                return;
            end

            if info.metricPositive && info.interiorNonnegative
                bounds.maxNegativeEigenvalueCount = info.endpointNumeratorNegativeDirections;
                bounds.reason = "The metric is positive and only endpoint numerator directions can make the form negative.";
                return;
            end
        end

        function selection = selectModes(self, eigenvalues, nModes, solver, A)
            % Select and label retained finite-real eigenmodes.
            %
            % Certified negative-count bounds decide when raw negative
            % discrete eigenvalues should be retained. With a positive metric
            % and nonnegative quadratic form, negative discrete eigenvalues
            % are ignored during mode selection.
            %
            % - Topic: Developer topics
            % - Declaration: selection = selectModes(evp,eigenvalues,nModes,solver,A)
            % - Parameter eigenvalues: finite real candidate eigenvalues
            % - Parameter nModes: number of retained modes
            % - Parameter solver: canonical solver
            % - Parameter A: assembled left matrix
            % - Returns selection: selected indices and mode numbers
            % - Developer: true
            arguments
                self IMEigenvalueProblem
                eigenvalues (:,1) double
                nModes (1,1) double {mustBeInteger, mustBePositive}
                solver IMSolver
                A double
            end

            tolerance = self.eigenvalueTolerance(eigenvalues);
            bounds = self.negativeEigenvalueBounds(solver, A);
            negativeCount = nnz(eigenvalues < -tolerance);
            if isnumeric(bounds.maxNegativeEigenvalueCount)
                negativeCount = min(negativeCount, bounds.maxNegativeEigenvalueCount);
            end

            negativeCandidates = find(eigenvalues < -tolerance);
            [~, negativeSort] = sort(eigenvalues(negativeCandidates), "ascend");
            negativeIndex = negativeCandidates(negativeSort(1:min(negativeCount,numel(negativeSort))));

            zeroIndex = zeros(0,1);
            if self.hasZeroMode
                zeroCandidates = find(abs(eigenvalues) <= tolerance);
                if isempty(zeroCandidates)
                    [~, zeroCandidate] = min(abs(eigenvalues));
                    zeroCandidates = zeroCandidate;
                end
                [~, zeroSort] = sort(abs(eigenvalues(zeroCandidates)), "ascend");
                zeroIndex = zeroCandidates(zeroSort(1));
            end

            positiveCandidates = find(eigenvalues > tolerance);
            [~, positiveSort] = sort(eigenvalues(positiveCandidates), "ascend");
            positiveIndex = positiveCandidates(positiveSort);

            sortIndex = [negativeIndex(:); zeroIndex(:); positiveIndex(:)];
            sortIndex = sortIndex(1:min(nModes,numel(sortIndex)));
            modeNumber = zeros(1,numel(sortIndex));
            nNegative = nnz(ismember(sortIndex, negativeIndex));
            modeNumber(1:nNegative) = -1:-1:-nNegative;
            nextIndex = nNegative + 1;
            if ~isempty(zeroIndex) && nextIndex <= numel(modeNumber) && sortIndex(nextIndex) == zeroIndex
                modeNumber(nextIndex) = 0;
                nextIndex = nextIndex + 1;
            end
            modeNumber(nextIndex:end) = 1:(numel(modeNumber) - nextIndex + 1);

            selection.sortIndex = sortIndex;
            selection.modeNumber = modeNumber;
            selection.index = bounds;
        end

        function basisSet = makeBasisSet(self, solver, nativeModes, eigenvalues, modeNumber, index)
            % Create the solved scalar basis set for this EVP.
            %
            % - Topic: Developer topics
            % - Declaration: basisSet = makeBasisSet(evp,solver,nativeModes,eigenvalues,modeNumber,index)
            % - Returns basisSet: solved scalar basis set
            % - Developer: true
            basisSet = IMBasisSet(solver=solver, evp=self, nativeModes=nativeModes, ...
                eigenvalues=eigenvalues, modeNumber=modeNumber, index=index, zDomain=self.zDomain);
        end
    end

    methods (Access = protected)
        function [pValues, qValues, rValues] = coefficientValues(self, z, context)
            pValues = IMEigenvalueProblem.evaluateCoefficient(self.p, z, context);
            qValues = IMEigenvalueProblem.evaluateCoefficient(self.q, z, context);
            rValues = IMEigenvalueProblem.evaluateCoefficient(self.r, z, context);
            pValues = IMEigenvalueProblem.expandCoefficient(pValues, z, "p");
            qValues = IMEigenvalueProblem.expandCoefficient(qValues, z, "q");
            rValues = IMEigenvalueProblem.expandCoefficient(rValues, z, "r");
        end
    end

    methods (Access = private)
        function [A, B] = applyBoundaryRow(~, A, B, solver, location, boundary, pValues)
            index = solver.boundaryIndex(location);
            D0 = solver.physicalDerivativeMatrix(0);
            D1 = solver.physicalDerivativeMatrix(1);
            pEndpoint = pValues(index);
            A(index,:) = -boundary.a*D0(index,:) + boundary.b*pEndpoint*D1(index,:);
            B(index,:) = boundary.c*D0(index,:) - boundary.d*pEndpoint*D1(index,:);
        end

        function weight = weightForBoundary(~, location, boundary)
            weight = struct("location", {}, "coefficient", {}, "c", {}, "d", {});
            if ~boundary.isEigenvalueDependent()
                return;
            end
            coefficient = boundary.metricWeight(location);
            if ~isfinite(coefficient)
                return;
            end
            weight(1,1).location = string(location);
            weight(1,1).coefficient = coefficient;
            weight(1,1).c = boundary.c;
            weight(1,1).d = boundary.d;
        end

        function count = endpointNegativeDirections(self)
            endpoints = [
                struct("location", "surface", "boundary", self.surfaceBoundary)
                struct("location", "bottom", "boundary", self.bottomBoundary)
            ];
            count = 0;
            for iEndpoint = 1:numel(endpoints)
                boundary = endpoints(iEndpoint).boundary;
                location = endpoints(iEndpoint).location;
                if boundary.isEigenvalueDependent()
                    if self.isDegenerateEndpointMetric(location, boundary)
                        continue;
                    end
                    H = boundary.endpointNumeratorMatrix(location);
                    count = count + nnz(eig(0.5*(H + H.')) < -sqrt(eps)*max(1,norm(H,"fro")));
                elseif boundary.b ~= 0
                    beta = boundary.robinEnergyCoefficient(location);
                    count = count + double(beta < -sqrt(eps)*max(1,abs(beta)));
                end
            end
        end

        function tf = hasDegenerateEndpointMetric(self)
            tf = self.isDegenerateEndpointMetric("surface", self.surfaceBoundary) ...
                || self.isDegenerateEndpointMetric("bottom", self.bottomBoundary);
        end

        function tf = isDegenerateEndpointMetric(~, location, boundary)
            tf = boundary.isEigenvalueDependent() && ~isfinite(boundary.metricWeight(location));
        end

        function status = zeroEigenvalueStatus(~, A)
            if isempty(A)
                status = "unchecked";
                return;
            end
            singularValues = svd(A);
            scale = max(1,norm(A,"fro"));
            if min(singularValues) > 1e-10*scale
                status = "absent";
            else
                status = "possible";
            end
        end

        function tolerance = eigenvalueTolerance(~, eigenvalues)
            finiteScale = abs(eigenvalues(isfinite(eigenvalues)));
            if isempty(finiteScale)
                scale = 1;
            else
                scale = max(1,min(finiteScale(finiteScale > 0), [], "omitnan"));
                if isempty(scale) || ~isfinite(scale)
                    scale = max(1,max(finiteScale));
                end
            end
            tolerance = 1e-8*scale;
        end
    end

    methods (Static)
        function values = evaluateCoefficient(coefficient, z, context)
            % Evaluate a scalar, vector, or coefficient function.
            %
            % - Topic: Developer topics
            % - Declaration: values = IMEigenvalueProblem.evaluateCoefficient(coefficient,z,context)
            % - Developer: true
            if isa(coefficient, "function_handle")
                try
                    values = coefficient(z, context);
                catch
                    try
                        values = coefficient(z);
                    catch
                        values = coefficient(context);
                    end
                end
            else
                values = coefficient;
            end
        end
    end

    methods (Static, Access = private)
        function values = expandCoefficient(values, z, name)
            if isscalar(values)
                values = values*ones(size(z));
            end
            values = values(:);
            if length(values) ~= length(z)
                error("IMEigenvalueProblem:InvalidCoefficientSize", ...
                    "Coefficient %s must evaluate to a scalar or one value per grid point.", name);
            end
        end

        function tolerance = signTolerance(values, tau)
            tolerance = tau*max(1,max(abs(values(:))));
        end

        function normalizations = resolveNormalizations(normalizations)
            if ~isfield(normalizations, "unity")
                normalizations.unity = @(basisSet,iMode) basisSet.innerProductNormFactor(iMode);
            end
        end
    end
end

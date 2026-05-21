classdef IMIndexPolicy
    % Specify the expected eigenvalue index of an EVP.
    %
    % The index records negative, zero, and positive eigenvalue counts. Mode
    % numbers are separate physical labels: `-1` identifies a surface
    % boundary branch, `-2` identifies a bottom boundary branch, `0`
    % identifies a true null mode, and positive labels identify interior
    % baroclinic modes.
    %
    % ```matlab
    % boundaryConditions = IMBoundary.partialDepthPE(boundarySign="negative");
    % policy = IMIndexPolicy.fromBoundaryConditions(boundaryConditions);
    % ```
    %
    % - Topic: Create index policies
    % - Topic: Validate index counts
    % - Declaration: classdef IMIndexPolicy

    properties
        % Expected number of negative directions.
        %
        % - Topic: Validate index counts
        expectedNegativeCountValue = NaN

        % Expected number of zero directions.
        %
        % - Topic: Validate index counts
        expectedZeroCountValue = NaN

        % Eigenvalue tolerance for zero-index classification.
        %
        % - Topic: Validate index counts
        indexTolerance = 1e-10

        % Validation behavior.
        %
        % Values are `"error"`, `"warning"`, or `"none"`.
        %
        % - Topic: Validate index counts
        validationMode = "error"

        % Declared endpoint boundary-mode slots.
        %
        % These modes are selected before null and interior modes. Each
        % entry has a physical `modeNumber` and an eigenvalue `indexSign`.
        %
        % - Topic: Validate index counts
        boundaryModes = struct("modeNumber", {}, "indexSign", {})
    end

    methods
        function self = IMIndexPolicy(options)
            % Create an index policy.
            %
            % - Topic: Create index policies
            % - Declaration: policy = IMIndexPolicy(options)
            % - Parameter options.expectedNegativeCount: expected negative count
            % - Parameter options.expectedZeroCount: expected zero count
            % - Parameter options.indexTolerance: tolerance for zero classification
            % - Parameter options.validationMode: `"error"`, `"warning"`, or `"none"`
            % - Parameter options.boundaryModes: endpoint boundary-mode slots
            % - Returns policy: initialized index policy
            arguments
                options.expectedNegativeCount (1,1) double = NaN
                options.expectedZeroCount (1,1) double = NaN
                options.indexTolerance (1,1) double {mustBePositive} = 1e-10
                options.validationMode {mustBeTextScalar} = "error"
                options.boundaryModes struct = struct("modeNumber", {}, "indexSign", {})
            end

            self.expectedNegativeCountValue = options.expectedNegativeCount;
            self.expectedZeroCountValue = options.expectedZeroCount;
            self.indexTolerance = options.indexTolerance;
            self.validationMode = string(options.validationMode);
            self.boundaryModes = IMIndexPolicy.sortBoundaryModes(options.boundaryModes);
        end

        function count = expectedNegativeCount(self, ~)
            % Return the expected negative index count.
            %
            % - Topic: Validate index counts
            % - Declaration: count = expectedNegativeCount(policy,context)
            % - Parameter context: solver or EVP context
            % - Returns count: expected negative count
            count = self.expectedNegativeCountValue;
        end

        function count = expectedZeroCount(self, ~)
            % Return the expected zero index count.
            %
            % - Topic: Validate index counts
            % - Declaration: count = expectedZeroCount(policy,context)
            % - Parameter context: solver or EVP context
            % - Returns count: expected zero count
            count = self.expectedZeroCountValue;
        end

        function count = expectedPositiveCount(self, context, nModes)
            % Return the expected positive index count.
            %
            % - Topic: Validate index counts
            % - Declaration: count = expectedPositiveCount(policy,context,nModes)
            % - Parameter context: solver or EVP context
            % - Parameter nModes: number of retained modes
            % - Returns count: expected positive count
            negativeCount = self.expectedNegativeCount(context);
            zeroCount = self.expectedZeroCount(context);
            if isnan(negativeCount) || isnan(zeroCount)
                count = NaN;
            else
                count = max(0, nModes - negativeCount - zeroCount);
            end
        end

        function index = classify(self, eigenvalues, context)
            % Classify eigenvalues and validate the observed index.
            %
            % - Topic: Validate index counts
            % - Declaration: index = classify(policy,eigenvalues,context)
            % - Parameter eigenvalues: eigenvalues to classify
            % - Parameter context: solver or EVP context
            % - Returns index: index summary structure
            lambda = real(eigenvalues(:));
            tolerance = self.indexTolerance*max(1,max(abs(lambda)));
            index.negativeCount = nnz(lambda < -tolerance);
            index.zeroCount = nnz(abs(lambda) <= tolerance);
            index.positiveCount = nnz(lambda > tolerance);
            index.expectedNegativeCount = self.expectedNegativeCount(context);
            index.expectedZeroCount = self.expectedZeroCount(context);
            index.expectedPositiveCount = self.expectedPositiveCount(context, length(lambda));
            index.validationPassed = true;

            if ~isnan(index.expectedNegativeCount) && index.negativeCount ~= index.expectedNegativeCount
                index.validationPassed = false;
            end
            if ~isnan(index.expectedZeroCount) && index.zeroCount ~= index.expectedZeroCount
                index.validationPassed = false;
            end

            if ~index.validationPassed
                message = sprintf("Observed index (%d negative, %d zero) does not match expected index (%g negative, %g zero).", ...
                    index.negativeCount, index.zeroCount, index.expectedNegativeCount, index.expectedZeroCount);
                switch self.validationMode
                    case "error"
                        error("IMIndexPolicy:IndexMismatch", "%s", message);
                    case "warning"
                        warning("IMIndexPolicy:IndexMismatch", "%s", message);
                    case "none"
                    otherwise
                        error("IMIndexPolicy:InvalidValidationMode", ...
                            "Unknown validation mode ""%s"".", self.validationMode);
                end
            end
        end

        function selection = selectModes(self, eigenvalues, nModes, context)
            % Select and label retained modes according to the index policy.
            %
            % The returned `modeNumber` uses `-1` for a surface boundary
            % mode, `-2` for a bottom boundary mode, `0` for a true null
            % mode with $$F_0(z)=1$$ and $$G_0(z)=0$$, and positive labels
            % for interior baroclinic modes.
            %
            % - Topic: Validate index counts
            % - Declaration: selection = selectModes(policy,eigenvalues,nModes,context)
            % - Parameter eigenvalues: candidate eigenvalues
            % - Parameter nModes: number of modes to retain
            % - Parameter context: solver or EVP context
            % - Returns selection: structure with `sortIndex`, `modeNumber`, and `index`
            lambda = real(eigenvalues(:));
            expectedNegativeCount = self.expectedNegativeCount(context);
            expectedZeroCount = self.expectedZeroCount(context);
            if isnan(expectedNegativeCount) || isnan(expectedZeroCount)
                [~, sortIndex] = sortrows([IMIndexPolicy.signWithZero(lambda, self.indexTolerance), abs(lambda), lambda]);
                sortIndex = sortIndex(1:min(nModes, length(sortIndex)));
                selection.sortIndex = sortIndex(:).';
                selection.modeNumber = 1:length(selection.sortIndex);
                selection.index = self.classify(lambda(sortIndex), context);
                return;
            end

            signs = IMIndexPolicy.signWithZero(lambda, self.indexTolerance);
            available = true(size(lambda));
            sortIndex = zeros(0,1);
            modeNumber = zeros(1,0);
            selectedBoundaryNegativeCount = 0;
            selectedBoundaryZeroCount = 0;

            for iBoundaryMode = 1:length(self.boundaryModes)
                if length(sortIndex) >= nModes
                    break;
                end
                boundaryMode = self.boundaryModes(iBoundaryMode);
                candidates = IMIndexPolicy.sortedCandidates(lambda, signs, available, boundaryMode.indexSign);
                if isempty(candidates)
                    candidates = IMIndexPolicy.sortedAvailable(lambda, signs, available);
                end
                if isempty(candidates)
                    continue;
                end
                sortIndex(end+1,1) = candidates(1);
                modeNumber(end+1) = boundaryMode.modeNumber;
                available(candidates(1)) = false;
                selectedBoundaryNegativeCount = selectedBoundaryNegativeCount + double(boundaryMode.indexSign < 0);
                selectedBoundaryZeroCount = selectedBoundaryZeroCount + double(boundaryMode.indexSign == 0);
            end

            remainingNegativeCount = max(0, expectedNegativeCount - selectedBoundaryNegativeCount);
            negativeCandidates = IMIndexPolicy.sortedCandidates(lambda, signs, available, -1);
            selectedNegativeCount = min([remainingNegativeCount, length(negativeCandidates), nModes - length(sortIndex)]);
            if selectedNegativeCount > 0
                selectedNegative = negativeCandidates(1:selectedNegativeCount);
                sortIndex = [sortIndex; selectedNegative(:)];
                modeNumber = [modeNumber, IMIndexPolicy.unusedNegativeLabels(selectedNegativeCount, modeNumber)];
                available(selectedNegative) = false;
            end

            nullCount = max(0, expectedZeroCount - selectedBoundaryZeroCount);
            zeroCandidates = IMIndexPolicy.sortedCandidates(lambda, signs, available, 0);
            selectedZeroCount = min([nullCount, length(zeroCandidates), nModes - length(sortIndex)]);
            if selectedZeroCount > 0
                selectedZero = zeroCandidates(1:selectedZeroCount);
                sortIndex = [sortIndex; selectedZero(:)];
                modeNumber = [modeNumber, zeros(1,selectedZeroCount)];
                available(selectedZero) = false;
            end

            missingZeroCount = min(nullCount - selectedZeroCount, nModes - length(sortIndex));
            if missingZeroCount > 0
                promotedZero = IMIndexPolicy.sortedAvailable(lambda, signs, available);
                promotedZero = promotedZero(1:min(missingZeroCount, length(promotedZero)));
                sortIndex = [sortIndex; promotedZero(:)];
                modeNumber = [modeNumber, zeros(1,length(promotedZero))];
                available(promotedZero) = false;
            end

            positiveCandidates = IMIndexPolicy.sortedCandidates(lambda, signs, available, 1);
            selectedPositiveCount = min(length(positiveCandidates), nModes - length(sortIndex));
            if selectedPositiveCount > 0
                selectedPositive = positiveCandidates(1:selectedPositiveCount);
                sortIndex = [sortIndex; selectedPositive(:)];
                modeNumber = [modeNumber, 1:selectedPositiveCount];
                available(selectedPositive) = false;
            end

            if length(sortIndex) < nModes
                extraIndex = IMIndexPolicy.sortedAvailable(lambda, signs, available);
                extraIndex = extraIndex(1:min(nModes - length(sortIndex), length(extraIndex)));
                nextModeNumber = nnz(modeNumber > 0) + 1;
                sortIndex = [sortIndex; extraIndex(:)];
                modeNumber = [modeNumber, nextModeNumber:(nextModeNumber + length(extraIndex) - 1)];
            end

            selection.sortIndex = sortIndex(:).';
            selection.modeNumber = modeNumber(1:length(selection.sortIndex));
            selection.index = self.classify(lambda(sortIndex), context);
        end
    end

    methods (Static)
        function policy = none(options)
            % Create a policy that records but does not validate index counts.
            %
            % - Topic: Create index policies
            % - Declaration: policy = IMIndexPolicy.none(options)
            % - Parameter options.indexTolerance: tolerance for zero classification
            % - Returns policy: initialized index policy
            arguments
                options.indexTolerance (1,1) double {mustBePositive} = 1e-10
            end

            policy = IMIndexPolicy(indexTolerance=options.indexTolerance, validationMode="none");
        end

        function policy = fixed(options)
            % Create a policy with fixed expected negative and zero counts.
            %
            % - Topic: Create index policies
            % - Declaration: policy = IMIndexPolicy.fixed(options)
            % - Parameter options.expectedNegativeCount: expected negative count
            % - Parameter options.expectedZeroCount: expected zero count
            % - Parameter options.validationMode: `"error"`, `"warning"`, or `"none"`
            % - Returns policy: initialized index policy
            arguments
                options.expectedNegativeCount (1,1) double {mustBeInteger, mustBeNonnegative} = 0
                options.expectedZeroCount (1,1) double {mustBeInteger, mustBeNonnegative} = 0
                options.validationMode {mustBeTextScalar} = "error"
            end

            policy = IMIndexPolicy(expectedNegativeCount=options.expectedNegativeCount, ...
                expectedZeroCount=options.expectedZeroCount, validationMode=options.validationMode);
        end

        function policy = fromBoundarySigns(signs, options)
            % Create an index policy from active-boundary signs.
            %
            % Negative signs add one negative-index direction each.
            %
            % - Topic: Create index policies
            % - Declaration: policy = IMIndexPolicy.fromBoundarySigns(signs,options)
            % - Parameter signs: active-boundary signs
            % - Parameter options.expectedZeroCount: expected zero count
            % - Parameter options.validationMode: `"error"`, `"warning"`, or `"none"`
            % - Returns policy: initialized index policy
            arguments
                signs (:,1) double
                options.expectedZeroCount (1,1) double {mustBeInteger, mustBeNonnegative} = 0
                options.validationMode {mustBeTextScalar} = "error"
            end

            policy = IMIndexPolicy.fixed(expectedNegativeCount=nnz(signs < 0), ...
                expectedZeroCount=options.expectedZeroCount, validationMode=options.validationMode);
        end

        function policy = fromBoundaryConditions(boundaryConditions, options)
            % Create an index policy from boundary-condition index metadata.
            %
            % Placed boundary conditions contribute their negative and zero
            % index counts. Conditions with unknown compatible inner-product
            % terms do not contribute expected counts.
            %
            % - Topic: Create index policies
            % - Declaration: policy = IMIndexPolicy.fromBoundaryConditions(boundaryConditions,options)
            % - Parameter boundaryConditions: boundary-condition array
            % - Parameter options.expectedZeroCount: additional expected zero count
            % - Parameter options.validationMode: `"error"`, `"warning"`, or `"none"`
            % - Returns policy: initialized index policy
            arguments
                boundaryConditions (:,1) IMBoundary
                options.expectedZeroCount (1,1) double {mustBeInteger, mustBeNonnegative} = 0
                options.validationMode {mustBeTextScalar} = "error"
            end

            expectedNegativeCount = 0;
            expectedZeroCount = options.expectedZeroCount;
            boundaryModes = struct("modeNumber", {}, "indexSign", {});
            for iBoundary = 1:length(boundaryConditions)
                expectedNegativeCount = expectedNegativeCount + boundaryConditions(iBoundary).expectedNegativeCount();
                expectedZeroCount = expectedZeroCount + boundaryConditions(iBoundary).expectedZeroCount();
                boundaryModes = [boundaryModes; boundaryConditions(iBoundary).boundaryModeDescriptors()]; %#ok<AGROW>
            end
            policy = IMIndexPolicy(expectedNegativeCount=expectedNegativeCount, ...
                expectedZeroCount=expectedZeroCount, validationMode=options.validationMode, boundaryModes=boundaryModes);
        end
    end

    methods (Static, Access = private)
        function boundaryModes = sortBoundaryModes(boundaryModes)
            if isempty(boundaryModes)
                boundaryModes = struct("modeNumber", {}, "indexSign", {});
                return;
            end
            modeNumbers = [boundaryModes.modeNumber];
            uniqueModeNumbers = unique(modeNumbers, "stable");
            uniqueBoundaryModes = struct("modeNumber", {}, "indexSign", {});
            for iMode = 1:length(uniqueModeNumbers)
                matching = find(modeNumbers == uniqueModeNumbers(iMode));
                signs = [boundaryModes(matching).indexSign];
                if any(signs ~= signs(1))
                    error("IMIndexPolicy:ConflictingBoundaryModeSigns", ...
                        "Duplicate endpoint boundary mode numbers must declare the same eigenvalue sign.");
                end
                uniqueBoundaryModes(end+1,1) = boundaryModes(matching(1)); %#ok<AGROW>
            end
            boundaryModes = uniqueBoundaryModes;
            modeNumbers = [boundaryModes.modeNumber];
            [~, order] = sort(modeNumbers, "descend");
            boundaryModes = boundaryModes(order);
        end

        function candidates = sortedCandidates(lambda, signs, available, signValue)
            switch signValue
                case -1
                    candidates = find(available & signs < 0);
                    [~, order] = sort(abs(lambda(candidates)), "ascend");
                case 0
                    candidates = find(available & signs == 0);
                    [~, order] = sort(abs(lambda(candidates)), "ascend");
                case 1
                    candidates = find(available & signs > 0);
                    [~, order] = sort(lambda(candidates), "ascend");
                otherwise
                    candidates = zeros(0,1);
                    order = zeros(0,1);
            end
            candidates = candidates(order);
        end

        function candidates = sortedAvailable(lambda, signs, available)
            candidates = find(available);
            [~, order] = sortrows([signs(candidates), abs(lambda(candidates)), lambda(candidates)]);
            candidates = candidates(order);
        end

        function labels = unusedNegativeLabels(count, usedLabels)
            labels = zeros(1,count);
            nextLabel = -1;
            for iLabel = 1:count
                while any(usedLabels == nextLabel) || any(labels == nextLabel)
                    nextLabel = nextLabel - 1;
                end
                labels(iLabel) = nextLabel;
            end
        end

        function signs = signWithZero(values, indexTolerance)
            nonzeroValues = abs(values(abs(values) > 0 & isfinite(values)));
            if isempty(nonzeroValues)
                scale = 1;
            else
                scale = max(1,median(nonzeroValues));
            end
            tolerance = indexTolerance*scale;
            signs = ones(size(values));
            signs(values < -tolerance) = -1;
            signs(abs(values) <= tolerance) = 0;
        end
    end
end

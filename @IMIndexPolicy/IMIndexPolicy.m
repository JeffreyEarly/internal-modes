classdef IMIndexPolicy
    % Specify the expected eigenvalue index of an EVP.
    %
    % The index records negative, zero, and positive eigenvalue counts.
    % Negative index directions represent active boundary modes, zero index
    % directions represent barotropic or null modes, and positive index
    % directions represent baroclinic modes.
    %
    % ```matlab
    % boundaryRows = IMBoundaryRow.partialDepthPE(boundarySign="negative");
    % policy = IMIndexPolicy.fromBoundaryRows(boundaryRows);
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
            % - Returns policy: initialized index policy
            arguments
                options.expectedNegativeCount (1,1) double = NaN
                options.expectedZeroCount (1,1) double = NaN
                options.indexTolerance (1,1) double {mustBePositive} = 1e-10
                options.validationMode {mustBeTextScalar} = "error"
            end

            self.expectedNegativeCountValue = options.expectedNegativeCount;
            self.expectedZeroCountValue = options.expectedZeroCount;
            self.indexTolerance = options.indexTolerance;
            self.validationMode = string(options.validationMode);
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
            % The returned `modeIndex` uses negative labels for boundary
            % modes, zero for barotropic or null modes, and positive labels
            % for baroclinic modes.
            %
            % - Topic: Validate index counts
            % - Declaration: selection = selectModes(policy,eigenvalues,nModes,context)
            % - Parameter eigenvalues: candidate eigenvalues
            % - Parameter nModes: number of modes to retain
            % - Parameter context: solver or EVP context
            % - Returns selection: structure with `sortIndex`, `modeIndex`, and `index`
            lambda = real(eigenvalues(:));
            expectedNegativeCount = self.expectedNegativeCount(context);
            expectedZeroCount = self.expectedZeroCount(context);
            if isnan(expectedNegativeCount) || isnan(expectedZeroCount)
                [~, sortIndex] = sortrows([IMIndexPolicy.signWithZero(lambda, self.indexTolerance), abs(lambda), lambda]);
                sortIndex = sortIndex(1:min(nModes, length(sortIndex)));
                selection.sortIndex = sortIndex(:).';
                selection.modeIndex = 1:length(selection.sortIndex);
                selection.index = self.classify(lambda(sortIndex), context);
                return;
            end

            signs = IMIndexPolicy.signWithZero(lambda, self.indexTolerance);
            negativeIndex = find(signs < 0);
            zeroIndex = find(signs == 0);
            positiveIndex = find(signs > 0);

            [~, negativeOrder] = sort(abs(lambda(negativeIndex)), "ascend");
            [~, zeroOrder] = sort(abs(lambda(zeroIndex)), "ascend");
            [~, positiveOrder] = sort(lambda(positiveIndex), "ascend");
            negativeIndex = negativeIndex(negativeOrder);
            zeroIndex = zeroIndex(zeroOrder);
            positiveIndex = positiveIndex(positiveOrder);

            selectedNegativeCount = min(expectedNegativeCount, length(negativeIndex));
            selectedZeroCount = min(expectedZeroCount, length(zeroIndex));
            selectedNegative = negativeIndex(1:selectedNegativeCount);
            selectedZero = zeroIndex(1:selectedZeroCount);
            remainingNegative = negativeIndex((selectedNegativeCount+1):end);
            remainingZero = zeroIndex((selectedZeroCount+1):end);

            sortIndex = [selectedNegative(:); selectedZero(:); positiveIndex(:); remainingZero(:); remainingNegative(:)];
            missingZeroCount = expectedZeroCount - selectedZeroCount;
            if missingZeroCount > 0
                unusedIndex = setdiff((1:length(lambda)).', sortIndex, "stable");
                [~, unusedOrder] = sort(abs(lambda(unusedIndex)), "ascend");
                promotedZero = unusedIndex(unusedOrder(1:min(missingZeroCount, length(unusedOrder))));
                sortIndex = [selectedNegative(:); selectedZero(:); promotedZero(:); positiveIndex(:); remainingZero(:); remainingNegative(:)];
            end
            if length(sortIndex) < nModes
                unusedIndex = setdiff((1:length(lambda)).', sortIndex, "stable");
                sortIndex = [sortIndex(:); unusedIndex(:)];
            end
            sortIndex = sortIndex(1:min(nModes, length(sortIndex)));

            selection.sortIndex = sortIndex(:).';
            selection.modeIndex = self.labelsForSelection(length(sortIndex), expectedNegativeCount, expectedZeroCount);
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

        function policy = fromBoundaryRows(boundaryRows, options)
            % Create an index policy from boundary-row index metadata.
            %
            % Resolved boundary rows contribute their negative and zero index
            % counts. Unresolved boundary rows do not contribute expected
            % counts.
            %
            % - Topic: Create index policies
            % - Declaration: policy = IMIndexPolicy.fromBoundaryRows(boundaryRows,options)
            % - Parameter boundaryRows: boundary-row array
            % - Parameter options.expectedZeroCount: additional expected zero count
            % - Parameter options.validationMode: `"error"`, `"warning"`, or `"none"`
            % - Returns policy: initialized index policy
            arguments
                boundaryRows (:,1) IMBoundaryRow
                options.expectedZeroCount (1,1) double {mustBeInteger, mustBeNonnegative} = 0
                options.validationMode {mustBeTextScalar} = "error"
            end

            expectedNegativeCount = 0;
            expectedZeroCount = options.expectedZeroCount;
            for iBoundary = 1:length(boundaryRows)
                expectedNegativeCount = expectedNegativeCount + boundaryRows(iBoundary).expectedNegativeCount();
                expectedZeroCount = expectedZeroCount + boundaryRows(iBoundary).expectedZeroCount();
            end
            policy = IMIndexPolicy.fixed(expectedNegativeCount=expectedNegativeCount, ...
                expectedZeroCount=expectedZeroCount, validationMode=options.validationMode);
        end
    end

    methods (Access = private)
        function modeIndex = labelsForSelection(~, nModes, expectedNegativeCount, expectedZeroCount)
            selectedNegativeCount = min(expectedNegativeCount, nModes);
            negativeLabels = -selectedNegativeCount:-1;
            remainingCount = nModes - selectedNegativeCount;
            selectedZeroCount = min(expectedZeroCount, remainingCount);
            zeroLabels = zeros(1,selectedZeroCount);
            positiveCount = nModes - selectedNegativeCount - selectedZeroCount;
            modeIndex = [negativeLabels, zeroLabels, 1:positiveCount];
        end
    end

    methods (Static, Access = private)
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

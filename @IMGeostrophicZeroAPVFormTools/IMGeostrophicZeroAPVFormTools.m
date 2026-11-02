classdef (Hidden, Sealed) IMGeostrophicZeroAPVFormTools
    % Share representation-independent zero-APV form and rotation operations.

    methods (Static)
        function [responseMetric,physicalEnergy,surfaceBuoyancy,bottomBuoyancy] = formMatrices(FSurface,FBottom,surfaceResponse,bottomResponse,selectedResponse,endpoints,k,f0,g)
            nEndpoints = numel(endpoints);
            nK = numel(k);
            responseMetric = zeros(nEndpoints,nEndpoints,nK);
            physicalEnergy = zeros(nEndpoints,nEndpoints,nK);
            surfaceBuoyancy = zeros(nEndpoints,nEndpoints,nK);
            bottomBuoyancy = zeros(nEndpoints,nEndpoints,nK);

            for iK = 1:nK
                response = selectedResponse(:,:,iK);
                responseMetric(:,:,iK) = IMGeostrophicZeroAPVFormTools.symmetrize(response.'*response);
                energyScale = f0^2/(2*g*k(iK)^4);
                boundaryForm = -FSurface(:,:,iK).'*surfaceResponse(:,:,iK) + FBottom(:,:,iK).'*bottomResponse(:,:,iK);
                physicalEnergy(:,:,iK) = energyScale*IMGeostrophicZeroAPVFormTools.symmetrize(boundaryForm);
                buoyancyScale = f0^2/(2*g^2*k(iK)^4);
                if ismember("surface",endpoints)
                    surfaceBuoyancy(:,:,iK) = buoyancyScale*IMGeostrophicZeroAPVFormTools.symmetrize(surfaceResponse(:,:,iK).'*surfaceResponse(:,:,iK));
                end
                if ismember("bottom",endpoints)
                    bottomBuoyancy(:,:,iK) = buoyancyScale*IMGeostrophicZeroAPVFormTools.symmetrize(bottomResponse(:,:,iK).'*bottomResponse(:,:,iK));
                end
            end
        end

        function pages = expandMatrixPages(matrix,nEndpoints,nK,name,errorPrefix)
            if ndims(matrix) > 3 || size(matrix,1) ~= nEndpoints || size(matrix,2) ~= nEndpoints
                error(errorPrefix + ":InvalidPencilShape", "%s must have nEndpoints rows and columns.",name);
            end
            if size(matrix,3) == 1
                pages = repmat(matrix,1,1,nK);
            elseif size(matrix,3) == nK
                pages = matrix;
            else
                error(errorPrefix + ":InvalidPencilShape", "%s must have one page or nK pages.",name);
            end
            if ~isreal(pages) || any(~isfinite(pages),"all")
                error(errorPrefix + ":InvalidPencilMatrix", "%s must be finite and real.",name);
            end
            for iK = 1:nK
                page = pages(:,:,iK);
                scale = max([abs(page(:));0]);
                if scale == 0, continue; end
                unitPage = page/scale;
                if norm(unitPage-unitPage.',"fro") > 1e-12*norm(unitPage,"fro")
                    error(errorPrefix + ":NonSymmetricPencil", "%s must be symmetric on every wavenumber page.",name);
                end
                pages(:,:,iK) = IMGeostrophicZeroAPVFormTools.symmetrize(page);
            end
        end

        function [vectors,eigenvalues,rightSignatures,residuals] = solveSymmetricPencil(leftMatrix,rightMatrix,nEndpoints,nK,ordering,errorPrefix)
            leftMatrix = IMGeostrophicZeroAPVFormTools.expandMatrixPages(leftMatrix,nEndpoints,nK,"leftMatrix",errorPrefix);
            rightMatrix = IMGeostrophicZeroAPVFormTools.expandMatrixPages(rightMatrix,nEndpoints,nK,"rightMatrix",errorPrefix);
            vectors = zeros(nEndpoints,nEndpoints,nK);
            eigenvalues = zeros(nEndpoints,nK);
            rightSignatures = zeros(nEndpoints,nK);
            relativeResidual = zeros(nEndpoints,nK);
            rightOrthogonalityError = zeros(1,nK);
            rotationConditionNumber = zeros(1,nK);
            qualificationTolerance = 1e-10;
            metricRelativeTolerance = 100*eps;

            for iK = 1:nK
                if nEndpoints == 0
                    rotationConditionNumber(iK) = 1;
                    continue
                end
                left = leftMatrix(:,:,iK);
                right = rightMatrix(:,:,iK);
                leftScale = max(abs(left),[],"all");
                rightScale = max(abs(right),[],"all");
                if rightScale == 0
                    error(errorPrefix + ":NullRightMetric", "The right metric must be nonsingular on the retained coordinate space.");
                end
                if leftScale == 0, leftScale = 1; end
                unitLeft = left/leftScale;
                unitRight = right/rightScale;
                if any(left ~= 0 & unitLeft == 0,"all") || any(right ~= 0 & unitRight == 0,"all")
                    error(errorPrefix + ":InvalidPencilSpectrum", "Scaling the pencil must not discard nonzero matrix entries outside floating-point range.");
                end
                metricEigenvalues = eig(unitRight);
                if min(abs(metricEigenvalues)) <= metricRelativeTolerance*max(abs(metricEigenvalues))
                    error(errorPrefix + ":NullRightMetric", "The right metric is numerically singular relative to its own scale.");
                end
                [C,Lambda] = eig(unitLeft,unitRight);
                lambda = diag(Lambda);
                imaginaryScale = max(abs(lambda));
                if any(~isfinite(lambda)) || any(abs(imag(lambda)) > qualificationTolerance*imaginaryScale)
                    error(errorPrefix + ":InvalidPencilSpectrum", "The matrix pencil must have a complete finite real spectrum on every wavenumber page.");
                end
                C = real(C);
                lambda = real(lambda);
                % Generalized eigensolvers may return a nonorthogonal or
                % even individually null basis inside a repeated eigenspace.
                % Normalize its signed restricted metric as a whole.
                unused = true(nEndpoints,1);
                for first = 1:nEndpoints
                    if ~unused(first), continue; end
                    clusterScale = max(abs(lambda),abs(lambda(first)));
                    cluster = find(unused & abs(lambda-lambda(first)) <= 1000*eps*clusterScale);
                    if isempty(cluster), cluster = first; end
                    unused(cluster) = false;
                    block = C(:,cluster);
                    singularValues = svd(block);
                    if min(singularValues) <= metricRelativeTolerance*max(singularValues)
                        error(errorPrefix + ":SingularCoordinateRotation", "The pencil eigenvectors do not span the retained coordinate space.");
                    end
                    [frame,~] = qr(block,0);
                    restricted = IMGeostrophicZeroAPVFormTools.symmetrize(frame.'*unitRight*frame);
                    [directions,metricValues] = eig(restricted,"vector");
                    if min(abs(metricValues)) <= metricRelativeTolerance*norm(unitRight,2)
                        error(errorPrefix + ":NullRightMetric", "A pencil eigenspace has a numerically null right-metric direction.");
                    end
                    [metricValues,metricOrder] = sort(metricValues,"descend");
                    C(:,cluster) = (frame*directions(:,metricOrder))./sqrt(abs(metricValues.'));
                end

                switch string(ordering)
                    case "nonzeroFirst"
                        zeroTolerance = 1e-11*max(abs(lambda));
                        isZero = abs(lambda) <= zeroTolerance;
                        [~,sortIndex] = sortrows([isZero(:),-abs(lambda(:)),-lambda(:)],[1 2 3]);
                    otherwise
                        [~,sortIndex] = sort(lambda,"descend");
                end
                lambda = lambda(sortIndex);
                C = C(:,sortIndex);

                for iMode = 1:nEndpoints
                    % Reorthogonalize distinct, close roots as well. Their
                    % exact right-metric cross products vanish; correcting
                    % roundoff is accepted only if the pencil still passes.
                    for pass = 1:2
                        previous = C(:,1:iMode-1);
                        crossProducts = previous.'*unitRight*C(:,iMode);
                        C(:,iMode) = C(:,iMode)-previous*(rightSignatures(1:iMode-1,iK).*crossProducts);
                    end
                    rightNorm = C(:,iMode).'*unitRight*C(:,iMode);
                    if abs(rightNorm) <= metricRelativeTolerance*norm(unitRight,2)*norm(C(:,iMode))^2
                        error(errorPrefix + ":NullRightMetric", "A pencil eigenvector has a numerically zero right-metric norm.");
                    end
                    C(:,iMode) = C(:,iMode)/sqrt(abs(rightNorm));
                    [~,referenceIndex] = max(abs(C(:,iMode)));
                    if C(referenceIndex,iMode) < 0
                        C(:,iMode) = -C(:,iMode);
                    end
                    rightSignatures(iMode,iK) = sign(C(:,iMode).'*unitRight*C(:,iMode));
                    residual = unitLeft*C(:,iMode)-lambda(iMode)*unitRight*C(:,iMode);
                    residualScale = (norm(unitLeft,"fro")+abs(lambda(iMode))*norm(unitRight,"fro"))*norm(C(:,iMode));
                    if residualScale == 0
                        relativeResidual(iMode,iK) = norm(residual);
                    else
                        relativeResidual(iMode,iK) = norm(residual)/residualScale;
                    end
                end

                targetRightGram = diag(rightSignatures(:,iK));
                rightOrthogonalityError(iK) = norm(C.'*unitRight*C-targetRightGram,"fro")/max(1,norm(targetRightGram,"fro"));
                if max(relativeResidual(:,iK)) > qualificationTolerance || rightOrthogonalityError(iK) > qualificationTolerance
                    error(errorPrefix + ":UnqualifiedPencil", "The rotation must satisfy both the pencil residual and signed right-metric orthogonality checks.");
                end
                rotationConditionNumber(iK) = cond(C);
                if rcond(C) <= sqrt(eps)
                    error(errorPrefix + ":SingularCoordinateRotation", "The coordinate rotation is singular or ill-conditioned on page %d.",iK);
                end
                physicalVectors = C/sqrt(rightScale);
                physicalEigenvalues = IMGeostrophicZeroAPVFormTools.restoreEigenvalueScale(lambda,leftScale,rightScale);
                if any(~isfinite(physicalVectors),"all") || any(~isfinite(physicalEigenvalues)) || any(lambda ~= 0 & physicalEigenvalues == 0)
                    error(errorPrefix + ":InvalidPencilSpectrum", "The normalized rotation and nonzero eigenvalues must be representable in the supplied physical units.");
                end
                vectors(:,:,iK) = physicalVectors;
                eigenvalues(:,iK) = physicalEigenvalues;
            end
            residuals = struct(relativePencilResidual=relativeResidual,rightOrthogonalityError=rightOrthogonalityError,rotationConditionNumber=rotationConditionNumber);
        end

        function tf = pagesEqual(left,right)
            tf = isequal(size(left),size(right));
            if ~tf, return; end
            scale = max([abs(left(:));abs(right(:));0]);
            if scale == 0, return; end
            left = left/scale; right = right/scale;
            tf = norm(left(:)-right(:)) <= 1e-12*max(norm(left(:)),norm(right(:)));
        end

        function matrix = symmetrize(matrix)
            if isequal(matrix,matrix.'), return; end
            originalDiagonal = diag(matrix);
            matrix = 0.5*matrix+0.5*matrix.';
            matrix(1:size(matrix,1)+1:end) = originalDiagonal;
        end

        function pages = symmetrizePages(pages)
            for iPage = 1:size(pages,3)
                pages(:,:,iPage) = IMGeostrophicZeroAPVFormTools.symmetrize(pages(:,:,iPage));
            end
        end

        function validateRotationPages(rotation,nEndpoints,nK,errorPrefix)
            if ~isequal(size(rotation),[nEndpoints nEndpoints nK]) && ~(nK == 1 && isequal(size(rotation),[nEndpoints nEndpoints]))
                error(errorPrefix + ":InvalidRotationShape", "A rotation must have nEndpoints rows and columns and one page per wavenumber.");
            end
            for iK = 1:nK
                if nEndpoints == 0, continue; end
                if rcond(rotation(:,:,iK)) <= sqrt(eps)
                    error(errorPrefix + ":SingularCoordinateRotation", "The coordinate rotation is singular or ill-conditioned on page %d.",iK);
                end
            end
        end
    end

    methods (Static, Access = private)
        function values = restoreEigenvalueScale(values,leftScale,rightScale)
            [leftFraction,leftExponent] = log2(leftScale);
            [rightFraction,rightExponent] = log2(rightScale);
            [valueFraction,valueExponent] = log2(abs(values));
            values = sign(values).*pow2(valueFraction*(leftFraction/rightFraction),valueExponent+leftExponent-rightExponent);
        end
    end
end

classdef IMProductAssessment
    % Report sampled product errors separately from acceptance policy.
    %
    % Each measurement is a maximum over individually examined products.
    % Neither a passing sparse plan nor a bounded all-products control proves
    % arbitrary-superposition accuracy or complete physical qualification.
    %
    % ```matlab
    % result = plan.assess(grids);
    % decision = result.applyPolicy(quadraticTolerance=0.1,referenceTolerance=1e-4);
    % ```
    %
    % - Topic: Inspect product evidence
    % - Topic: Apply product policies
    % - Declaration: classdef IMProductAssessment
    properties (SetAccess = private)
        % Immutable reservation and original scientific inventory.
        % - Topic: Inspect product evidence
        plan
        % Exact grids and prepared numerical output/reference identities.
        % - Topic: Inspect product evidence
        assessmentIdentity
        % One error/status/limiting-identity row per requested count.
        % - Topic: Inspect product evidence
        measurements
        % Individual products, signs, labels, zeros, and reference evidence.
        % - Topic: Inspect product evidence
        evidence
        % Examined and omitted inventory families and interactions.
        % - Topic: Inspect product evidence
        coverage
        % Planning, preparation, evaluation, assessment, and array memory costs.
        % - Topic: Inspect product evidence
        costs
    end
    methods
        function self = IMProductAssessment(plan,evidence,costs,identity)
            % Aggregate complete executed evidence without applying tolerances.
            % - Topic: Inspect product evidence
            arguments (Input)
                plan (1,1) IMProductPlan
                evidence (1,:) cell
                costs (1,1) struct
                identity (1,1) struct
            end
            counts=plan.prefixCounts; n=numel(counts);
            if numel(evidence)~=numel(plan.rows)
                error("IMProductAssessment:InvalidEvidence","Evidence must include every reserved product-family row.");
            end
            value=zeros(n,1); referenceError=zeros(n,1); quadratureError=zeros(n,1); independentSolveError=zeros(n,1);
            status=repmat("measured",n,1); referenceStatus=repmat("qualified",n,1);
            examinedCount=zeros(n,1); zeroCount=zeros(n,1); limitingInputI=strings(n,1); limitingInputJ=strings(n,1);
            limitingInteraction=strings(n,1); limitingChannel=strings(n,1); limitingOutput=strings(n,1);
            for q=1:n
                hasValue=false;
                for j=1:numel(evidence)
                    e=evidence{j}; selected=e.firstCounts<=counts(q); indices=find(selected);
                    if isempty(indices), continue; end
                    examinedCount(q)=examinedCount(q)+nnz(selected); zeroCount(q)=zeroCount(q)+nnz(e.isZero(selected));
                    [worst,i]=max(e.error(q,selected));
                    if ~hasValue || worst>value(q)
                        value(q)=worst; index=indices(i); limitingInputI(q)=e.inputLabels(1,index); limitingInputJ(q)=e.inputLabels(2,index);
                        limitingInteraction(q)=e.interactionId; limitingChannel(q)=e.channel; limitingOutput(q)=e.outputId; hasValue=true;
                    end
                    referenceError(q)=max(referenceError(q),max(e.referenceError(q,selected)));
                    quadratureError(q)=max(quadratureError(q),max(e.quadratureError(q,selected)));
                    independentSolveError(q)=max(independentSolveError(q),max(e.independentSolveError(q,selected)));
                    if e.referenceStatus~="qualified", referenceStatus(q)="inconclusive"; status(q)="inconclusive"; end
                end
                if ~hasValue, value(q)=NaN; status(q)="inconclusive"; end
            end
            columnCount=counts(:); quantity=repmat("quadraticAliasing",n,1);
            self.measurements=table(columnCount,quantity,value,status,referenceStatus,referenceError,quadratureError,independentSolveError,examinedCount,zeroCount,limitingInputI,limitingInputJ,limitingInteraction,limitingChannel,limitingOutput);
            products=plan.inventory.products; families=strings(numel(evidence),3); productCounts=zeros(numel(evidence),1);
            for j=1:numel(evidence)
                e=evidence{j}; families(j,:)=[e.inputFamilies e.outputFamily]; productCounts(j)=numel(e.isZero);
            end
            [uniqueFamilies,~,group]=unique(families,"rows","stable");
            examinedFamilies=table(uniqueFamilies(:,1),uniqueFamilies(:,2),uniqueFamilies(:,3),accumarray(group,productCounts),VariableNames=["inputA","inputB","output","productCount"]);
            availableDeclaredProducts=zeros(n,1); selectedFamilyAvailableProducts=zeros(n,1);
            for q=1:n
                for row=1:height(products)
                    a=plan.inventory.factors{products.factorA(row)}; b=plan.inventory.factors{products.factorB(row)};
                    na=numel(a.labels); nb=numel(b.labels);
                    if a.countRole=="retained", na=nnz(a.ordinals<=counts(q)); end
                    if b.countRole=="retained", nb=nnz(b.ordinals<=counts(q)); end
                    availableDeclaredProducts(q)=availableDeclaredProducts(q)+na*nb;
                    if ismember(row,plan.rows), selectedFamilyAvailableProducts(q)=selectedFamilyAvailableProducts(q)+na*nb; end
                end
            end
            prefixCoverage=table(counts(:),availableDeclaredProducts,selectedFamilyAvailableProducts,examinedCount,selectedFamilyAvailableProducts-examinedCount,availableDeclaredProducts-examinedCount,VariableNames=["retainedCount","availableDeclaredProducts","selectedFamilyAvailableProducts","examinedProducts","omittedWithinSelectedFamilies","omittedProducts"]);
            self.coverage=struct(selection=plan.selection,selectedRows=plan.rows,omittedRows=setdiff(1:height(products),plan.rows),selectedInteractions=unique(products.interactionId(plan.rows),"stable"),omittedInteractions=setdiff(unique(products.interactionId,"stable"),unique(products.interactionId(plan.rows),"stable"),"stable"),examinedFamilies=examinedFamilies,prefixCoverage=prefixCoverage,selectedFamilyRowCount=numel(plan.rows),totalFamilyRowCount=height(products),exhaustiveGuarantee=false,superpositionGuarantee=false,scope="Only supplied valid interactions and explicit physical source recipes; no model qualification");
            self.plan=plan; self.evidence=evidence; self.costs=costs; self.assessmentIdentity=identity;
        end

        function decision = applyPolicy(self,options)
            % Apply explicit product/reference tolerances without reducing counts.
            %
            % Reference failures remain inconclusive. Qualified product errors
            % exceeding tolerance reject. Strict requests refer to an examined
            % count and retain all scientific columns supplied by the caller.
            %
            % - Topic: Apply product policies
            % - Parameter options.quadraticTolerance: nonnegative product-error tolerance
            % - Parameter options.referenceTolerance: maximum qualified reference discrepancy
            % - Parameter options.requestedCount: optional strict examined retained count
            arguments (Input)
                self (1,1) IMProductAssessment
                options.quadraticTolerance (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
                options.referenceTolerance (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
                options.requestedCount (1,:) double {mustBeInteger,mustBePositive} = self.plan.prefixCounts(end)
            end
            if numel(options.requestedCount)~=1 || ~ismember(options.requestedCount,self.plan.prefixCounts)
                error("IMProductAssessment:InvalidRequestedCount","A strict requestedCount must be one explicitly examined count.");
            end
            m=self.measurements; statuses=repmat("accepted",height(m),1);
            for q=1:height(m)
                inconclusive=false; trustedFailure=false; examined=false;
                for j=1:numel(self.evidence)
                    e=self.evidence{j}; selected=e.firstCounts<=m.columnCount(q);
                    if ~any(selected), continue; end
                    examined=true;
                    qualified=e.referenceStatus=="qualified" & e.referenceError(q,selected)<=options.referenceTolerance;
                    errors=e.error(q,selected);
                    trustedFailure=trustedFailure || any(qualified & errors>options.quadraticTolerance);
                    inconclusive=inconclusive || any(~qualified | isnan(errors));
                end
                if inconclusive || ~examined, statuses(q)="inconclusive"; end
                if trustedFailure, statuses(q)="rejected"; end
            end
            consecutive=cumprod(double(statuses=="accepted"))>0;
            selected=find(m.columnCount<=options.requestedCount);
            status="accepted";
            if any(statuses(selected)=="inconclusive"), status="inconclusive"; end
            if any(statuses(selected)=="rejected"), status="rejected"; end
            passing=m.columnCount(consecutive & m.columnCount<=options.requestedCount);
            largest=0; if ~isempty(passing), largest=passing(end); end
            decision=struct(status=status,requestedCount=options.requestedCount,requestedCountAccepted=status=="accepted",largestExaminedAcceptedPrefix=largest,prefixDecisions=table(m.columnCount,statuses,consecutive,VariableNames=["columnCount","status","consecutiveAccepted"]),tolerances=options,countDescription="Largest count passing the explicitly sampled products; scientific columns remain unchanged");
        end

        function comparison = compareCoverage(self,control,options)
            % Detect failures missed by this plan using a separate bounded control.
            %
            % A control must use the same declared inventory and counts, with
            % qualified references. This reports observed missed products;
            % no failure can be detected outside both examined inventories.
            %
            % - Topic: Inspect product evidence
            % - Parameter control: independently selected or all-products assessment
            % - Parameter options.quadraticTolerance: explicit failure threshold
            % - Parameter options.referenceTolerance: required reference qualification
            arguments (Input)
                self (1,1) IMProductAssessment
                control (1,1) IMProductAssessment
                options.quadraticTolerance (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
                options.referenceTolerance (1,1) double {mustBeReal,mustBeFinite,mustBeNonnegative}
            end
            if ~isequal(self.plan.prefixCounts,control.plan.prefixCounts) || ~isequal(self.plan.inventory.products,control.plan.inventory.products) || ~isequaln(self.plan.inventory.factors,control.plan.inventory.factors) || ~isequaln(self.plan.inventory.outputs,control.plan.inventory.outputs) || ~isequaln(self.assessmentIdentity.grids,control.assessmentIdentity.grids) || self.costs.minimumReciprocalCondition~=control.costs.minimumReciprocalCondition
                error("IMProductAssessment:IncompatibleControl","Coverage comparison requires the same scientific inventory, grids, and examined counts.");
            end
            shared=intersect(self.plan.inventory.products.output(self.plan.rows),control.plan.inventory.products.output(control.plan.rows));
            for index=shared.'
                left=self.assessmentIdentity.outputs{index}; right=control.assessmentIdentity.outputs{index};
                if ~isequaln(left.projection,right.projection) || ~isequaln(left.references,right.references)
                    error("IMProductAssessment:IncompatibleControl","Coverage comparison requires identical prepared projections and reference recipes on shared outputs.");
                end
            end
            sampledDecision=self.applyPolicy(quadraticTolerance=options.quadraticTolerance,referenceTolerance=options.referenceTolerance);
            controlDecision=control.applyPolicy(quadraticTolerance=options.quadraticTolerance,referenceTolerance=options.referenceTolerance);
            worst=0; limiting=struct(); missed=0;
            for j=1:numel(control.evidence)
                e=control.evidence{j}; own=find(self.plan.rows==e.inventoryRow,1);
                ownKeys=strings(1,0); ownFirst=zeros(1,0);
                if ~isempty(own), ownKeys=pairKeys(self.evidence{own}); ownFirst=self.evidence{own}.firstCounts; end
                [present,location]=ismember(pairKeys(e),ownKeys);
                for k=1:numel(e.isZero)
                    checked=e.firstCounts(k)<=self.plan.prefixCounts;
                    if present(k), checked=checked & self.plan.prefixCounts<ownFirst(location(k)); end
                    q=find(checked); if isempty(q), continue; end
                    [candidate,i]=max(e.error(q,k));
                    if candidate>options.quadraticTolerance, missed=missed+1; end
                    if candidate>worst
                        worst=candidate; limiting=struct(interactionId=e.interactionId,channel=e.channel,inputLabels=e.inputLabels(:,k),frequencySigns=e.frequencySigns(:,k),outputId=e.outputId,columnCount=self.plan.prefixCounts(q(i)));
                    end
                end
            end
            qualified=all(control.measurements.referenceStatus=="qualified" & control.measurements.referenceError<=options.referenceTolerance) && all(self.measurements.referenceStatus=="qualified" & self.measurements.referenceError<=options.referenceTolerance);
            comparison=struct(referenceQualified=qualified,observedMissedFailures=missed,worstMissedError=worst,limitingMissedProduct=limiting,falseAcceptance=qualified && sampledDecision.requestedCountAccepted && ~controlDecision.requestedCountAccepted,sampledDecision=sampledDecision,controlDecision=controlDecision,scope="Only failures exposed by this separate supplied control; unseen products remain unqualified");
            if ~qualified, comparison.falseAcceptance=NaN; comparison.observedMissedFailures=NaN; end
        end
    end
end

function keys = pairKeys(evidence)
keys=strings(1,size(evidence.inputLabels,2));
for j=1:numel(keys)
    keys(j)=jsonencode(struct(labels=evidence.inputLabels(:,j),signs=evidence.frequencySigns(:,j)));
end
end

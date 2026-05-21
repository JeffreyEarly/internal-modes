classdef IMOperator
    % Represent a physical-coordinate linear differential operator.
    %
    % `IMOperator` stores structured terms of the form
    % $$a_p(z)\partial_z^p$$. Coordinate-aware solvers pull these terms
    % back to their native coordinate when assembling EVP matrices.
    % Coefficient handles receive the EVP-built framework context `ctx`,
    % including `ctx.N2(z)`, `ctx.dzLogN2(z)`, `ctx.g`, `ctx.f0`,
    % `ctx.zDomain`, and `ctx.coordinateKind`.
    %
    % ```matlab
    % op = IMOperator.strong().plus(coefficient=@(z,ctx) ctx.N2(z), derivativeOrder=0);
    % ```
    %
    % - Topic: Create operators
    % - Topic: Inspect operators
    % - Topic: Assemble operators
    % - Declaration: classdef IMOperator

    properties (SetAccess = private)
        % Operator representation.
        %
        % The initial implementation supports `"strong"` operators.
        %
        % - Topic: Create operators
        form

        % Structured operator terms.
        %
        % Each term has a coefficient and a physical derivative order.
        %
        % - Topic: Inspect operators
        terms
    end

    methods
        function self = IMOperator(options)
            % Create an empty physical-coordinate operator.
            %
            % - Topic: Create operators
            % - Declaration: op = IMOperator(options)
            % - Parameter options.form: operator representation
            % - Returns op: initialized operator
            arguments
                options.form {mustBeTextScalar} = "strong"
            end

            self.form = string(options.form);
            self.terms = struct("coefficient", {}, "derivativeOrder", {});
        end

        function self = plus(self, options)
            % Add a structured term to the operator.
            %
            % - Topic: Create operators
            % - Declaration: op = op.plus(options)
            % - Parameter options.coefficient: scalar, vector, or function handle coefficient
            % - Parameter options.derivativeOrder: physical derivative order
            % - Returns op: updated operator
            arguments
                self IMOperator
                options.coefficient = 1
                options.derivativeOrder (1,1) double {mustBeInteger, mustBeNonnegative} = 0
            end

            term.coefficient = options.coefficient;
            term.derivativeOrder = options.derivativeOrder;
            self.terms(end+1) = term;
        end

        function M = matrix(self, solver, options)
            % Assemble the operator on a solver's native basis.
            %
            % - Topic: Assemble operators
            % - Declaration: M = matrix(op,solver)
            % - Parameter solver: coordinate-aware internal-mode solver
            % - Parameter options.context: framework coefficient context
            % - Returns M: assembled matrix
            arguments
                self IMOperator
                solver
                options.context struct = struct()
            end

            if self.form ~= "strong"
                error("IMOperator:UnsupportedForm", ...
                    "Only strong-form operators are currently supported.");
            end

            n = solver.nEVP;
            M = zeros(n,n);
            context = IMOperator.resolveContext(solver, options.context);
            z = solver.zNative;
            for iTerm = 1:length(self.terms)
                coefficient = IMOperator.evaluateCoefficient(self.terms(iTerm).coefficient, z, context);
                D = solver.physicalDerivativeMatrix(self.terms(iTerm).derivativeOrder);
                M = M + diag(coefficient(:))*D;
            end
        end

        function row = boundaryRow(self, solver, location, options)
            % Assemble a boundary functional row.
            %
            % - Topic: Assemble operators
            % - Declaration: row = boundaryRow(op,solver,location,options)
            % - Parameter solver: coordinate-aware internal-mode solver
            % - Parameter location: boundary location, `"surface"` or `"bottom"`
            % - Parameter options.context: framework coefficient context
            % - Returns row: assembled boundary row
            arguments
                self IMOperator
                solver
                location {mustBeTextScalar}
                options.context struct = struct()
            end

            if self.form ~= "strong"
                error("IMOperator:UnsupportedForm", ...
                    "Only strong-form boundary functionals are currently supported.");
            end

            index = solver.boundaryIndex(location);
            z = solver.zNative(index);
            context = IMOperator.resolveContext(solver, options.context);
            row = zeros(1,solver.nEVP);
            for iTerm = 1:length(self.terms)
                coefficient = IMOperator.evaluateCoefficient(self.terms(iTerm).coefficient, z, context);
                D = solver.physicalDerivativeMatrix(self.terms(iTerm).derivativeOrder);
                row = row + coefficient*D(index,:);
            end
        end

        function values = evaluate(self, solver, nativeModes, z, options)
            % Evaluate an operator applied to native mode columns.
            %
            % - Topic: Assemble operators
            % - Developer: true
            % - Declaration: values = evaluate(op,solver,nativeModes,z,options)
            % - Parameter solver: coordinate-aware internal-mode solver
            % - Parameter nativeModes: native mode columns
            % - Parameter z: physical-coordinate evaluation points
            % - Parameter options.context: framework coefficient context
            % - Returns values: evaluated operator values
            arguments
                self IMOperator
                solver
                nativeModes double
                z (:,1) double
                options.context struct = struct()
            end

            if self.form ~= "strong"
                error("IMOperator:UnsupportedForm", ...
                    "Only strong-form operator evaluation is currently supported.");
            end

            values = zeros(length(z), size(nativeModes,2));
            context = IMOperator.resolveContext(solver, options.context);
            for iTerm = 1:length(self.terms)
                coefficient = IMOperator.evaluateCoefficient(self.terms(iTerm).coefficient, z, context);
                derivativeValues = solver.evaluatePhysicalDerivative(nativeModes, z, self.terms(iTerm).derivativeOrder);
                values = values + coefficient(:).*derivativeValues;
            end
        end
    end

    methods (Static)
        function op = strong()
            % Create a strong-form physical-coordinate operator.
            %
            % - Topic: Create operators
            % - Declaration: op = IMOperator.strong()
            % - Returns op: empty strong-form operator
            op = IMOperator(form="strong");
        end

        function values = evaluateCoefficient(coefficient, z, context)
            % Evaluate an operator coefficient on a grid.
            %
            % - Topic: Assemble operators
            % - Developer: true
            % - Declaration: values = IMOperator.evaluateCoefficient(coefficient,z,context)
            % - Parameter coefficient: scalar, vector, or function handle coefficient
            % - Parameter z: physical-coordinate evaluation points
            % - Parameter context: framework coefficient context
            % - Returns values: coefficient values matching `z`
            if isa(coefficient, "function_handle")
                try
                    values = coefficient(z, context);
                catch
                    values = coefficient(z);
                end
            else
                values = coefficient;
            end

            if isscalar(values)
                values = values*ones(size(z));
            end
            values = values(:);
        end

        function context = resolveContext(solver, context)
            if isempty(fieldnames(context))
                context = solver.context();
            end
        end
    end
end

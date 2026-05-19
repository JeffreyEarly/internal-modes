classdef IMBoundaryCondition
    % Represent passive and active physical boundary conditions.
    %
    % Passive boundary conditions replace rows in an EVP matrix pair. Active
    % boundary conditions also carry index metadata for indefinite
    % active-boundary quadratic forms.
    %
    % ```matlab
    % bc = IMBoundaryCondition.dirichlet(location="surface", component="G");
    % ```
    %
    % - Topic: Create boundary conditions
    % - Topic: Apply boundary conditions
    % - Topic: Inspect index metadata
    % - Declaration: classdef IMBoundaryCondition

    properties
        % Boundary-condition role.
        %
        % Values are `"passive"` or `"active"`.
        %
        % - Topic: Inspect index metadata
        role = "passive"

        % Physical boundary location.
        %
        % Supported values are `"surface"` and `"bottom"`.
        %
        % - Topic: Apply boundary conditions
        location = "surface"

        % Component name to which the condition applies.
        %
        % - Topic: Apply boundary conditions
        component = "G"

        % Left-side boundary functional.
        %
        % - Topic: Apply boundary conditions
        operator = IMOperator.strong().plus()

        % Right-side eigenvalue boundary functional.
        %
        % - Topic: Apply boundary conditions
        eigenvalueOperator = IMOperator.strong()

        % Boundary value for passive constraints.
        %
        % Homogeneous EVPs currently require this to be zero.
        %
        % - Topic: Apply boundary conditions
        value = 0

        % Sign associated with an active boundary contribution.
        %
        % Negative signs contribute to the negative index.
        %
        % - Topic: Inspect index metadata
        signatureSign = 0

        % Number of independent active boundary directions.
        %
        % - Topic: Inspect index metadata
        rank = 0
    end

    methods
        function self = IMBoundaryCondition(options)
            % Create a physical boundary condition.
            %
            % - Topic: Create boundary conditions
            % - Declaration: bc = IMBoundaryCondition(options)
            % - Parameter options.role: `"passive"` or `"active"`
            % - Parameter options.location: boundary location
            % - Parameter options.component: component name
            % - Parameter options.operator: left-side boundary functional
            % - Parameter options.eigenvalueOperator: right-side eigenvalue functional
            % - Parameter options.value: passive boundary value
            % - Parameter options.signatureSign: active-boundary sign
            % - Parameter options.rank: active-boundary rank
            % - Returns bc: initialized boundary condition
            arguments
                options.role {mustBeTextScalar} = "passive"
                options.location {mustBeTextScalar} = "surface"
                options.component {mustBeTextScalar} = "G"
                options.operator IMOperator = IMOperator.strong().plus()
                options.eigenvalueOperator IMOperator = IMOperator.strong()
                options.value (1,1) double = 0
                options.signatureSign (1,1) double = 0
                options.rank (1,1) double {mustBeInteger, mustBeNonnegative} = 0
            end

            self.role = string(options.role);
            self.location = string(options.location);
            self.component = string(options.component);
            self.operator = options.operator;
            self.eigenvalueOperator = options.eigenvalueOperator;
            self.value = options.value;
            self.signatureSign = options.signatureSign;
            self.rank = options.rank;
        end

        function [A, B] = apply(self, A, B, solver)
            % Apply the boundary condition to a matrix pair.
            %
            % - Topic: Apply boundary conditions
            % - Declaration: [A,B] = apply(bc,A,B,solver)
            % - Parameter A: left EVP matrix
            % - Parameter B: right EVP matrix
            % - Parameter solver: coordinate-aware internal-mode solver
            % - Returns A: boundary-conditioned left matrix
            % - Returns B: boundary-conditioned right matrix
            if self.role == "active"
                return;
            end
            if self.value ~= 0
                error("IMBoundaryCondition:NonhomogeneousEVP", ...
                    "The v2 EVP assembler currently supports only homogeneous boundary values.");
            end

            index = solver.boundaryIndex(self.location);
            A(index,:) = self.operator.boundaryRow(solver, self.location);
            B(index,:) = self.eigenvalueOperator.boundaryRow(solver, self.location);
        end
    end

    methods (Static)
        function bc = dirichlet(options)
            % Create a homogeneous Dirichlet boundary condition.
            %
            % - Topic: Create boundary conditions
            % - Declaration: bc = IMBoundaryCondition.dirichlet(options)
            % - Parameter options.location: boundary location
            % - Parameter options.component: component name
            % - Returns bc: initialized boundary condition
            arguments
                options.location {mustBeTextScalar}
                options.component {mustBeTextScalar} = "G"
            end

            op = IMOperator.strong().plus(derivativeOrder=0);
            bc = IMBoundaryCondition(location=options.location, component=options.component, operator=op);
        end

        function bc = neumann(options)
            % Create a homogeneous Neumann boundary condition.
            %
            % - Topic: Create boundary conditions
            % - Declaration: bc = IMBoundaryCondition.neumann(options)
            % - Parameter options.location: boundary location
            % - Parameter options.component: component name
            % - Returns bc: initialized boundary condition
            arguments
                options.location {mustBeTextScalar}
                options.component {mustBeTextScalar} = "G"
            end

            op = IMOperator.strong().plus(derivativeOrder=1);
            bc = IMBoundaryCondition(location=options.location, component=options.component, operator=op);
        end

        function bc = robin(options)
            % Create a homogeneous Robin boundary condition.
            %
            % - Topic: Create boundary conditions
            % - Declaration: bc = IMBoundaryCondition.robin(options)
            % - Parameter options.location: boundary location
            % - Parameter options.component: component name
            % - Parameter options.operator: left-side boundary functional
            % - Parameter options.eigenvalueOperator: right-side eigenvalue functional
            % - Returns bc: initialized boundary condition
            arguments
                options.location {mustBeTextScalar}
                options.component {mustBeTextScalar} = "G"
                options.operator IMOperator
                options.eigenvalueOperator IMOperator = IMOperator.strong()
            end

            bc = IMBoundaryCondition(location=options.location, component=options.component, ...
                operator=options.operator, eigenvalueOperator=options.eigenvalueOperator);
        end

        function bc = active(options)
            % Create an active boundary contribution for index accounting.
            %
            % - Topic: Create boundary conditions
            % - Declaration: bc = IMBoundaryCondition.active(options)
            % - Parameter options.location: boundary location
            % - Parameter options.component: component name
            % - Parameter options.signatureSign: active-boundary sign
            % - Parameter options.rank: number of active directions
            % - Returns bc: initialized active boundary condition
            arguments
                options.location {mustBeTextScalar}
                options.component {mustBeTextScalar} = "G"
                options.signatureSign (1,1) double {mustBeMember(options.signatureSign, [-1 1])}
                options.rank (1,1) double {mustBeInteger, mustBePositive} = 1
            end

            bc = IMBoundaryCondition(role="active", location=options.location, ...
                component=options.component, signatureSign=options.signatureSign, rank=options.rank);
        end
    end
end

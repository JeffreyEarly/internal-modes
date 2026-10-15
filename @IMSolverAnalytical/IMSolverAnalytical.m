classdef IMSolverAnalytical
    % Request the exact backend declared by an `IMModeBackground` source.
    %
    % This solver marker owns no stratification or physical parameters.
    % `IMModes` obtains them from its background when selecting a supported
    % analytical solution.
    %
    % ```matlab
    % solver = IMSolverAnalytical();
    % ```
    %
    % - Topic: Create an analytical solver
    % - Declaration: classdef IMSolverAnalytical

    methods
        function self = IMSolverAnalytical()
            % Request analytical formulas for a source-specific background.
            %
            % - Topic: Create an analytical solver
            % - Declaration: solver = IMSolverAnalytical()
            % - Returns self: analytical backend marker
        end
    end
end

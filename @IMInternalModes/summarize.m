function summarize(self, solver)
% Print a readable internal-mode EVP summary.
%
% `summarize` first prints the canonical EVP report, then appends the
% internal-mode interpretation: solved formulation, physical constants,
% equivalent-depth mapping, and factory-specific parameter names.
%
% ```matlab
% evp.summarize()
% ```
%
% prints output like:
%
% ```text
% Internal-mode context
%   formulation: G
%   f0: 0 s^-1
%   g: 9.81 m s^-2
%   equivalent depth: h = hFromEigenvalue(lambda)
%   factory parameters: none
% ```
%
% - Topic: Summarize internal-mode EVPs
% - Declaration: summarize(evp,solver)
% - Parameter solver: optional solver for grid-level coefficient and mode-selection assessment
arguments
    self IMInternalModes
    solver = []
end

if isempty(solver)
    summarize@IMEigenvalueProblem(self);
else
    summarize@IMEigenvalueProblem(self, solver);
end

fprintf('\nInternal-mode context\n');
fprintf('  formulation: %s\n', self.formulation);
fprintf('  f0: %s s^-1\n', formatNumber(self.f0));
fprintf('  g: %s m s^-2\n', formatNumber(self.g));
fprintf('  equivalent depth: h = hFromEigenvalue(lambda)\n');
fprintf('  factory parameters: %s\n', factoryParameterText(self.parameters));
end

function text = factoryParameterText(parameters)
names = setdiff(string(fieldnames(parameters)), ["f0", "g", "formulation"], "stable");
if isempty(names)
    text = "none";
else
    text = join(names, ", ");
end
end

function text = formatNumber(value)
if value == 0
    text = "0";
else
    text = string(sprintf('%.6g', value));
end
end

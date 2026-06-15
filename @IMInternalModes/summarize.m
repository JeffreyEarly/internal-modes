function summarize(self, solver)
% Print a readable internal-mode EVP summary.
%
% `summarize` first prints the canonical EVP report, then appends the
% internal-mode interpretation: solved formulation, physical constants,
% equivalent-depth mapping, diagnostic-variable relation, factory-specific
% parameter names, and the available `F` and `G` inner products.
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
%
% Internal-mode variables
%   solved formulation: G
%   diagnostic variable: F
%   diagnostic relation: F_j(z) = h_j dG_j/dz(z)
%
% Internal-mode inner products
%   F
%     interior weight: 1
%     endpoint terms: none
%     availability: interiorOnly
%     reason: For hydrostatic G modes with G=0 at both endpoints, the diagnostic F inner product is the interior F integral.
%   G
%     interior weight: N^2(z)/g
%     endpoint terms: none
%     availability: interiorOnly
%     reason: The solved formulation has no endpoint metric terms, so the inner product is the interior integral only.
% ```
%
% For a fixed-frequency wave EVP, the diagnostic `F` inner product is
% reported as unavailable until a fixed diagnostic catalog entry is added:
%
% ```text
% Internal-mode inner products
%   F
%     interior weight: 1
%     endpoint terms: none
%     availability: unknown
%     reason: No fixed diagnostic inner-product catalog entry is installed for this EVP and boundary combination.
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

fprintf('\nInternal-mode variables\n');
fprintf('  solved formulation: %s\n', self.formulation);
fprintf('  diagnostic variable: %s\n', diagnosticVariable(self));
fprintf('  diagnostic relation: %s\n', diagnosticRelationText(self));

fprintf('\nInternal-mode inner products\n');
printInnerProductSpec(self.innerProduct("F"));
printInnerProductSpec(self.innerProduct("G"));
end

function text = factoryParameterText(parameters)
names = setdiff(string(fieldnames(parameters)), ["f0", "g", "formulation"], "stable");
if isempty(names)
    text = "none";
else
    text = join(names, ", ");
end
end

function variable = diagnosticVariable(evp)
if evp.formulation == "G"
    variable = "F";
else
    variable = "G";
end
end

function text = diagnosticRelationText(evp)
if evp.formulation == "G"
    text = "F_j(z) = h_j dG_j/dz(z)";
elseif string(evp.name) == "hydrostaticFModes"
    text = "G_j(z) = -g/N^2(z) dF_j/dz(z)";
else
    text = "G_j(z) = GfromFz(z,dF_j/dz,h_j,ctx)";
end
end

function printInnerProductSpec(spec)
fprintf('  %s\n', spec.variable);
fprintf('    interior weight: %s\n', interiorWeightText(spec.variable));
fprintf('    endpoint terms: %s\n', endpointTermsText([spec.surfaceWeights; spec.bottomWeights]));
fprintf('    availability: %s\n', spec.status);
fprintf('    reason: %s\n', spec.reason);
end

function text = interiorWeightText(variable)
if string(variable) == "F"
    text = "1";
else
    text = "N^2(z)/g";
end
end

function text = endpointTermsText(weights)
if isempty(weights)
    text = "none";
    return;
end

parts = strings(1,numel(weights));
for iWeight = 1:numel(weights)
    weight = weights(iWeight);
    expression = endpointExpression(weight.c, -weight.d, weight.location);
    parts(iWeight) = weight.location + ": " + formatSignedNumber(weight.coefficient) + " * (" + expression + ")^2";
end
text = join(parts, "; ");
end

function expression = endpointExpression(valueCoefficient, fluxCoefficient, location)
labels = ["u(" + location + ")", "p(" + location + ")*du/dz(" + location + ")"];
expression = signedExpression([valueCoefficient fluxCoefficient], labels);
end

function expression = signedExpression(coefficients, labels)
expression = "";
tolerance = 100*eps*max(1,max(abs(coefficients)));
for iCoefficient = 1:numel(coefficients)
    coefficient = coefficients(iCoefficient);
    if abs(coefficient) <= tolerance
        continue;
    end
    term = formatTerm(abs(coefficient), labels(iCoefficient));
    if strlength(expression) == 0
        if coefficient < 0
            expression = "-" + term;
        else
            expression = term;
        end
    elseif coefficient < 0
        expression = expression + " - " + term;
    else
        expression = expression + " + " + term;
    end
end
if strlength(expression) == 0
    expression = "0";
end
end

function text = formatTerm(coefficient, label)
if coefficient == 1
    text = label;
else
    text = formatNumber(coefficient) + "*" + label;
end
end

function text = formatNumber(value)
if value == 0
    text = "0";
else
    text = string(sprintf('%.6g', value));
end
end

function text = formatSignedNumber(value)
if value == 0
    text = "+0";
else
    text = string(sprintf('%+.6g', value));
end
end

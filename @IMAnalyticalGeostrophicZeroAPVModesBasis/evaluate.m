function values = evaluate(self,z,options)
% Evaluate canonical or rotated profiles and their regular diagnostic derivatives.
%
% Uses Fz=-N2 G/g and Gz=-g k^2 F/f0^2. The same relations apply
% after boundary-coordinate rotations, including on neutral intervals.
% Evaluation preserves the basis's canonical or rotated normalization.
% A different normalization requires constructing the explicitly rotated basis.
%
% - Topic: Evaluate geostrophic zero-APV modes
% - Declaration: values = evaluate(self,z,options)
% - Parameter z: physical sample column
% - Parameter options.variable: required channel F or G
% - Parameter options.derivativeOrder: zero or one for F, through two for G, zero by default
% - Parameter options.pages: source pages, preserving order and repeats
% - Parameter options.normalization: current normalization convention; different conventions require an explicit basis rotation
% - Returns values: nZ-by-nColumns-by-nPages profiles
arguments (Input)
    self IMAnalyticalGeostrophicZeroAPVModesBasis
    z (:,1) double {mustBeReal,mustBeFinite}
    options.variable (1,1) string = ""
    options.derivativeOrder (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.pages (1,:) double {mustBeInteger,mustBePositive} = 1:numel(self.k)
    options.normalization (1,1) string = self.normalizationConvention
end
arguments (Output)
    values double
end
variable = options.variable;
derivativeOrder = options.derivativeOrder;
if variable == ""
    error("IMAnalyticalGeostrophicZeroAPVModesBasis:MissingVariable","Specify the channel with variable= followed by a supported channel name.");
end
if ~ismember(variable,["F","G"])
    error("IMAnalyticalGeostrophicZeroAPVModesBasis:UnsupportedVariable","Supported channels are F and G.");
end
if options.normalization ~= self.normalizationConvention
    error("IMAnalyticalGeostrophicZeroAPVModesBasis:UnsupportedNormalization","This basis uses normalization %s. Construct an explicitly rotated basis to change its normalization convention.",self.normalizationConvention);
end
if any(options.pages > numel(self.k))
    error("IMAnalyticalGeostrophicZeroAPVModesBasis:InvalidPage","pages must select existing source wavenumber pages.");
end
if derivativeOrder == 0
    values = self.(variable)(z,pages=options.pages);
elseif variable == "F" && derivativeOrder == 1
    values = -self.N2(z)/self.g.*self.G(z,pages=options.pages);
elseif variable == "G" && derivativeOrder <= 2
    factor = reshape(-self.g*self.k(options.pages).^2/self.f0^2,1,1,[]);
    values = factor.*self.evaluate(z,variable="F",derivativeOrder=derivativeOrder-1,pages=options.pages);
else
    error("IMAnalyticalGeostrophicZeroAPVModesBasis:UnsupportedDerivative","F supports derivative orders zero and one; G supports orders zero through two.");
end
end

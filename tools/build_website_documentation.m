function build_website_documentation(options)
arguments
    options.rootDir = ".."
    options.rebuildTutorials (1,1) logical = false
end

rootDir = char(java.io.File(char(options.rootDir)).getCanonicalPath());
buildFolder = fullfile(rootDir, "docs");
sourceFolder = fullfile(rootDir, "Documentation", "WebsiteDocumentation");
tutorialSources = {
    fullfile(rootDir, "Examples", "Tutorials", "InternalModesBasics.m")
};
previousTutorialBuildFolder = "";

if ~isfolder(sourceFolder)
    error("build_website_documentation:SourceMissing", ...
        "Could not find source documentation at %s", sourceFolder);
end

ensureDocumentationToolingOnPath(rootDir);

if ~exist("ClassDocumentation", "class") || ~exist("TutorialDocumentation", "class") || ~exist("rebuildWebsiteDocumentationFromSource", "file")
    error("build_website_documentation:MissingDocumentationTools", ...
        "Could not find ClassDocumentation, TutorialDocumentation, or rebuildWebsiteDocumentationFromSource on the MATLAB path.");
end

if isfolder(fullfile(buildFolder, "tutorials"))
    previousTutorialBuildFolder = tempname();
    mkdir(previousTutorialBuildFolder);
    copyfile(fullfile(buildFolder, "tutorials"), fullfile(previousTutorialBuildFolder, "tutorials"));
end

if isfolder(buildFolder)
    rmdir(buildFolder, "s");
end
mkdir(buildFolder);

addPackageToPath(rootDir);

tutorialDocumentation = TutorialDocumentation.documentationFromSourceFiles( ...
    tutorialSources, ...
    buildFolder=buildFolder, ...
    websiteRootURL="internal-modes/", ...
    websiteFolder="tutorials", ...
    sourceRoot=rootDir, ...
    previousBuildFolder=previousTutorialBuildFolder, ...
    rebuildTutorials=options.rebuildTutorials);
preservedTutorialDirectories = unique(string({tutorialDocumentation.preservedAssetDirectoryRelativeToBuildFolder}))';

rebuildWebsiteDocumentationFromSource( ...
    sourceFolder, ...
    buildFolder, ...
    string.empty(0, 1), ...
    preservedRelativeDirectories=preservedTutorialDirectories);

changelogPath = fullfile(rootDir, "CHANGELOG.md");
if isfile(changelogPath)
    header = "---" + newline + ...
             "layout: default" + newline + ...
             "title: Version History" + newline + ...
             "nav_order: 100" + newline + ...
             "---" + newline + newline;
    versionHistoryText = header + fileread(changelogPath);
    versionHistoryFilePath = fullfile(buildFolder, "version-history.md");
    fid = fopen(versionHistoryFilePath, "w");
    assert(fid ~= -1, "Could not open version-history.md for writing");
    fwrite(fid, versionHistoryText);
    fclose(fid);
end

TutorialDocumentation.writeMarkdownIndex( ...
    tutorialDocumentation, ...
    buildFolder=buildFolder, ...
    websiteFolder="tutorials", ...
    nav_order=3, ...
    description="Start here for a minimal workflow with InternalModesSpectral before moving on to the full class reference.");
arrayfun(@(a) a.writeToFile(), tutorialDocumentation)
clear tutorialDocumentation
if previousTutorialBuildFolder ~= "" && isfolder(previousTutorialBuildFolder)
    rmdir(previousTutorialBuildFolder, "s");
end

evalin("base", "clear classes");
evalin("base", "rehash");

websiteRootURL = "internal-modes/";
classFolderName = "Class documentation";
classGroups = {
    struct( ...
        "parent", "Core classes", ...
        "grandparent", classFolderName, ...
        "websiteFolder", "classes/core-classes", ...
        "classes", {{"InternalModes", "InternalModesBase"}} ...
    )
    struct( ...
        "parent", "Numerical solvers", ...
        "grandparent", classFolderName, ...
        "websiteFolder", "classes/numerical-solvers", ...
        "classes", {{"InternalModesSpectral", "InternalModesWKBSpectral", "InternalModesAdaptiveSpectral", "InternalModesDensitySpectral", "InternalModesFiniteDifference"}} ...
    )
    struct( ...
        "parent", "Analytical and asymptotic models", ...
        "grandparent", classFolderName, ...
        "websiteFolder", "classes/analytical-models", ...
        "classes", {{"InternalModesConstantStratification", "InternalModesExponentialStratification", "InternalModesWKB", "InternalModesWKBHydrostatic"}} ...
    )
    struct( ...
        "parent", "Supporting types", ...
        "grandparent", classFolderName, ...
        "websiteFolder", "classes/supporting-types", ...
        "classes", {{"Normalization", "UpperBoundary", "LowerBoundary"}} ...
    )
};

for iGroup = 1:numel(classGroups)
    group = classGroups{iGroup};
    classDocumentation = ClassDocumentation.empty(numel(group.classes), 0);
    for iName = 1:numel(group.classes)
        className = group.classes{iName};
        classDocumentation(iName) = ClassDocumentation( ...
            className, ...
            nav_order=iName, ...
            websiteRootURL=websiteRootURL, ...
            buildFolder=buildFolder, ...
            websiteFolder=group.websiteFolder, ...
            parent=group.parent, ...
            grandparent=group.grandparent, ...
            excludedSuperclasses=excludedSuperclassesForClass(className), ...
            excludedMethodNames=excludedMethodNamesForClass(className));
    end
    arrayfun(@(a) a.writeToFile(), classDocumentation)
end

trimTrailingWhitespaceInMarkdown(buildFolder)
end

function excludedSuperclasses = excludedSuperclassesForClass(className)
switch string(className)
    case {"InternalModes", "InternalModesBase"}
        excludedSuperclasses = {"handle"};
    case {"InternalModesConstantStratification", "InternalModesExponentialStratification", "InternalModesFiniteDifference"}
        excludedSuperclasses = {"handle", "InternalModesBase"};
    case "InternalModesSpectral"
        excludedSuperclasses = {"handle", "InternalModesBase"};
    case {"InternalModesDensitySpectral", "InternalModesWKBSpectral", "InternalModesWKB", "InternalModesWKBHydrostatic"}
        excludedSuperclasses = {"handle", "InternalModesBase", "InternalModesSpectral"};
    case "InternalModesAdaptiveSpectral"
        excludedSuperclasses = {"handle", "InternalModesBase", "InternalModesSpectral", "InternalModesWKBSpectral"};
    otherwise
        excludedSuperclasses = {};
end
end

function excludedMethodNames = excludedMethodNamesForClass(className)
excludedMethodNames = string.empty(0, 1);
classMetadata = meta.class.fromName(className);
if isempty(classMetadata)
    return;
end

methodMetadata = classMetadata.MethodList;
methodNames = reshape(string({methodMetadata.Name}), [], 1);

if ismember(string(className), ["Normalization", "UpperBoundary", "LowerBoundary"])
    excludedMethodNames = methodNames;
    return;
end

definingClasses = reshape(arrayfun(@(method) string(method.DefiningClass.Name), methodMetadata), [], 1);
hiddenMethods = reshape([methodMetadata.Hidden], [], 1);

excludedMethodNames = unique(methodNames( ...
    startsWith(methodNames, "get.") | ...
    startsWith(methodNames, "set.") | ...
    hiddenMethods | ...
    definingClasses ~= string(className)));
end

function trimTrailingWhitespaceInMarkdown(rootFolder)
markdownFiles = dir(fullfile(rootFolder, "**", "*.md"));
for iFile = 1:numel(markdownFiles)
    filePath = fullfile(markdownFiles(iFile).folder, markdownFiles(iFile).name);
    fileText = fileread(filePath);
    trimmedText = regexprep(fileText, "[ \t]+(\r?\n)", "$1");
    if ~strcmp(fileText, trimmedText)
        fid = fopen(filePath, "w");
        assert(fid ~= -1, "Could not open markdown file for writing");
        fwrite(fid, trimmedText);
        fclose(fid);
    end
end
end

function addPackageToPath(repoRoot)
folderPaths = string.empty(0, 1);
metadataPath = fullfile(repoRoot, "resources", "mpackage.json");
if isfile(metadataPath)
    metadata = jsondecode(fileread(metadataPath));
    if isfield(metadata, "folders")
        packageFolders = strings(length(metadata.folders), 1);
        for iFolder = 1:length(metadata.folders)
            packageFolders(iFolder) = string(fullfile(repoRoot, metadata.folders(iFolder).path));
        end
        folderPaths = [folderPaths; packageFolders];
    end
end
folderPaths(end+1, 1) = string(repoRoot);

siblingCandidates = [
    fullfile(repoRoot, "..", "spline-core")
    fullfile(repoRoot, "..", "distributions")
    fullfile(repoRoot, "..", "chebfun")
];
for iCandidate = 1:numel(siblingCandidates)
    if isfolder(siblingCandidates(iCandidate))
        folderPaths(end+1, 1) = string(char(java.io.File(siblingCandidates(iCandidate)).getCanonicalPath())); %#ok<AGROW>
    end
end

folderPaths = unique(folderPaths, "stable");
existingFolders = isfolder(folderPaths);
if any(existingFolders)
    foldersToAdd = cellstr(folderPaths(existingFolders));
    addpath(foldersToAdd{:});
end
end

function ensureDocumentationToolingOnPath(repoRoot)
if exist("ClassDocumentation", "class") && exist("rebuildWebsiteDocumentationFromSource", "file")
    return;
end

candidateFolders = [
    fullfile(repoRoot, "..", "class-docs")
    latestClassDocumentationFolder(repoRoot)
];
candidateFolders = unique(candidateFolders(strlength(candidateFolders) > 0), "stable");
existingFolders = candidateFolders(isfolder(candidateFolders));
if ~isempty(existingFolders)
    foldersToAdd = cellstr(existingFolders);
    addpath(foldersToAdd{:});
end
end

function folderPath = latestClassDocumentationFolder(repoRoot)
folderPath = "";
matches = dir(fullfile(repoRoot, "..", "OceanKit", "ClassDocumentation-*"));
if isempty(matches)
    return;
end

[~, sortIndices] = sort(string({matches.name}), "descend");
folderPath = string(fullfile(matches(sortIndices(1)).folder, matches(sortIndices(1)).name));
end

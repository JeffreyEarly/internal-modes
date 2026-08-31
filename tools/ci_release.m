function ci_release(options)
arguments
    options.rootDir = ".."
    options.bumpType = "none"
    options.notes string = ""
    options.shouldBuildWebsiteDocumentation = false
    options.shouldPackageForDistribution = false
    options.dist_folder = "dist"
    options.excluded_dist_folders = [".git", ".github", "docs", "tools", "Documentation"]
end
%CI_RELEASE CI entry point for MPM release.
%   CI_RELEASE(options) where
%       options.bumpType is "patch", "minor", or "major". If it is left
%
%   Steps:
%     1) Bump version in resources/mpackage.json
%     2) Run custom documentation build
%     3) Export package root to dist/<name>-<version> for MPM repo

if ~isfield(options,"rootDir")
    options.rootDir = pwd;
end

%% bump the version number, if requested
package = matlab.mpm.Package(options.rootDir);
oldVersion = package.Version;
newVersion = oldVersion;

if options.bumpType ~= "none"
    major = oldVersion.Major;
    minor = oldVersion.Minor;
    patch = oldVersion.Patch;
    switch options.bumpType
        case "major"
            major = major + 1;
            minor = 0;
            patch = 0;
        case "minor"
            minor = minor + 1;
            patch = 0;
        case "patch"
            patch = patch + 1;
    end

    newVersion = matlab.mpm.Version(major,minor,patch);
    package.Version = newVersion;

    fprintf('Bumping version: %s -> %s (%s)\n', string(oldVersion), string(newVersion), options.bumpType);

    %% If we bumped the version, and were handed notes, record that
    changelogPath = fullfile(options.rootDir, "CHANGELOG.md");
    update_changelog(changelogPath,options.notes,string(newVersion));
end

%% 2) Run your custom documentation build
if options.shouldBuildWebsiteDocumentation
    % Replace this with your actual doc build entry point
    % e.g. waveVortexDiagnostics_build_docs, or build_docs
    if exist("build_website_documentation","file")
        fprintf('Running documentation builder\n');
        build_website_documentation(rootDir=options.rootDir);
    end
end

%% 3) Export package root to dist/<name>-<version> for MPM repo

if options.shouldPackageForDistribution == true
    pkgName = package.Name;
    distDir = string(options.dist_folder);
    if ~java.io.File(char(distDir)).isAbsolute()
        distDir = fullfile(options.rootDir,distDir);
    end
    if ~isfolder(distDir)
        mkdir(distDir);
    end

    pkgFolderName = pkgName + "-" + string(newVersion);
    targetRoot    = fullfile(distDir, pkgFolderName);

    % Clean any stale output
    if isfolder(targetRoot)
        rmdir(targetRoot, "s");
    end

    fprintf('Exporting package root to %s\n', targetRoot);
    copyfile(options.rootDir, targetRoot);

    % Strip CI-only junk from the exported package
    % (best-effort: ignore errors if these don't exist)
    for iFolder = 1:length(options.excluded_dist_folders)
        try
            rmdir(fullfile(targetRoot, options.excluded_dist_folders(iFolder)), "s");
        catch
        end
    end

    try
        rmdir(fullfile(targetRoot, "dist"), "s");
    catch
    end
    finderFiles = dir(fullfile(targetRoot,"**",".DS_Store"));
    for iFile = 1:numel(finderFiles)
        delete(fullfile(finderFiles(iFile).folder,finderFiles(iFile).name));
    end

    % Write a small metadata file for the GitHub Action
    metaPath = fullfile(distDir, "mpm_release_metadata.txt");
    fid = fopen(metaPath, "w");
    assert(fid ~= -1, "Could not open metadata file for writing");
    fprintf(fid, "NAME=%s\nVERSION=%s\nFOLDER=%s\n", ...
        pkgName, string(newVersion), pkgFolderName);
    fclose(fid);

    fprintf('ci_release complete: %s %s\n', pkgName, string(newVersion));
end
end

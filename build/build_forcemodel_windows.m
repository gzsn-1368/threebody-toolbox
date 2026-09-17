%BUILD_FORCEMODEL_WINDOWS Build forcemodel on 64-bit Windows.
%   Microsoft Visual C++ must be selected with MEX -SETUP C. If CL.EXE and
%   LIB.EXE are not on PATH, start MATLAB from an x64 Native Tools Command
%   Prompt for Visual Studio before running this script.

buildDir = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(buildDir);
cspiceRoot = fullfile(repositoryRoot, 'thirdparty', 'cspice');
cspiceSource = fullfile(cspiceRoot, 'src', 'cspice');
cspiceInclude = fullfile(cspiceRoot, 'include');
artifactDir = fullfile(buildDir, 'artifacts', computer('arch'));
objectDir = fullfile(artifactDir, 'cspice_objects');
cspiceLibrary = fullfile(artifactDir, 'cspice.lib');
outputDir = fullfile(repositoryRoot, '+dynamics', 'private');
forcemodelSource = fullfile(repositoryRoot, 'src', 'forcemodel', 'forcemodel.c');
builder = fullfile(buildDir, 'private', 'build_cspice_windows.ps1');

assert(ispc, 'build_forcemodel:WrongPlatform', ...
    'This build script is for 64-bit Windows.');
assert(isfolder(cspiceSource) && isfolder(cspiceInclude), ...
    'build_forcemodel:MissingCSPICE', 'Bundled CSPICE source is missing.');

if ~isfolder(artifactDir), mkdir(artifactDir); end
if ~isfolder(outputDir), mkdir(outputDir); end

command = sprintf(['powershell -NoProfile -ExecutionPolicy Bypass -File "%s" ' ...
    '-SourceDir "%s" -IncludeDir "%s" -OutputLibrary "%s" -ObjectDir "%s"'], ...
    builder, cspiceSource, cspiceInclude, cspiceLibrary, objectDir);
[status, message] = system(command, '-echo');
if status ~= 0
    error('build_forcemodel:CSPICEBuildFailed', ...
        'CSPICE compilation failed:\n%s', message);
end

mex('-R2018a', ['-I', cspiceInclude], forcemodelSource, cspiceLibrary, ...
    '-outdir', outputDir, '-output', 'forcemodel');
fprintf('Built %s\n', fullfile(outputDir, ['forcemodel.', mexext]));

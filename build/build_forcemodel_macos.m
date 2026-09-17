%BUILD_FORCEMODEL_MACOS Build forcemodel on Apple Silicon or Intel macOS.

buildDir = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(buildDir);
cspiceRoot = fullfile(repositoryRoot, 'thirdparty', 'cspice');
cspiceSource = fullfile(cspiceRoot, 'src', 'cspice');
cspiceInclude = fullfile(cspiceRoot, 'include');
artifactDir = fullfile(buildDir, 'artifacts', computer('arch'));
objectDir = fullfile(artifactDir, 'cspice_objects');
cspiceLibrary = fullfile(artifactDir, 'libcspice.a');
outputDir = fullfile(repositoryRoot, '+dynamics', 'private');
forcemodelSource = fullfile(repositoryRoot, 'src', 'forcemodel', 'forcemodel.c');
builder = fullfile(buildDir, 'private', 'build_cspice_unix.sh');

assert(ismac, 'build_forcemodel:WrongPlatform', ...
    'This build script is for macOS.');
assert(isfolder(cspiceSource) && isfolder(cspiceInclude), ...
    'build_forcemodel:MissingCSPICE', 'Bundled CSPICE source is missing.');

if ~isfolder(artifactDir), mkdir(artifactDir); end
if ~isfolder(outputDir), mkdir(outputDir); end

command = sprintf('CC=clang /bin/sh "%s" "%s" "%s" "%s" "%s"', ...
    builder, cspiceSource, cspiceInclude, cspiceLibrary, objectDir);
[status, message] = system(command, '-echo');
if status ~= 0
    error('build_forcemodel:CSPICEBuildFailed', ...
        'CSPICE compilation failed:\n%s', message);
end

mex('-R2018a', ['-I', cspiceInclude], forcemodelSource, cspiceLibrary, ...
    '-outdir', outputDir, '-output', 'forcemodel');
fprintf('Built %s\n', fullfile(outputDir, ['forcemodel.', mexext]));

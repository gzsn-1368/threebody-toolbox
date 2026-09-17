%BUILD_FORCEMODEL Build the forcemodel MEX file for the current platform.
%   Run this script from any folder. It dispatches to the macOS, Linux, or
%   Windows build script and writes the MEX file to +dynamics/private.

buildDir = fileparts(mfilename('fullpath'));

switch computer('arch')
    case {'maca64', 'maci64'}
        run(fullfile(buildDir, 'build_forcemodel_macos.m'));
    case 'glnxa64'
        run(fullfile(buildDir, 'build_forcemodel_linux.m'));
    case 'win64'
        run(fullfile(buildDir, 'build_forcemodel_windows.m'));
    otherwise
        error('build_forcemodel:UnsupportedPlatform', ...
            'Unsupported MATLAB architecture: %s.', computer('arch'));
end

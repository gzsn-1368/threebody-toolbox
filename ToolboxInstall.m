%TOOLBOXINSTALL Install this Threebody Toolbox open-source release.
%   Run this script from any current folder. The script resolves its own
%   location, adds the required toolbox and data paths, configures the SPICE
%   kernel path, attempts to build the ephemeris MEX file when necessary,
%   and saves the updated MATLAB search path.

toolboxRoot = fileparts(mfilename('fullpath'));
addpath(toolboxRoot);
addpath(fullfile(toolboxRoot, 'kernels'));
addpath(fullfile(toolboxRoot, '+orbits', 'crtbp_data'));
run(fullfile(toolboxRoot, 'startup.m'));

if isempty(which('dynamics.crtbp')) || ...
        isempty(which('dynamics.ephem')) || ...
        isempty(which('dynamics.ephemoptions')) || ...
        isempty(which('orbits.crtbp_initstates'))
    error('ToolboxInstall:VerificationFailed', ...
        'The Threebody Toolbox open-source release could not be resolved on the MATLAB path.');
end

forcemodelMex = fullfile(toolboxRoot, '+dynamics', 'private', ...
    ['forcemodel.', mexext]);
if ~isfile(forcemodelMex)
    fprintf('The ephemeris MEX file was not found; attempting to build it.\n');
    try
        run(fullfile(toolboxRoot, 'build', 'build_forcemodel.m'));
        if ~isfile(forcemodelMex)
            warning('ToolboxInstall:ForcemodelBuildFailed', ...
                ['The platform build script completed without producing %s. ' ...
                'dynamics.ephem will remain unavailable.'], forcemodelMex);
        end
    catch buildError
        warning('ToolboxInstall:ForcemodelBuildFailed', ...
            ['The ephemeris MEX file could not be built automatically. ' ...
            'dynamics.ephem will remain unavailable until it is built manually.\n' ...
            'Reason: %s'], buildError.message);
    end
end

savepath;

fprintf('Threebody Toolbox open-source release installed.\n');
fprintf('Root: %s\n', toolboxRoot);

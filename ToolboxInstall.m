%TOOLBOXINSTALL Add this Threebody Toolbox open-source release to the MATLAB path.
%   Run this script once per MATLAB session from any current folder. The
%   script resolves its own location and adds only the repository root. It
%   intentionally does not call SAVEPATH, which would change global MATLAB
%   configuration.

toolboxRoot = fileparts(mfilename('fullpath'));
addpath(toolboxRoot);

if isempty(which('dynamics.crtbp')) || ...
        isempty(which('orbits.crtbp_initstates'))
    error('ToolboxInstall:VerificationFailed', ...
        'The Threebody Toolbox open-source release could not be resolved on the MATLAB path.');
end

fprintf('Threebody Toolbox open-source release added to the MATLAB path.\n');
fprintf('Root: %s\n', toolboxRoot);

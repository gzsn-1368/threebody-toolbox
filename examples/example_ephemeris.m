%% Moon-centred ephemeris dynamics
% Build the forcemodel MEX file before running this example.

clear; clc;
toolboxRoot = fileparts(fileparts(mfilename('fullpath')));
run(fullfile(toolboxRoot, 'ToolboxInstall.m'));

forcemodelMex = fullfile(toolboxRoot, '+dynamics', 'private', ...
    ['forcemodel.', mexext]);
if ~isfile(forcemodelMex)
    error('example_ephemeris:MissingMEX', ...
        ['Build the ephemeris MEX file first with ' ...
        'run(fullfile(toolboxRoot, ''build'', ''build_forcemodel.m'')).']);
end

epoch = 0; % SPICE ephemeris seconds past J2000 TDB
state = [2000; 0; 0; 0; 1.6; 0]; % km and km/s, Moon-centred J2000
options = dynamics.ephemoptions( ...
    'EarthHarmonics', [2, 0], ...
    'LunarHarmonics', [10, 10], ...
    'FourthBody', 'sun', ...
    'SolarRadiationPressure', 'none');

[stateDerivative, stateJacobian] = dynamics.ephem( ...
    epoch, state, zeros(3, 1), options);

disp(stateDerivative);
disp(stateJacobian);

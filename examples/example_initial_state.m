%EXAMPLE_INITIAL_STATE Generate and propagate a CRTBP halo-orbit seed.
%   This example demonstrates the public initial-state generator and an
%   independent propagation with the public CRTBP dynamics function.

clear; clc; close all; format long;

exampleRoot = fileparts(mfilename('fullpath'));
toolboxRoot = fileparts(exampleRoot);
addpath(toolboxRoot);

odeOptions = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
[initialStates, periods, info] = orbits.crtbp_initstates( ...
    'halo', 'L2', 'northern', 15000, 0, ...
    'UsePlot', false, 'OdeOptions', odeOptions);

state0 = initialStates(1, :).';
[epochs, states] = ode113(@dynamics.crtbp, [0, periods(1)], state0, odeOptions);

fprintf('Generated %d initial state(s).\n', size(initialStates, 1));
fprintf('Period: %.15g normalized time units\n', periods(1));
fprintf('Jacobi constant: %.15g\n', info.jacobiConstants(1));
fprintf('Stability index: %.15g\n', info.stabilityIndices(1));
fprintf('Endpoint closure residual: %.3e\n', ...
    norm(states(end, :).'-state0));

figure;
plot3(states(:, 1), states(:, 2), states(:, 3), 'LineWidth', 1.0);
grid on; axis equal;
xlabel('x [LU]'); ylabel('y [LU]'); zlabel('z [LU]');
title('Earth--Moon CRTBP halo orbit');

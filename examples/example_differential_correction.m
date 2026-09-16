%EXAMPLE_DIFFERENTIAL_CORRECTION Correct a perturbed CRTBP periodic-orbit seed.
%   This example perturbs an interpolated Halo-orbit seed and recovers a
%   periodic solution with single-shooting differential correction.

clear; clc; close all; format long;

exampleRoot = fileparts(mfilename('fullpath'));
toolboxRoot = fileparts(exampleRoot);
addpath(toolboxRoot);

odeOptions = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
[seedStates, seedPeriods] = orbits.crtbp_initstates( ...
    'halo', 'L2', 'northern', 15000, 0, ...
    'UsePlot', false, 'OdeOptions', odeOptions);

seedState = seedStates(1, :);
seedPeriod = seedPeriods(1);
perturbedState = seedState;
perturbedState(1) = perturbedState(1) + 1e-5;
perturbedPeriod = seedPeriod * (1 + 1e-6);

[correctedState, correctedPeriod, residual, exitflag] = ...
    orbits.crtbp_diffcorrection(perturbedState, perturbedPeriod, ...
    [1, 5], [2, 4, 6], 'IsPeriodFree', true, ...
    'UseHalfPeriod', true, 'Tolerance', 1e-9, ...
    'MaxIterations', 100, 'OdeOptions', odeOptions);

[epochs, states] = ode113(@dynamics.crtbp, [0, correctedPeriod], ...
    correctedState(:), odeOptions);
closureResidual = norm(states(end, :).'-correctedState(:));

fprintf('Differential-correction exit flag: %d\n', exitflag);
fprintf('Corrected period: %.15g normalized time units\n', correctedPeriod);
fprintf('Correction residual: %.3e\n', residual);
fprintf('Full-period closure residual: %.3e\n', closureResidual);

figure;
plot3(states(:, 1), states(:, 2), states(:, 3), 'LineWidth', 1.0);
grid on; axis equal;
xlabel('x [LU]'); ylabel('y [LU]'); zlabel('z [LU]');
title('Corrected Earth--Moon CRTBP halo orbit');

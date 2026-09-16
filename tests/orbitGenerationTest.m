classdef orbitGenerationTest < matlab.unittest.TestCase
    %ORBITGENERATIONTEST Test the public orbit-generation interfaces.

    methods (TestClassSetup)
        function addSourceFolder(testCase)
            testFolder = fileparts(mfilename('fullpath'));
            toolboxRoot = fileparts(testFolder);
            testCase.applyFixture(matlab.unittest.fixtures.PathFixture(toolboxRoot));
        end
    end

    methods (Test)
        function testBundledDataAndHaloGeneration(testCase)
            odeOptions = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
            amplitudeLimits = orbits.crtbp_initstates('halo', 'L2', 'northern');
            [initialStates, periods, info] = orbits.crtbp_initstates( ...
                'halo', 'L2', 'northern', 15000, 0, ...
                'UsePlot', false, 'OdeOptions', odeOptions);

            testCase.verifySize(amplitudeLimits, [1, 2]);
            testCase.verifySize(initialStates, [1, 6]);
            testCase.verifySize(periods, [1, 1]);
            testCase.verifyEqual(numel(info.fullstates), 1);
            testCase.verifyGreaterThan(periods(1), 0);
            testCase.verifyTrue(all(isfinite(initialStates), 'all'));
        end

        function testDifferentialCorrection(testCase)
            odeOptions = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
            [seedStates, seedPeriods] = orbits.crtbp_initstates( ...
                'halo', 'L2', 'northern', 15000, 0, ...
                'UsePlot', false, 'OdeOptions', odeOptions);
            perturbedState = seedStates(1, :);
            perturbedState(1) = perturbedState(1) + 1e-5;

            [correctedState, correctedPeriod, residual, exitflag] = ...
                orbits.crtbp_diffcorrection(perturbedState, ...
                seedPeriods(1) * (1 + 1e-6), [1, 5], [2, 4, 6], ...
                'IsPeriodFree', true, 'UseHalfPeriod', true, ...
                'Tolerance', 1e-9, 'OdeOptions', odeOptions);

            testCase.verifyEqual(exitflag, 1);
            testCase.verifySize(correctedState, [1, 6]);
            testCase.verifyGreaterThan(correctedPeriod, 0);
            testCase.verifyLessThan(residual, 1e-9);
        end

        function testConstants(testCase)
            [mu, earthRadius, moonRadius] = constants( ...
                'MassUnitEM', 'EarthRadius', 'MoonRadius');

            testCase.verifyGreaterThan(mu, 0);
            testCase.verifyGreaterThan(earthRadius, moonRadius);
            testCase.verifyGreaterThan(moonRadius, 0);
        end
    end
end

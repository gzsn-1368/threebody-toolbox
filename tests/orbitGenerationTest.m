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

        function testHaloL1AmplitudeGapRegression(testCase)
            odeOptions = odeset('AbsTol', 1e-12, 'RelTol', 1e-12);
            requestedAmplitudes = [73000; 76000; 79000];
            [initialStates, periods, info] = orbits.crtbp_initstates( ...
                'halo', 'L1', 'northern', requestedAmplitudes, 0, ...
                'UseWarnings', false, 'UsePlot', false, ...
                'OdeOptions', odeOptions);

            testCase.verifySize(initialStates, [3, 6]);
            testCase.verifySize(periods, [3, 1]);
            testCase.verifyGreaterThan(periods, zeros(3, 1));
            testCase.verifyLessThan(periods, 5 * ones(3, 1));
            testCase.verifyLessThan(max(abs(initialStates(:, 1:3)), [], 2), ...
                2 * ones(3, 1));

            lengthUnitKm = constants('LengthUnitEM') / 1000;
            for ii = 1:numel(requestedAmplitudes)
                states = info.fullstates{ii};
                closureError = norm(states(end, 1:6) - states(1, 1:6));
                actualAmplitude = max(abs(states(:, 3))) * lengthUnitKm;
                testCase.verifyLessThan(closureError, 1e-7);
                testCase.verifyEqual(actualAmplitude, ...
                    requestedAmplitudes(ii), 'AbsTol', 1);
            end
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

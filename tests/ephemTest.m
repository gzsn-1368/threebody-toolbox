classdef ephemTest < matlab.unittest.TestCase
    methods (Test)
        function startupConfiguresKernelPath(testCase)
            testFile = mfilename('fullpath');
            toolboxRoot = fileparts(fileparts(testFile));
            run(fullfile(toolboxRoot, 'startup.m'));
            testCase.verifyEqual(getenv('SPICE_KERNEL_PATH'), ...
                fullfile(toolboxRoot, 'kernels'));
        end

        function optionsAreSelfContained(testCase)
            options = dynamics.ephemoptions( ...
                'EarthGravityModel', 'EGM96', ...
                'LunarGravityModel', 'LP100K', ...
                'SolarRadiationPressure', 'none');
            testCase.verifyEqual(options.EarthGravityField.modelVer, 2);
            testCase.verifyEqual(options.LunarGravityField.modelVer, 3);
            testCase.verifySize(options.FourthBody, [1, 9]);
            testCase.verifyFalse(isfield(options, 'Mass'));
            testCase.verifyFalse(isfield(options, 'InertiaMatrix'));
            testCase.verifyFalse(isfield(options, 'Plates'));
        end

        function attitudeOptionsAreNotSupported(testCase)
            testCase.verifyError(@() dynamics.ephemoptions( ...
                'SolarRadiationPressure', 'n-plate'), ...
                'checkEphemOptions:InvalidParam');
        end

        function attitudeStateIsNotSupported(testCase)
            options = dynamics.ephemoptions('SolarRadiationPressure', 'none');
            testCase.verifyError(@() dynamics.ephem( ...
                0, zeros(13, 1), zeros(3, 1), options), ...
                'ephem:InvalidInput');
        end

        function derivativeAndJacobian(testCase)
            testFile = mfilename('fullpath');
            toolboxRoot = fileparts(fileparts(testFile));
            forcemodelMex = fullfile(toolboxRoot, '+dynamics', 'private', ...
                ['forcemodel.', mexext]);
            testCase.assumeTrue(isfile(forcemodelMex), ...
                'The forcemodel MEX file has not been built.');
            options = dynamics.ephemoptions( ...
                'LunarHarmonics', [2, 0], ...
                'SolarRadiationPressure', 'none');
            [stateDerivative, jacobian] = dynamics.ephem( ...
                0, [2000; 0; 0; 0; 1.6; 0], zeros(3, 1), options);
            testCase.verifySize(stateDerivative, [6, 1]);
            testCase.verifySize(jacobian, [6, 6]);
            testCase.verifyTrue(all(isfinite(stateDerivative)));
            testCase.verifyTrue(all(isfinite(jacobian), 'all'));
        end
    end
end

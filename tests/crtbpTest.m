classdef crtbpTest < matlab.unittest.TestCase
    %CRTBPTEST Test the public CRTBP dynamics interface.

    methods (Test)
        function testSixElementState(testCase)
            state = [1.05; 0.02; 0.01; 0.0; 0.1; -0.02];
            [stateDerivative, jacobian] = dynamics.crtbp(0, state);

            testCase.verifySize(stateDerivative, [6, 1]);
            testCase.verifySize(jacobian, [6, 6]);
            testCase.verifyTrue(all(isfinite(stateDerivative)));
            testCase.verifyTrue(all(isfinite(jacobian), 'all'));
        end

        function testFortyTwoElementVariationalState(testCase)
            state = [1.05; 0.02; 0.01; 0.0; 0.1; -0.02];
            variationalState = [state; reshape(eye(6), 36, 1)];
            stateDerivative = dynamics.crtbp(0, variationalState);

            testCase.verifySize(stateDerivative, [42, 1]);
            testCase.verifyEqual(stateDerivative(1:6), ...
                dynamics.crtbp(0, state), AbsTol=1e-14);
        end

        function testOptionalControlAcceleration(testCase)
            state = [1.05; 0.02; 0.01; 0.0; 0.1; -0.02];
            control = [1e-4; -2e-4; 3e-4];
            uncontrolled = dynamics.crtbp(0, state);
            controlled = dynamics.crtbp(0, state, control);

            testCase.verifyEqual(controlled - uncontrolled, ...
                [zeros(3, 1); control], AbsTol=1e-14);
        end

        function testInvalidStateLength(testCase)
            testCase.verifyError(@() dynamics.crtbp(0, zeros(5, 1)), ...
                'crtbp:InvalidInput');
        end
    end
end

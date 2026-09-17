function [initstates, period, residual, exitflag] = crtbp_diffcorrection( ...
        initstates, period, freeIndex, targetIndex, varargin)
%CRTBP_DIFFCORRECTION Correct a periodic-orbit initial state in the CRTBP.
%   [INITSTATES, PERIOD, RESIDUAL, EXITFLAG] = CRTBP_DIFFCORRECTION( ...
%   INITSTATES, PERIOD, FREEINDEX) corrects a CRTBP orbit initial state by
%   single shooting and a state transition matrix. Newton iterations adjust
%   the initial state and, optionally, the period to satisfy periodicity.
%
%   Name-value options may follow TARGETINDEX. 'Tolerance' sets the
%   periodicity tolerance. 'IsPeriodFree' selects whether the period is
%   corrected, and 'UseHalfPeriod' selects half-period symmetry correction.
%   'MaxIterations' and 'OdeOptions' control the iteration limit and ODE
%   solver settings.
%
%   Inputs:
%     INITSTATES     - Six-element initial state to be corrected.
%     PERIOD         - Positive initial period in normalized time units.
%     FREEINDEX      - Indices in 1:6 of the variables that may be corrected.
%     TARGETINDEX    - Indices in 1:6 of the periodicity constraints.
%     Name-value options:
%       'IsPeriodFree'  - Whether to correct the period (default: true).
%       'UseHalfPeriod' - Whether to impose half-period symmetry (default: false).
%       'Tolerance'     - Positive periodicity tolerance (default: 1e-9).
%       'MaxIterations' - Maximum number of iterations (default: 100).
%       'OdeOptions'    - ODE solver options (default: AbsTol and RelTol equal
%                        to 1e-12).
%
%   Outputs:
%     INITSTATES     - Corrected six-element initial state, returned as a row.
%     PERIOD         - Corrected period in normalized time units.
%     RESIDUAL       - Final periodicity residual.
%     EXITFLAG       - 1 if the correction converged; 0 if it stopped before
%                      convergence.
%
%   See also ORBITS.CRTBP_INITSTATES.
%
%   Copyright 2025 Guo Zisen.

    % Parse name-value options.
    p = inputParser;
    addParameter(p, 'IsPeriodFree', true);
    addParameter(p, 'UseHalfPeriod', false);
    addParameter(p, 'Tolerance', 1e-9);
    addParameter(p, 'MaxIterations', 100);
    addParameter(p, 'OdeOptions', odeset('AbsTol', 1e-12, 'RelTol', 1e-12));

    parse(p, varargin{:});
    options = p.Results;

    isPeriodFree = options.IsPeriodFree;
    useHalfPeriod = options.UseHalfPeriod;
    tolerance = options.Tolerance;
    maxiter = options.MaxIterations;
    odeOptions = options.OdeOptions;

    % Validate inputs.
    if ~isnumeric(initstates) || ~isvector(initstates) || ...
            length(initstates) ~= 6
        error('crtbp_diffcorrection:InvalidInput', ...
            'initstates must be a six-element vector.');
    end
    if ~isnumeric(period) || ~isscalar(period) || period < 0
        error('crtbp_diffcorrection:InvalidInput', 'period must be a positive scalar.');
    end
    if ~isnumeric(freeIndex) || ~isvector(freeIndex) || ...
            any(~ismember(freeIndex, 1:6))
        error('crtbp_diffcorrection:InvalidInput', ...
            'freeIndex must be a subset of the integers 1 through 6.');
    end
    if ~isnumeric(targetIndex) || ~isvector(targetIndex) || ...
            any(~ismember(targetIndex, 1:6))
        error('crtbp_diffcorrection:InvalidInput', ...
            'targetIndex must be a subset of the integers 1 through 6.');
    end
    if ~islogical(isPeriodFree)
        error('crtbp_diffcorrection:InvalidInput', ...
            'IsPeriodFree must be a logical scalar.');
    end
    if ~islogical(useHalfPeriod)
        error('crtbp_diffcorrection:InvalidInput', ...
            'UseHalfPeriod must be a logical scalar.');
    end
    if ~isnumeric(tolerance) || ~isscalar(tolerance) || tolerance <= 0
        error('crtbp_diffcorrection:InvalidInput', ...
            'Tolerance must be a positive scalar.');
    end
    if ~isnumeric(maxiter) || ~isscalar(maxiter) || maxiter <= 0
        error('crtbp_diffcorrection:InvalidInput', ...
            'MaxIterations must be a positive scalar.');
    end
    if ~isstruct(odeOptions)
        error('crtbp_diffcorrection:InvalidInput', ...
            'OdeOptions must be a structure.');
    end

    % Start differential correction.
    initstates = initstates(:);
    nfree = length(freeIndex);
    numiter = 1;
    while true
        if useHalfPeriod
            tspan = [0, period/2];
        else
            tspan = [0, period];
        end
        y0 = [initstates; reshape(eye(6), 36, 1)];
        [~, tmpstates] = ode113(@dynamics.crtbp, tspan, y0, odeOptions);
        finalstates = tmpstates(end, 1:6)';
        if isPeriodFree
            tmp = dynamics.crtbp(tspan(end), finalstates); % Endpoint derivative.
        end

        % Compute the residual and check the stopping condition.
        residual = norm(finalstates(targetIndex) - initstates(targetIndex));
        if residual < tolerance
            exitflag = 1;
            break;
        end

        % Compute the correction vector.
        STM = reshape(tmpstates(end, 7:42), 6, 6);
        if isPeriodFree
            if useHalfPeriod
                periodSensitivity = 0.5 * tmp(targetIndex);
            else
                periodSensitivity = tmp(targetIndex);
            end
            L = [STM(targetIndex, freeIndex), periodSensitivity];
        else
            L = STM(targetIndex, freeIndex);
        end 
        if useHalfPeriod
            b = -finalstates(targetIndex);
        else
            b = -finalstates(targetIndex) + initstates(targetIndex);
        end
        correction = pinv(L) * b;

        % Backtrack the Newton step to keep the correction on the local
        % solution branch. This is especially important near folds, where
        % an undamped step can converge to a distant periodic solution.
        stepAccepted = false;
        stepScale = 1;
        for lineIter = 1:24
            trialStates = initstates;
            trialStates(freeIndex) = trialStates(freeIndex) + ...
                stepScale * correction(1:nfree);
            trialPeriod = period;
            if isPeriodFree
                trialPeriod = trialPeriod + stepScale * correction(end);
            end

            if trialPeriod > 0
                trialResidual = computeResidual(trialStates, trialPeriod, ...
                    targetIndex, useHalfPeriod, odeOptions);
                if isfinite(trialResidual) && trialResidual < residual
                    initstates = trialStates;
                    period = trialPeriod;
                    stepAccepted = true;
                    break;
                end
            end
            stepScale = stepScale / 2;
        end

        if ~stepAccepted
            exitflag = 0;
            break;
        end

        numiter = numiter + 1;

        % Stop when the maximum iteration count is exceeded.
        if numiter > maxiter
            exitflag = 0;
            break;
        end
    end
    % Return the state as a row vector.
    initstates = initstates.';
end

function residual = computeResidual(initstates, period, targetIndex, ...
        useHalfPeriod, odeOptions)
    if useHalfPeriod
        tspan = [0, period / 2];
    else
        tspan = [0, period];
    end
    [~, states] = ode113(@dynamics.crtbp, tspan, initstates, odeOptions);
    finalstates = states(end, 1:6).';
    residual = norm(finalstates(targetIndex) - initstates(targetIndex));
end

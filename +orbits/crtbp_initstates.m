function [initstates, periods, info] = crtbp_initstates(orbit, location, ...
        direction, targetamps, targetphases, varargin)
%CRTBP_INITSTATES Generate CRTBP periodic-orbit initial states and periods.
%   [INITSTATES, PERIOD] = CRTBP_INITSTATES(ORBIT, LOCATION/RESONANCE, ...
%   DIRECTION, TARGETAMPS, TARGETPHASES) interpolates the bundled CRTBP orbit
%   data to generate initial states and periods for the requested amplitudes
%   and phases.
%
%   [INITSTATES, PERIODS, INFO] = CRTBP_INITSTATES(...) also returns orbit
%   metadata, including the Jacobi constant, stability index, and sampled
%   states over one period. Some families may return multiple states for the
%   same amplitude; INFO.targetamps and INFO.targetphases preserve the mapping
%   between requested amplitudes, phases, states, and periods.
%
%   Name-value options may follow TARGETPHASES. Use 'UsePlot' to preview the
%   interpolated and corrected orbits, or set 'UseWarnings' to false to silence
%   warnings. Differential-correction options can also be passed through.
%
%   ALIM = CRTBP_INITSTATES(ORBIT, LOCATION, DIRECTION) returns the supported
%   amplitude limits when exactly three inputs and one output are requested.
%   ALIM is a two-element row vector [minimum maximum].
%
%   Inputs:
%     ORBIT        - Orbit-family name, one of:
%                      'lyapunov'      Lyapunov orbit
%                      'halo'          Halo orbit
%                      'nrho'          Near-rectilinear halo orbit (NRHO)
%                      'vertical'      Vertical orbit
%                      'axial'         Axial orbit
%                      'shortperiod'   Short-period orbit
%                      'longperiod'    Long-period orbit
%                      'butterfly'     Butterfly orbit
%                      'dragonfly'     Dragonfly orbit
%                      'dro'           Distant retrograde orbit (DRO)
%                      'dro3d'         Three-dimensional DRO
%                      'dpo'           Distant prograde orbit (DPO)
%                      'lpo'           Low prograde orbit (LoPO)
%                      'resonant'      Resonant orbit
%     LOCATION     - Libration-point or resonance-ratio identifier. For
%     (RESONANCE)    Lyapunov, Halo, NRHO, and Axial families, use 'L1', 'L2',
%                    or 'L3' as applicable. Vertical and Axial families also
%                    support 'L4' and 'L5'; short- and long-period families
%                    use 'L4' or 'L5'. Resonant families use one of '1:1',
%                    '1:2', '1:3', '1:4', '2:1', '2:3', '3:1', '3:2', '3:4',
%                    '4:1', or '4:3'. Other families do not use LOCATION.
%     DIRECTION    - Orbit direction. Halo, NRHO, Axial, Butterfly, Dragonfly,
%                    and 3D DRO families use 'northern' or 'southern'. LoPO
%                    uses 'eastern' or 'western'. Other families do not use
%                    DIRECTION.
%     TARGETAMPS   - Row or column vector of target amplitudes. Multiple
%                    amplitudes may be requested. The supported ranges are:
%                                  Libration point
%                         Orbit       (or ratio)    Minimum       Maximum
%                       lyapunov         L1          3     109741
%                                        L2         17      90503
%                                        L3         35     361870
%                           halo         L1          1     382485
%                                        L2         40      77787
%                                        L3        193     739527
%                           nrho         L1        683       6566
%                                        L2         29      16672
%                                        L3       2293      51378
%                       vertical         L1         75     390176
%                                        L2        324     378459
%                                        L3        164     382856
%                                        L4     218718     384762
%                                        L5     218718     384762
%                          axial         L1          2      28376
%                                        L2         10      36609
%                                        L3       1673     343488
%                                        L4      70741     338689
%                                        L5      70741     338689
%                    shortperiod         L4         14     286635
%                                        L5         14     286635
%                     longperiod         L4         85     210068
%                                        L5         85     210068
%                      butterfly         []      54993     195841
%                      dragonfly         []      24282      78958
%                            dro         []       2797     374360
%                          dro3d         []       3117     306146
%                            dpo         []       3421      26271
%                                        []      36542     287969
%                    lpo(eastern)        []       1330      28705
%                    lpo(western)        []       1961      22408
%                       resonant        1:1      36780      85163
%                                       1:2      33919     374556
%                                       1:3      14356     357793
%                                       1:4      23038     548903
%                                       2:1      32227     198347
%                                       2:3     130719     382656
%                                       3:1       5503     120140
%                                       3:2     125577     278087
%                                       3:4      96410     379247
%                                       4:1      43308     150347
%                                       4:3     489274     562751
%     TARGETPHASES - Row or column vector of initial phases with the same
%                    length as TARGETAMPS. Phases must lie in [0, 2*pi]. A
%                    scalar phase is applied to every requested amplitude.
%     Name-value options:
%       'UseWarnings'            Whether to issue warnings (default: true).
%       'UsePlot'                Whether to plot the generated orbits
%                                (default: false).
%       'InterpolationMethod'    Interpolation method: 'linear', 'makima',
%                                'pchip', or 'spline' (default: 'pchip').
%       'DiffCorrFreeIndex'      Indices of differential-correction free
%                                variables; default depends on ORBIT.
%       'DiffCorrTargetIndex'    Indices of differential-correction target
%                                variables; default depends on ORBIT.
%       'DiffCorrIsPeriodFree'   Whether to correct the period (default: true).
%       'DiffCorrIsHalfPeriod'   Whether to use half-period correction;
%                                default depends on ORBIT.
%       'DiffCorrTolerance'      Positive differential-correction tolerance
%                                (default: 1e-9).
%       'DiffCorrMaxIterations'  Maximum differential-correction iterations
%                                (default: 100).
%       'OdeOptions'             ODE solver options; default AbsTol and RelTol
%                                are 1e-12.
%
%   Outputs:
%   With exactly three inputs and one output:
%     ALIM         - 1-by-2 vector [minimum maximum] of supported amplitudes.
%   Otherwise:
%     INITSTATES   - M-by-6 matrix of initial states for the requested orbits.
%     PERIODS      - M-by-1 vector of corresponding orbital periods.
%     INFO         - Structure with fields:
%                      - fullepochs:       M-by-1 cell array of epochs over
%                                          one period for each orbit.
%                      - fullstates:       M-by-1 cell array of state histories.
%                      - jacobiConstants:  M-by-1 vector of Jacobi constants.
%                      - stabilityIndices: M-by-1 vector of stability indices.
%                      - targetamps:       M-by-1 vector of matched amplitudes.
%                      - targetphases:     M-by-1 vector of matched phases.
%
%   See also ORBITS.CRTBP_DIFFCORRECTION.
%
%   Copyright 2025-2026 Guo Zisen.

    persistent mu
    if isempty(mu)
        mu = constants('MassUnitEM');
    end

    checkOrbitValidity(orbit, location, direction);

    % Load the bundled orbit data.
    [X, T, A, S] = crtbp_loadData(orbit, location, direction);

    minamp = ceil(min(A));
    maxamp = floor(max(A));

    if nargin < 4 && nargout == 1
        initstates = [minamp, maxamp];
        return;
    end

    if ~isnumeric(targetamps)
        error('crtbp_initstates:InvalidInput', 'targetamps must be numeric.');
    end

    if ~isnumeric(targetphases)
        error('crtbp_initstates:InvalidInput', 'targetphases must be numeric.');
    end

    numcases = numel(targetamps);
    if length(targetphases) ~= numcases
        if numcases > 1 && isscalar(targetphases)
            targetphases = targetphases * ones(numcases, 1);
        else
            error('crtbp_initstates:InvalidInput', ...
                'targetamps and targetphases must have the same length.');
        end
    end
    if strcmp(orbit, 'dpo')
        if any(targetamps < minamp) || ...
                any(targetamps < ceil(A(634)) & targetamps > floor(A(635))) || ...
                any(targetamps > maxamp)
            error('crtbp_initstates:InvalidInput', ['targetamps must lie in ' ...
                '%d to %d or %d to %d.'], minamp, floor(A(635)), ...
                ceil(A(634)), maxamp);
        end
    elseif any(targetamps < minamp) || any(targetamps > maxamp)
        error('crtbp_initstates:InvalidInput', ['targetamps must lie between %d and ' ...
            '%d.'], minamp, maxamp);
    end
    if any(targetphases < 0) || any(targetphases > 2*pi)
        error('crtbp_initstates:InvalidInput', ...
            'targetphases must lie between 0 and 2*pi.');
    end

    % Parse name-value options.
    p = inputParser;
    addParameter(p, 'UseWarnings', true);
    addParameter(p, 'UsePlot', false);
    addParameter(p, 'InterpolationMethod', 'pchip');
    addParameter(p, 'DiffCorrFreeIndex', []);
    addParameter(p, 'DiffCorrTargetIndex', []);
    addParameter(p, 'DiffCorrIsPeriodFree', true);
    addParameter(p, 'DiffCorrIsHalfPeriod', []);
    addParameter(p, 'DiffCorrTolerance', 1e-9);
    addParameter(p, 'DiffCorrMaxIterations', 100);
    addParameter(p, 'OdeOptions', odeset('AbsTol', 1e-20, 'RelTol', 1e-13));

    parse(p, varargin{:});
    options = p.Results;
    [usewarnings, useplot, interpMethod, freeIndex, targetIndex, ...
        isPeriodFree, useHalfPeriod, tolerance, maxiter, odeOptions] = ...
        handleOptions(options);

    % Build interpolants over the continuation parameter.
    F.S = S;
    for ii = 1:6
        F.X{ii} = griddedInterpolant(S, X(:, ii), interpMethod);
    end
    F.T = griddedInterpolant(S, T, interpMethod);
    F.A = griddedInterpolant(S, A, interpMethod);

    % Select differential-correction defaults for the orbit family.
    if isempty(freeIndex)
        switch orbit
            case {'lyapunov', 'shortperiod', 'longperiod', 'dro', ...
                    'dpo', 'lpo', 'resonant'}
                freeIndex = 5;
            case {'axial'}
                freeIndex = [5, 6];
            otherwise
                freeIndex = [1, 5];
        end
    end
    if isempty(targetIndex) || isempty(useHalfPeriod)
        switch orbit
            case {'lyapunov', 'dro', 'dpo', 'lpo', 'resonant'}
                targetIndex = [2, 4];
                useHalfPeriod = true;
            case {'shortperiod', 'longperiod', 'dragonfly'}
                targetIndex = 1:6;
                useHalfPeriod = false;
            case {'axial'}
                switch location
                    case {'L1', 'L2', 'L3'}
                        targetIndex = [2, 3, 4];
                        useHalfPeriod = true;
                    otherwise
                        targetIndex = 1:6;
                        useHalfPeriod = false;
                end
            case {'vertical'}
                switch location
                    case {'L1', 'L2', 'L3'}
                        targetIndex = [2, 4, 6];
                        useHalfPeriod = true;
                    otherwise
                        targetIndex = 1:6;
                        useHalfPeriod = false;
                end
            otherwise
                targetIndex = [2, 4, 6];
                useHalfPeriod = true;
        end
    end

    initstates = [];
    periods = [];
    allamps = [];
    allphases = [];
    ss = F.S;
    for ii = 1:numcases
        amp = targetamps(ii);
        phase = targetphases(ii);

        % Interpolate the initial state and period for the requested amplitude.
        f = F.A(ss) - amp;
        idx = find(f(1:end-1) .* f(2:end) <= 0);
        nRoots = numel(idx);
        tmpstates = zeros(nRoots, 6);
        tmpperiods = zeros(nRoots, 1);

        rootCnt = 0;
        for k = idx(:).'
            sL = ss(k);
            sH = ss(k+1);
            if abs(f(k)) < 1e-10
                sR = sL;
            elseif abs(f(k+1)) < 1e-10
                sR = sH;
            else
                sR = fzero(@(s) F.A(s) - amp, [sL, sH]);
            end
            rootCnt = rootCnt + 1;
            for j = 1:6
                tmpstates(rootCnt, j) = F.X{j}(sR);
            end
            tmpperiods(rootCnt) = F.T(sR);
        end

        if usewarnings && nRoots > 1
            warning('crtbp_initstates:MultipleReturns', ['Amplitude %s ' ...
                'returned %d initial states (rows %d to %d). Keep the ' ...
                'state, period, and amplitude mapping in mind.'], num2str(amp), nRoots, ...
                size(initstates, 1)+1, size(initstates, 1) + nRoots);
        end

        % Correct each interpolated initial state.
        for jj = 1:nRoots
            [tmpstates(jj, :), tmpperiods(jj), residual, exitflag] = ...
                orbits.crtbp_diffcorrection( ...
                tmpstates(jj, :), tmpperiods(jj), ...
                freeIndex, targetIndex, ...
                'IsPeriodFree', isPeriodFree, ...
                'UseHalfPeriod', useHalfPeriod, ...
                'Tolerance', tolerance, ...
                'MaxIterations', maxiter, ...
                'OdeOptions', odeOptions);
            if usewarnings && exitflag <= 0
                warning('crtbp_initstates:CorrectionFailure', [ ...
                    'Correction for amplitude %s reached the iteration limit; ' ...
                    'the residual is %.15e.'], ...
                    num2str(amp), residual);
            end
        end

        if phase ~= 0
            for jj = 1:nRoots
                tspan = [0, phase * tmpperiods(jj) / (2*pi)];
                [~, allstates] = ode113(@dynamics.crtbp, tspan, ...
                    tmpstates(jj, :), odeOptions);
                tmpstates(jj, :) = allstates(end, :);
            end
        end

        initstates = [initstates; tmpstates]; %#ok<AGROW>
        periods = [periods; tmpperiods]; %#ok<AGROW>
        allamps = [allamps; amp * ones(nRoots, 1)];  %#ok<AGROW>
        allphases = [allphases; phase * ones(nRoots, 1)];  %#ok<AGROW>
    end

    targetamps = allamps;
    targetphases = allphases;

    % Compute additional orbit metadata.
    allcases = length(periods);
    fullepochs = cell(allcases, 1);
    fullstates = cell(allcases, 1);
    jacobiConstants = zeros(allcases, 1);
    stabilityIndices = zeros(allcases, 1);
    for ii = 1:allcases
        tspan = [0, periods(ii)];
        x0 = [initstates(ii, :)'; reshape(eye(6), 36, 1)];
        [tmpT, tmpX] = ode113(@dynamics.crtbp, tspan, x0, odeOptions);
        jacobiConstants(ii) = computeJacobiConstant(initstates(ii, :), mu);
        stabilityIndices(ii) = computeStabilityIndex(tmpX);
        fullepochs{ii} = tmpT;
        fullstates{ii} = tmpX(:, 1:6);
    end

    % Build the output structure.
    info.initstates = initstates;
    info.periods = periods;
    info.fullepochs = fullepochs;
    info.fullstates = fullstates;
    info.jacobiConstants = jacobiConstants;
    info.stabilityIndices = stabilityIndices;
    info.targetamps = targetamps;
    info.targetphases = targetphases;

    % Plot the generated orbits when requested.
    if useplot
        lp = constants('EMLibrationPoints');
        % planarorbit = ismember(orbit, {'lyapunov', 'shortperiod', ...
        %     'longperiod', 'dro', 'dpo', 'lpo', 'resonant'});
        figure;
        for ii = 1:length(periods)
            label = sprintf('Amplitude %s, initial phase %s', ...
                num2str(targetamps(ii)), num2str(targetphases(ii)));
            plot3(fullstates{ii}(:, 1), fullstates{ii}(:, 2), ...
                fullstates{ii}(:, 3), 'LineWidth', 0.8, ...
                'DisplayName', label); hold on;
            % Mark the initial and final points.
            plot3(fullstates{ii}(1, 1), fullstates{ii}(1, 2), ...
                fullstates{ii}(1, 3), 'r*', 'HandleVisibility', 'off');
            plot3(fullstates{ii}(end, 1), fullstates{ii}(end, 2), ...
                fullstates{ii}(end, 3), 'bo', 'HandleVisibility', 'off');
        end
        grid on;
        % Mark the requested libration point.
        if ~isempty(location)
            switch location
                case 'L1'
                    plot3(lp(1, 1), 0, 0, 'ko', ...
                        'HandleVisibility', 'off');
                    text(lp(1, 1)+1e-3, 0, 0, 'L1', ...
                        'HandleVisibility', 'off');
                case 'L2'
                    plot3(lp(2, 1), 0, 0, 'ko', ...
                        'HandleVisibility', 'off');
                    text(lp(2, 1)+1e-3, 0, 0, 'L2', ...
                        'HandleVisibility', 'off');
                case 'L3'
                    plot3(lp(3, 1), 0, 0, 'ko', ...
                        'HandleVisibility', 'off');
                    text(lp(3, 1)+1e-3, 0, 0, 'L3', ...
                        'HandleVisibility', 'off');
                case 'L4'
                    plot3(lp(4, 1), lp(4, 2), 0, 'ko', ...
                        'HandleVisibility', 'off');
                    text(lp(4, 1)+1e-3, lp(4, 2), 0, 'L4', ...
                        'HandleVisibility', 'off');
                case 'L5'
                    plot3(lp(5, 1), lp(5, 2), 0, 'ko', ...
                        'HandleVisibility', 'off');
                    text(lp(5, 1)+1e-3, lp(5, 2), 0, 'L5', ...
                        'HandleVisibility', 'off');
            end
        end

        % Mark L1, L2, and the Moon for families defined around the Moon.
        if ismember(orbit, {'butterfly', 'dragonfly', 'dro', 'dro3d', 'dpo', ...
                'lpo', 'resonant'})
            plot3(lp(1, 1), 0, 0, 'ko', 'HandleVisibility', 'off');
            text(lp(1, 1)+1e-3, 0, 0, 'L1', 'HandleVisibility', 'off');
            plot3(lp(2, 1), 0, 0, 'ko', 'HandleVisibility', 'off');
            text(lp(2, 1)+1e-3, 0, 0, 'L2', 'HandleVisibility', 'off');
            plot3(1-mu, 0, 0, 'ko', 'HandleVisibility', 'off');
            text(1-mu+1e-3, 0, 0, 'Moon', 'HandleVisibility', 'off');
        end

        % Resonant families also show Earth and the remaining libration points.
        if strcmp(orbit, 'resonant')
            plot3(lp(3, 1), 0, 0, 'ko', 'HandleVisibility', 'off');
            text(lp(3, 1)+1e-3, 0, 0, 'L3', 'HandleVisibility', 'off');
            plot3(lp(4, 1), lp(4, 2), 0, 'ko', 'HandleVisibility', 'off');
            text(lp(4, 1)+1e-3, lp(4, 2)+1e-3, 0, 'L4', ...
                'HandleVisibility', 'off');
            plot3(lp(5, 1), lp(5, 2), 0, 'ko', 'HandleVisibility', 'off');
            text(lp(5, 1)+1e-3, lp(5, 2)+1e-3, 0, 'L5', ...
                'HandleVisibility', 'off');
            plot3(-mu, 0, 0, 'ko', 'HandleVisibility', 'off');
            text(-mu+1e-3, 0, 0, 'Earth', 'HandleVisibility', 'off');
        end

        % Lyapunov and NRHO families show the nearby primary body.
        if ismember(orbit, {'lyapunov', 'nrho'})
            if ismember(location, {'L1', 'L2'})
                plot3(1-mu, 0, 0, 'ko', 'HandleVisibility', 'off');
                text(1-mu+1e-3, 0, 0, 'Moon', 'HandleVisibility', 'off');
            else
                plot3(-mu, 0, 0, 'ko', 'HandleVisibility', 'off');
                text(-mu+1e-3, 0, 0, 'Earth', 'HandleVisibility', 'off');
            end
        end
        
        % Configure the title, legend, and axes.
        title('Earth--Moon CRTBP periodic orbit');
        legend;
        xlabel('X [LU]');
        ylabel('Y [LU]');
        zlabel('Z [LU]');
        daspect([1 1 1]);
        pbaspect([1 1 1]);
    end
end

% Load CRTBP orbit data.
function [X, T, A, S] = crtbp_loadData(orbit, location, direction)
    % Resolve the data directory relative to this function file so the
    % package works from a clean MATLAB session and any current directory.
    filepath = fullfile(fileparts(mfilename('fullpath')), 'crtbp_data');
    if ismember(orbit, {'butterfly', 'dragonfly', 'dro', 'dro3d', 'dpo'})
        filename = [char(orbit), '.mat'];
    elseif strcmp(orbit, 'resonant')
        filename = [char(orbit), '_', ...
            strrep(char(location), ':', ''), '.mat'];
    elseif strcmp(orbit, 'lpo')
        if isempty(direction)
            direction = 'eastern';
        end
        filename = [char(orbit), '_', char(direction), '.mat'];
    else
        filename = [char(orbit), '_', char(location), '.mat'];
    end
    load(fullfile(filepath, filename), 'X', 'T', 'A', 'S');

    % Reflect the z position and z velocity for southern families.
    if ismember(orbit, {'halo', 'nrho', 'axial', 'butterfly', 'dragonfly', ...
            'dro3d'}) && strcmp(direction, 'southern')
        X(:, [3, 6]) = -X(:, [3, 6]);
    end
    % % Southern axial-family z velocity reflection can be enabled here if
    % % the underlying data convention requires it.
    % if strcmp(orbit, 'axial') && strcmp(direction, 'southern')
    %     X(:, 6) = -X(:, 6);
    % end
end

% Validate and unpack options.
function varargout = handleOptions(options)
    if ~islogical(options.UseWarnings)
            error('checkOptions:InvalidOption', ...
                'UseWarnings must be a logical scalar.');
    end
    if ~islogical(options.UsePlot)
        error('checkOptions:InvalidOption', ...
            'UsePlot must be a logical scalar.');
    end
    if ~ismember(options.InterpolationMethod, {'linear', 'makima', ...
            'pchip', 'spline'})
        error('checkOptions:InvalidOption', ['InterpolationMethod must be ' ...
            '''linear'', ''makima'', ''pchip'', or ''spline''.']);
    end
    if ~isempty(options.DiffCorrFreeIndex) && ...
            (~isnumeric(options.DiffCorrFreeIndex) || ...
            ~isvector(options.DiffCorrFreeIndex) || ...
            any(~ismember(options.DiffCorrFreeIndex, 1:6)))
        error('checkOptions:InvalidOption', ...
            'DiffCorrFreeIndex must be a subset of the integers 1 through 6.');
    end
    if ~isempty(options.DiffCorrTargetIndex) && ...
            (~isnumeric(options.DiffCorrTargetIndex) || ...
            ~isvector(options.DiffCorrTargetIndex) || ...
            any(~ismember(options.DiffCorrTargetIndex, 1:6)))
        error('checkOptions:InvalidOption', ...
            'DiffCorrTargetIndex must be a subset of the integers 1 through 6.');
    end
    if ~islogical(options.DiffCorrIsPeriodFree)
        error('checkOptions:InvalidOption', ...
            'DiffCorrIsPeriodFree must be a logical scalar.');
    end
    if ~isempty(options.DiffCorrIsHalfPeriod) && ...
            ~islogical(options.DiffCorrIsHalfPeriod)
        error('checkOptions:InvalidOption', ...
            'DiffCorrIsHalfPeriod must be a logical scalar.');
    end
    if ~isnumeric(options.DiffCorrTolerance) || ...
            ~isscalar(options.DiffCorrTolerance) || ...
            options.DiffCorrTolerance <= 0
        error('checkOptions:InvalidOption', ...
            'DiffCorrTolerance must be a positive scalar.');
    end
    if ~isnumeric(options.DiffCorrMaxIterations) || ...
            ~isscalar(options.DiffCorrMaxIterations) || ...
            options.DiffCorrMaxIterations <= 0
        error('checkOptions:InvalidOption', ...
            'DiffCorrMaxIterations must be a positive scalar.');
    end
    if ~isstruct(options.OdeOptions)
        error('checkOptions:InvalidOption', ...
            'OdeOptions must be a structure.');
    end
    varargout = {
        options.UseWarnings;
        options.UsePlot;
        options.InterpolationMethod;
        options.DiffCorrFreeIndex;
        options.DiffCorrTargetIndex;
        options.DiffCorrIsPeriodFree;
        options.DiffCorrIsHalfPeriod;
        options.DiffCorrTolerance;
        options.DiffCorrMaxIterations;
        options.OdeOptions
    };
end

% Compute the Jacobi constant.
function jacobiConstant = computeJacobiConstant(initstates, mu)
    x = initstates(1);
    y = initstates(2);
    z = initstates(3);
    vx = initstates(4);
    vy = initstates(5);
    vz = initstates(6);
    r1 = sqrt((x + mu)^2 + y^2 + z^2);
    r2 = sqrt((x - 1 + mu)^2 + y^2 + z^2);
    U = (1 - mu) / r1 + mu / r2 + 0.5 * (x^2 + y^2);
    jacobiConstant = 2 * U - (vx^2 + vy^2 + vz^2);
end

% Compute the stability index.
function stabilityIndex = computeStabilityIndex(allstates)
    monodromy = reshape(allstates(end, 7:42), 6, 6);
    maxeig = max(abs(eig(monodromy)));
    stabilityIndex = 0.5 * (maxeig + 1 / maxeig);
end

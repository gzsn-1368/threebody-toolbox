function [Xdot, jac] = ephem(t, X, u, options)
%EPHEM Compute the state derivative and Jacobian in a high-fidelity ephemeris model. The model uses a Moon-centred inertial frame.
%   XDOT = EPHEM(T, X) computes the state derivative XDOT of state X at time T
%   in an Earth-Moon high-fidelity ephemeris model and the Moon-centred J2000
%   inertial frame. By default, the model includes the central accelerations
%   of Earth and the Moon, Earth's J2 perturbation, solar gravitational
%   perturbation, and solar-radiation-pressure perturbation.
%
%   [XDOT, JAC] = EPHEM(T, X) also returns the state Jacobian JAC for
%   linearized dynamics or numerical integration of the state transition
%   matrix.
%
%   [XDOT, JAC] = EPHEM(T, X, U, OPTIONS) specifies the current control input
%   U in km/s^2. If U is not provided, it defaults to a three-element zero
%   vector. OPTIONS specifies the higher-order perturbations included in the
%   model; see DYNAMICS.EPHEMOPTIONS. Supported perturbations include lunar
%   nonspherical gravity, high-order Earth gravity, multi-body gravity,
%   general-relativistic corrections, and Earth albedo.
%
%   When X has 42 elements and contains a state transition matrix, the
%   function also computes the state-transition-matrix derivative and appends
%   its vectorized form to XDOT.
%
%   Coordinate frame and units:
%     - All position, velocity, and acceleration quantities use the
%       Moon-centred J2000 inertial frame;
%     - External-body positions, including Earth and the Sun, are obtained
%       from SPICE SPK files in the Moon-centred J2000 frame;
%     - The first three elements of X are position in km;
%     - Elements 4 through 6 of X are velocity in km/s;
%     - Elements 4 through 6 of XDOT are acceleration in km/s^2.
%
%   Input parameters:
%     T        - Scalar ephemeris time in seconds.
%     X        - Row or column state vector with 6 or 42 elements:
%                  - For 6 elements, X = [x; y; z; dx; dy; dz], where the
%                    first three elements are in km and the last three are
%                    in km/s;
%                  - For 42 elements, the first six elements are as above and
%                    the remaining 36 elements are the column-major form of a
%                    6-by-6 state transition matrix PHI. The function computes
%                    dot(PHI) = JAC * PHI and appends its column-major form to
%                    XDOT.
%     U        - Row or column vector with 3 elements, giving control
%                acceleration along the three coordinate axes in km/s^2.
%                The default is a zero vector.
%     OPTIONS  - Structure specifying the higher-order perturbations included
%                in the model.
%
%   Output parameters:
%     XDOT     - Column state-derivative vector with the same length as X:
%                  - For a 6-element X, XDOT = [dx; dy; dz; ddx; ddy; ddz];
%                  - For a 42-element X, the first six elements are the state
%                    derivative and the remaining 36 elements are the
%                    column-major state-transition-matrix derivative.
%     JAC      - 6-by-6 state Jacobian corresponding to X.
%
%   Note: For integration or shooting solvers that call this function
%   repeatedly, explicitly pass OPTIONS through an anonymous function as
%   @(T, X) DYNAMICS.EPHEM(T, X, U, OPTIONS). Otherwise, every call to this
%   function also calls DYNAMICS.EPHEMOPTIONS, which significantly reduces
%   performance.
%
%   See also DYNAMICS.EPHEMOPTIONS, DYNAMICS.CRTBP.
%
%   Copyright 2025 Guo Zisen.

    configureKernelPath;

    if nargin < 4
        options = dynamics.ephemoptions;
        if nargin < 3
            u = zeros(3, 1);
        end
    end

    % Validate inputs.
    if ~isnumeric(t) || ~isscalar(t)
        error('ephem:InvalidInput', 't must be a numeric scalar.');
    end
    nX = length(X);
    u = u(:);
    if ~ismember(nX, [6, 42])
        error('ephem:InvalidInput', 'X must contain 6 or 42 elements.');
    end
    if isempty(u) || nargin < 3
        u = zeros(3, 1);
    end
    if length(u) ~= 3
        error('ephem:InvalidInput', 'u must contain 3 elements.');
    end

    try
        % Parse the OPTIONS structure.
        earthdeg = options.EarthGravityField.deg;
        earthord = options.EarthGravityField.ord;
        earthGM = options.EarthGravityField.GM;
        earthRe = options.EarthGravityField.Re;
        earthC = options.EarthGravityField.C;
        earthS = options.EarthGravityField.S;

        lunardeg = options.LunarGravityField.deg;
        lunarord = options.LunarGravityField.ord;
        lunarGM = options.LunarGravityField.GM;
        lunarRe = options.LunarGravityField.Re;
        lunarC = options.LunarGravityField.C;
        lunarS = options.LunarGravityField.S;

        fb = options.FourthBody;
        doRel = options.Relativity;
        doSrp = options.SolarRadiationPressure;
        doAlb = options.EarthAlbedo;

        rpcm = options.RPCM;
    catch
        error('ephem:InvalidOption', ...
            'OPTIONS must be a structure returned by dynamics.ephemoptions.');
    end

    % Translational state vector.
    X = X(:);
    p = X(1:3);
    v = X(4:6);

    pEarth = forcemodel('spkpos', 'EARTH', t, 'J2000', 'NONE', 'MOON');
    pSun = forcemodel('spkpos', 'SUN', t, 'J2000', 'NONE', 'MOON');
    doJac = (nX == 42) || (nargout > 1);

    % Earth gravity acceleration.
    earthTM = forcemodel('pxform', 'J2000', 'ITRF93', t);
    [aEarth1, jacEarth] = forcemodel('gravityfield', p-pEarth, ...
        earthdeg, earthord, earthTM(1:3, 1:3), earthGM, earthRe, earthC, ...
        earthS, doJac, false);
    [aEarth2, ~] = forcemodel('gravityfield', -pEarth, ...
        earthdeg, earthord, earthTM(1:3, 1:3), earthGM, earthRe, earthC, ...
        earthS, doJac, false);
    aEarth = aEarth1 - aEarth2;

    % Lunar gravity acceleration.
    lunarTM = forcemodel('pxform', 'J2000', 'IAU_MOON', t);
    [aLuna, jacLuna] = forcemodel('gravityfield', p, lunardeg, ...
        lunarord, lunarTM(1:3, 1:3), lunarGM, lunarRe, lunarC, lunarS, ...
        doJac, false);

    % Fourth-body gravity acceleration.
    [aFourthbody, jacFourthbody] = forcemodel('fourthbody', t, p, fb, doJac);

    % Relativistic acceleration correction.
    if doRel
        [aRel, jacRel] = forcemodel('relativity', p, v, doJac);
    else
        aRel = zeros(3, 1);
        jacRel = zeros(3, 6);
    end

    % Solar-radiation-pressure acceleration.
    switch doSrp
        case 'cannonball'
            [nu, nugrad] = forcemodel('dualcone', p, pSun, pEarth);
            [aSrpRaw, jacSrpRaw] = forcemodel('srp', p, pSun, rpcm, doJac);
            aSrp = nu * aSrpRaw;
            jacSrp = [aSrpRaw * nugrad, zeros(3)] + nu * jacSrpRaw;
        otherwise
            aSrp = zeros(3, 1);
            jacSrp = zeros(3, 6);
    end

    % Earth-albedo acceleration.
    if doAlb
        [aAlb, jacAlb] = forcemodel('earthalbedo', p, pEarth, pSun, ...
            rpcm, doJac);
    else
        aAlb = zeros(3, 1);
        jacAlb = zeros(3, 6);
    end

    % Translational dynamics.
    a = aEarth + aLuna + aFourthbody + aRel + aSrp + aAlb + u;
    Xdot = [v; a];

    % First-order variational equations for the translational dynamics.
    if nX == 42 || nargout > 1
        jac = jacEarth + jacLuna + jacFourthbody + jacRel + jacSrp + jacAlb;
        jacobian = [zeros(3), eye(3); jac];

        % Variational equations and output.
        if nX == 42
            STM = jacobian * reshape(X(7:42), 6, 6);
            Xdot = [Xdot; reshape(STM, 36, 1)];
        end
        if nargout > 1
            jac = jacobian;
        end
    end
end

function configureKernelPath
    persistent isConfigured
    if ~isempty(isConfigured)
        return;
    end

    toolboxRoot = fileparts(fileparts(mfilename('fullpath')));
    kernelPath = fullfile(toolboxRoot, 'kernels');
    requiredKernel = fullfile(kernelPath, 'de440.bsp');
    if ~isfile(requiredKernel)
        error('ephem:MissingKernels', ...
            'Required SPICE kernels were not found in %s.', kernelPath);
    end
    setenv('SPICE_KERNEL_PATH', kernelPath);
    isConfigured = true;
end

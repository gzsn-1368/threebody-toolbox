function varargout = constants(varargin)
%CONSTANTS  Return predefined physical and astrodynamics constants.
%   VAL = CONSTANTS(NAME) returns the constant associated with NAME, for
%   example, gm = CONSTANTS('SunGM').
%
%   [VAL1, VAL2, ...] = CONSTANTS(NAME1, NAME2, ...) returns one value for
%   each requested name in the same order. Unknown names return NaN, for
%   example, [c, mu] = CONSTANTS('LightSpeed', 'MassUnitEM').
%
%   [VAL1, VAL2, ...] = CONSTANTS({NAME1, NAME2, ...}) accepts a cell array
%   of names and has the same behavior.
%
%   S = CONSTANTS(...) returns a structure with one field per requested name
%   when multiple names are supplied with a single output.
%
%   Supported names and units:
%       'SunGM'                 Solar gravitational parameter [km^3/s^2]
%       'MercuryGM'             Mercury gravitational parameter [km^3/s^2]
%       'VenusGM'               Venus gravitational parameter [km^3/s^2]
%       'EarthGM'               Earth gravitational parameter [km^3/s^2]
%       'MoonGM'                Moon gravitational parameter [km^3/s^2]
%       'MarsBarycenterGM'      Mars-system barycenter parameter [km^3/s^2]
%       'JupiterBarycenterGM'   Jupiter-system barycenter parameter [km^3/s^2]
%       'SaturnBarycenterGM'    Saturn-system barycenter parameter [km^3/s^2]
%       'UranusBarycenterGM'    Uranus-system barycenter parameter [km^3/s^2]
%       'NeptuneBarycenterGM'   Neptune-system barycenter parameter [km^3/s^2]
%       'PlutoBarycenterGM'     Pluto-system barycenter parameter [km^3/s^2]
%       'LengthUnitEM'          Earth--Moon distance [m]
%       'TimeUnitEM'            Earth--Moon time unit [s]
%       'MassUnitEM'            Earth--Moon mass parameter [dimensionless]
%       'EarthRadius'           Earth radius [km]
%       'MoonRadius'            Moon radius [km]
%       'LengthUnitSE'          Sun--Earth distance [m]
%       'MassUnitSE'            Sun--Earth mass parameter [dimensionless]
%       'LightSpeed'            Speed of light [km/s]
%       'SolarLuminosity'       Solar luminosity [W]
%       'EarthReflectionRate'   Earth albedo [dimensionless]
%       'EML1'...'EML5'         Earth--Moon libration-point coordinates
%                              [dimensionless]
%       'EML1d'...'EML5d'       Earth--Moon libration-point coordinates [km]
%       'EMLibrationPoints'     All Earth--Moon libration points [dimensionless]
%       'EMLibrationPointsD'    All Earth--Moon libration points [km]
%       'EarthGravity'          Standard Earth gravity [m/s^2]
%       'EccentricityEM'        Earth--Moon orbital eccentricity [dimensionless]
%
%   Copyright 2025 Guo Zisen.

    persistent C
    if isempty(C)
        EML1 = [0.836915132364303; 0; 0];
        EML2 = [1.155682160292340; 0; 0];
        EML3 = [-1.005062645252109; 0; 0];
        EML4 = [0.487849415730060; 0.866025403784439; 0];
        EML5 = [0.487849415730060; -0.866025403784439; 0];
        keys = {
            'SunGM';
            'MercuryGM';
            'VenusGM';
            'EarthGM';
            'MoonGM';
            'MarsBarycenterGM';
            'JupiterBarycenterGM';
            'SaturnBarycenterGM';
            'UranusBarycenterGM';
            'NeptuneBarycenterGM';
            'PlutoBarycenterGM';
            'LengthUnitEM';
            'TimeUnitEM';
            'MassUnitEM';
            'EarthRadius';
            'MoonRadius';
            'LengthUnitSE';
            'MassUnitSE';
            'LightSpeed';
            'SolarLuminosity';
            'EarthReflectionRate';
            'EML1';
            'EML2';
            'EML3';
            'EML4';
            'EML5';
            'EML1d';
            'EML2d';
            'EML3d';
            'EML4d';
            'EML5d';
            'EMLibrationPoints';
            'EMLibrationPointsD';
            'EarthGravity';
            'EccentricityEM';
        };
        vals = {
            1.327124400419393e11;
            2.203178000000002e4;
            3.248585920000000e5;
            3.986004354360959e5;
            4.902800066163796e3;
            4.282837521400002e4;
            1.267127648000002e8;
            3.794058520000000e7;
            5.794548600000008e6;
            6.836527100580023e6;
            9.770000000000007e2;
            3.844e8;
            3.751902619517228e5;
            0.012150584269940;
            6.3781363e3;
            1738;
            1.495978707000000e11;
            3.003480593992993e-06;
            299792.458;
            3.839e26;
            0.3;
            EML1;
            EML2;
            EML3;
            EML4;
            EML5;
            EML1 * 3.844e8 / 1e3;
            EML2 * 3.844e8 / 1e3;
            EML3 * 3.844e8 / 1e3;
            EML4 * 3.844e8 / 1e3;
            EML5 * 3.844e8 / 1e3;
            [EML1, EML2, EML3, EML4, EML5]';
            [EML1, EML2, EML3, EML4, EML5]' * 3.844e8 / 1e3;
            9.80665;
            0.0549006
        };
        C = containers.Map(keys, vals);
    end
    n = nargin;
    if n < 1
        varargout = repmat({0}, nargout, 1);
        return;
    end
    if isscalar(varargin) && iscell(varargin{1})
        varargin = varargin{1};
        n = numel(varargin);
    end
    
    if nargout == 1 && n > 1
        S = struct();
        for ii = 1:n
            key = varargin{ii};
            if ~ischar(key) && ~isstring(key)
                error('constants:InvalidInput', ...
                    'Input %d must be a character vector or string.', ii);
            end
            if C.isKey(key)
                S.(matlab.lang.makeValidName(key)) = C(key);
            else
                warning('constants:UnknownKey', 'Unknown constant name: ''%s''.', key);
                S.(matlab.lang.makeValidName(key)) = nan;
            end
        end
        varargout{1} = S;
    else
        varargout = cell(n, 1);
        for ii = 1:n
            key = varargin{ii};
            if ~ischar(key) && ~isstring(key)
                error('constants:InvalidInput', ...
                    'Input %d must be a character vector or string.', ii);
            end
            if C.isKey(key)
                varargout{ii} = C(key);
            else
                warning('constants:UnknownKey', 'Unknown constant name: ''%s''.', key);
                varargout{ii} = nan;
            end
        end
    end
end

function options = ephemoptions(varargin)
%EPHEMOPTIONS Configure high-fidelity model options.
%   OPTIONS = EPHEMOPTIONS(NAME, VALUE, ...) constructs and returns an options
%   structure for the high-fidelity model from the specified name-value pairs.
%
%   Input parameters:
%   All input parameters are optional name-value pairs:
%     'EarthGravityModel'        String or character vector specifying the
%                                Earth gravity-field model. Valid values are
%                                'EGM2008', 'EGM96', and 'EIGENGL04C'. The
%                                default is 'EGM2008'.
%     'EarthHarmonics'           Two-element vector whose first and second
%                                elements are the degree and order of the
%                                Earth gravity-field expansion. The default
%                                is [2, 0].
%     'LunarGravityModel'        String or character vector specifying the
%                                lunar gravity-field model. Valid values are
%                                'LP100K' and 'LP165P'. The default is
%                                'LP165P'.
%     'LunarHarmonics'           Two-element vector whose first and second
%                                elements are the degree and order of the
%                                lunar gravity-field expansion. The default
%                                is [1, 0].
%     'FourthBody'               String, character vector, or cell array of
%                                them specifying perturbing fourth bodies.
%                                Valid values are:
%                                  'all'                 All fourth bodies
%                                  'sun'                 Sun
%                                  'mercury'             Mercury
%                                  'venus'               Venus
%                                  'mars-barycenter'     Mars barycentre
%                                  'jupiter-barycenter'  Jupiter barycentre
%                                  'saturn-barycenter'   Saturn barycentre
%                                  'uranus-barycenter'   Uranus barycentre
%                                  'neptune-barycenter'  Neptune barycentre
%                                  'pluto-barycenter'    Pluto barycentre
%                                  'none'                No fourth-body term
%                                The default is 'sun'.
%     'Relativity'               Logical value specifying whether to apply a
%                                relativistic correction. The default is
%                                false.
%     'SolarRadiationPressure'   String or character vector specifying the
%                                solar-radiation-pressure model. Valid values
%                                are 'cannonball' and 'none'. The default is
%                                'cannonball'.
%     'EarthAlbedo'              Logical value specifying whether to include
%                                Earth-albedo acceleration. The default is
%                                false.
%     'RPCM'                     Positive scalar combining solar-radiation
%                                pressure and area-to-mass ratio. The default
%                                is 0.005.
%
%   Output parameters:
%     OPTIONS   - Options structure for the high-fidelity model.
%
%   See also DYNAMICS.EPHEM.
%
%   Copyright 2025 Guo Zisen.

    p = inputParser;

    % Earth gravity field.
    addParameter(p, 'EarthGravityModel', 'EGM2008');
    addParameter(p, 'EarthHarmonics', [2, 0]);

    % Lunar gravity field.
    addParameter(p, 'LunarGravityModel', 'LP165P');
    addParameter(p, 'LunarHarmonics', [1, 0]);

    % Fourth-body gravity.
    addParameter(p, 'FourthBody', 'sun');

    % Relativistic correction.
    addParameter(p, 'Relativity', false);

    % Solar radiation pressure.
    addParameter(p, 'SolarRadiationPressure', 'cannonball');
    addParameter(p, 'RPCM', 0.005);

    % Earth albedo.
    addParameter(p, 'EarthAlbedo', false);

    % Parse all inputs.
    parse(p, varargin{:});
    options = p.Results;
    
    % Validate the options structure.
    options = checkEphemOptions(options);

    toolboxRoot = fileparts(fileparts(mfilename('fullpath')));
    kernelPath = fullfile(toolboxRoot, 'kernels');

    % Load gravity-field models.
    loadGravityModel(options.EarthGravityModel);
    loadGravityModel(options.LunarGravityModel);

    % Load a gravity-field model.
    function loadGravityModel(model)
        switch lower(model)
            case 'egm2008'
                % Earth gravity-field model: 1.
                centralbody = 'Earth';
                modelVer = 1;
                load(fullfile(kernelPath, 'egm2008.mat'), ...
                    'GM', 'Re', 'degree', 'C', 'S');
                default = 120;
            case 'egm96'
                % Earth gravity-field model: 2.
                centralbody = 'Earth';
                modelVer = 2;
                load(fullfile(kernelPath, 'egm96.mat'), ...
                    'GM', 'Re', 'degree', 'C', 'S');
                default = 70;
            case 'lp100k'
                % Lunar gravity-field model: 3.
                centralbody = 'Lunar';
                modelVer = 3;
                load(fullfile(kernelPath, 'lp100k.mat'), ...
                    'GM', 'Re', 'degree', 'C', 'S');
                default = 60;
            case 'lp165p'
                % Lunar gravity-field model: 4.
                centralbody = 'Lunar';
                modelVer = 4;
                load(fullfile(kernelPath, 'lp165p.mat'), ...
                    'GM', 'Re', 'degree', 'C', 'S');
                default = 60;
            case 'eigengl04c'
                % Earth gravity-field model: 5.
                centralbody = 'Earth';
                modelVer = 5;
                load(fullfile(kernelPath, 'eigengl04c.mat'), ...
                    'GM', 'Re', 'degree', 'C', 'S');
                default = 70;
        end

        GravityField = [centralbody, 'GravityField'];
        Harmonics = [centralbody, 'Harmonics'];

        options.(GravityField) = struct( ...
            'modelVer', modelVer, ...
            'deg', options.(Harmonics)(1), ...
            'ord', options.(Harmonics)(2), ...
            'GM', GM, ...
            'Re', Re, ...
            'degree', degree, ...
            'C', C, ...
            'S', S, ...
            'default', default);

        if options.(GravityField).deg > options.(GravityField).degree
            warning('ephemoptions:MaxDegExceeded', ...
                'Requested degree exceeds the model limit; using %d.', ...
                options.(GravityField).default);
            options.(GravityField).deg = options.(GravityField).default;
        end

        options = rmfield(options, {[centralbody, 'GravityModel'], ...
            [centralbody, 'Harmonics']});
    end
end

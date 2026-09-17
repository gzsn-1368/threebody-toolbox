function options = checkEphemOptions(options)
%CHECKEPHEMOPTIONS Validate and normalize ephemeris-model options.
%   OPTIONS = CHECKEPHEMOPTIONS(OPTIONS) validates fields used by
%   DYNAMICS.EPHEM and returns normalized values where appropriate.
%
%   See also DYNAMICS.EPHEMOPTIONS, DYNAMICS.EPHEM.
%
%   Copyright 2025 Guo Zisen.

    if ~isstruct(options)
        error('checkEphemOptions:NotAStruct', 'options must be a structure.');
    end

    % Option fields.
    validFields = {'EarthGravityModel', 'EarthHarmonics', ...
        'LunarGravityModel', 'LunarHarmonics', 'FourthBody', ...
        'Relativity', 'SolarRadiationPressure', 'EarthAlbedo', 'RPCM'};

    optionFields = fieldnames(options);

    % Check for missing fields.
    if ~all(ismember(validFields, optionFields))
        missing = validFields(~ismember(validFields, optionFields));
        error('checkEphemOptions:MissingFields', ...
            'options is missing required field(s): %s.', ...
            strjoin(missing, ', '));
    end

    % Check for unsupported fields.
    for ii = 1:numel(optionFields)
        field = optionFields{ii};
        if ~ismember(field, validFields)
            error('checkEphemOptions:InvalidField', ...
                '%s is not a valid option for the high-fidelity model.', field);
        end
    end

    % Validate individual parameters.
    if ~ismember(lower(options.EarthGravityModel), ...
            {'egm2008', 'egm96', 'eigengl04c'})
        error('checkEphemOptions:InvalidParam', ...
            ['EarthGravityModel must be ''EGM2008'', ''EGM96'', ' ...
            'or ''EIGENGL04C''.']);
    end
    if ~isnumeric(options.EarthHarmonics) || ...
            ~isvector(options.EarthHarmonics) || ...
            length(options.EarthHarmonics) ~= 2
        error('checkEphemOptions:InvalidParam', ...
            ['EarthHarmonics must be a two-element vector containing ' ...
            'the spherical-harmonic degree and order.']);
    end
    if options.EarthHarmonics(1) < 0 || options.EarthHarmonics(2) < 0
        error('checkEphemOptions:InvalidParam', ...
            'EarthHarmonics elements must be nonnegative.');
    end
    if options.EarthHarmonics(1) < options.EarthHarmonics(2)
        error('checkEphemOptions:InvalidParam', ...
            'EarthHarmonics order cannot exceed its degree.');
    end

    if ~ismember(lower(options.LunarGravityModel), {'lp100k', 'lp165p'})
        error('checkEphemOptions:InvalidParam', ...
            'LunarGravityModel must be ''LP100K'' or ''LP165P''.');
    end
    if ~isnumeric(options.LunarHarmonics) || ...
            ~isvector(options.LunarHarmonics) || ...
            length(options.LunarHarmonics) ~= 2
        error('checkEphemOptions:InvalidParam', ...
            ['LunarHarmonics must be a two-element vector containing ' ...
            'the spherical-harmonic degree and order.']);
    end
    if options.LunarHarmonics(1) < 0 || options.LunarHarmonics(2) < 0
        error('checkEphemOptions:InvalidParam', ...
            'LunarHarmonics elements must be nonnegative.');
    end
    if options.LunarHarmonics(1) < options.LunarHarmonics(2)
        error('checkEphemOptions:InvalidParam', ...
            'LunarHarmonics order cannot exceed its degree.');
    end

    fb = options.FourthBody;
    fbList = {'all', 'sun', 'mercury', 'venus', 'mars-barycenter', ...
        'jupiter-barycenter', 'saturn-barycenter', 'uranus-barycenter', ...
        'neptune-barycenter', 'pluto-barycenter', 'none'};
    if ischar(fb) || isstring(fb)
        fb = {char(fb)};
    end
    if ~iscell(fb)
        error('checkEphemOptions:InvalidParam', ...
            ['FourthBody must be a character vector, string, or cell ' ...
            'array of character vectors or strings.']);
    else
        fb = cellfun(@char, fb, 'UniformOutput', false);
        fb = lower(fb);
    end
    if any(ismember(fb, {'all', 'none'})) && numel(fb) > 1
        error('checkEphemOptions:InvalidParam', ...
            'FourthBody values ''all'' and ''none'' cannot be combined with other values.');
    end
    for ii = 1:numel(fb)
        if ~ismember(lower(fb{ii}), fbList)
            error('checkEphemOptions:InvalidParam', ...
                ['FourthBody value %d is invalid. See the function help ' ...
                'for the supported body names.'], ii);
        end
    end
    if ismember('all', fb)
        fb = true(1, 9);
    elseif ismember('none', fb)
        fb = false(1, 9);
    else
        fb = ismember(fbList(2:end-1), fb);
    end
    options.FourthBody = fb;
    if ~islogical(options.Relativity)
        error('checkEphemOptions:InvalidParam', ...
            'Relativity must be logical true or false.');
    end
    if ~ismember(options.SolarRadiationPressure, {'cannonball', 'none'})
        error('checkEphemOptions:InvalidParam', ...
            'SolarRadiationPressure must be ''cannonball'' or ''none''.');
    end
    if ~islogical(options.EarthAlbedo)
        error('checkEphemOptions:InvalidParam', ...
            'EarthAlbedo must be logical true or false.');
    end
    if strcmpi(options.SolarRadiationPressure, 'cannonball') || options.EarthAlbedo
        if ~isnumeric(options.RPCM) || ...
                ~isscalar(options.RPCM) || ...
                options.RPCM <= 0
            error('checkEphemOptions:InvalidParam', ...
                'RPCM must be a positive scalar.');
        end
    end
end 

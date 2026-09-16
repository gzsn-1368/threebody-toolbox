function checkOrbitValidity(orbit, location, direction)
%CHECKORBITVALIDITY Validate an orbit-family, location, and direction.
%   CHECKORBITVALIDITY(ORBIT, LOCATION, DIRECTION) checks whether the input
%   combination is supported. Invalid combinations raise an error; missing
%   direction values for certain families generate a warning.
%
%   Supported combinations of ORBIT, LOCATION, and DIRECTION are:
%
%   ORBIT          LOCATION/RESONANCE                 DIRECTION
%   'lyapunov'     'L1', 'L2', 'L3'                   []
%   'halo'         'L1', 'L2', 'L3'                   'northern', 'southern'
%   'nrho'         'L1', 'L2', 'L3'                   'northern', 'southern'
%   'vertical'     'L1', 'L2', 'L3', 'L4', 'L5'       []
%   'axial'        'L1', 'L2', 'L3', 'L4', 'L5'       'northern', 'southern'
%   'shortperiod'  'L4', 'L5'                         []
%   'longperiod'   'L4', 'L5'                         []
%   'butterfly'    []                                 'northern', 'southern'
%   'dragonfly'    []                                 'northern', 'southern'
%   'dro'          []                                 'northern', 'southern'
%   'dro3d'        []                                 'northern', 'southern'
%   'dpo'          []                                 []
%   'lpo'          []                                 'eastern', 'western'
%   'resonant'     '1:1', '1:2', '1:3', '1:4', ...    []
%                 '2:1', '2:3', '3:1', '3:2', ...
%                 '3:4', '4:1', '4:3'
%
%   Copyright 2025 Guo Zisen.

    orbitList = {'lyapunov', 'halo', 'nrho', 'vertical', 'axial', ...
        'shortperiod', 'longperiod', 'butterfly', 'dragonfly', ...
        'dro', 'dro3d', 'dpo', 'lpo', 'resonant'};
    % Validate the orbit family.
    if ~ismember(orbit, orbitList)
        error('checkOrbitValidity:InvalidInput', ['Invalid orbit family.\n' ...
            'Expected ''lyapunov'', ''halo'', ''nrho'', ''vertical'', ' ...
            '''axial'', ''shortperiod'', ''longperiod'', ''butterfly'', ' ...
            '''dragonfly'', ''dro'', ''dro3d'', ''dpo'', ''lpo'', or ' ...
            '''resonant''.']);
    end
    % Validate the location or resonance ratio.
    if ismember(orbit, {'shortperiod', 'longperiod'})
        if ~ismember(location, {'L4', 'L5'})
            error('checkOrbitValidity:InvalidInput', ['Invalid location.\n' ...
                'Expected triangular libration point ''L4'' or ''L5''.']);
        end
    elseif ismember(orbit, orbitList(1:3))
        if ~ismember(location, {'L1', 'L2', 'L3'})
            error('checkOrbitValidity:InvalidInput', ['Invalid location.\n' ...
                'Expected collinear libration point ''L1'', ''L2'', or ''L3''.']);
        end
    elseif ismember(orbit, orbitList(4:5))
        if ~ismember(location, {'L1', 'L2', 'L3', 'L4', 'L5'})
            error('checkOrbitValidity:InvalidInput', ['Invalid location.\n' ...
                'Expected ''L1'', ''L2'', ''L3'', ''L4'', or ''L5''.']);
        end
    elseif strcmp(orbit, 'resonant')
        if ~ismember(location, {'1:1', '1:2', '1:3', '1:4', '2:1', ...
                '2:3', '3:1', '3:2', '3:4', '4:1', '4:3'})
            error('checkOrbitValidity:InvalidInput', ['Invalid resonance ratio.\n' ...
                'Expected ''1:1'', ''1:2'', ''1:3'', ''1:4'', ' ...
                '''2:1'', ''2:3'', ''3:1'', ''3:2'', ''3:4'', ''4:1'', or ' ...
                '''4:3''.']);
        end
    end
    % Validate the orbit direction.
    if ismember(orbit, {'halo', 'nrho', 'axial', 'butterfly', 'dragonfly', 'dro3d'})
        if isempty(direction)
            warning('checkOrbitValidity:NoDirection', ...
                'No direction specified; using ''northern''.');
        elseif ~ismember(direction, {'northern', 'southern'})
            error('checkOrbitValidity:InvalidInput', ['Invalid direction.\n' ...
                'Expected ''northern'' or ''southern''.']);
        end
    elseif strcmp(orbit, 'lpo')
        if isempty(direction)
            warning('checkOrbitValidity:NoDirection', ...
                'No direction specified; using ''eastern''.');
        elseif ~ismember(direction, {'eastern', 'western'})
            error('checkOrbitValidity:InvalidInput', ['Invalid direction.\n' ...
                'Expected ''eastern'' or ''western''.']);
        end
    end
end

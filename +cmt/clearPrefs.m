function clearPrefs
% Remove appdata prefs storage.
%
% Needed since "clear classes" complains if any objects are stored
% in appdata.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

try
    rmappdata(0, 'cmt_prefs')
catch err
    if strcmp(err.identifier, ...
            'MATLAB:HandleGraphics:Appdata:InvalidNameAppdata')
        % Nothing there, that's ok.
        return
    else
        rethrow(err)
    end
end

end

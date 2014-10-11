function dispPrefs
% Preference dump to console.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

prefs = getappdata(0, 'cmt_prefs');
if isempty(prefs)
    fprintf('No CMT preferences found.\n')
    return
end

fn = fieldnames(prefs);
for k = 1:numel(fn)
    fprintf(['\n==========================\n' ...
        ' For class %s:\n' ...
        '--------------------------\n'], fn{k})
    disp(prefs.(fn{k}))
end

end

function varargout = optargs(default,arg)

% You give a cell array of default values. Nonempty values in the second
% cell array override the defaults, and the result is dealt to output.

varargout = default;
idx = find( ~cellfun('isempty',arg) );
varargout(idx) = arg(idx);


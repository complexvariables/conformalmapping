classdef optset
% OPTSET base abstract class for options structures.
%
% Provides a base set of functions for providing options structures.
% Subclasses must provide a protected property called 'proplist', along
% with public properties matching the names in the property list.
%
% opt = optset(varargin)
%   Constructor does name/value pair matching/setting based on proplist
%   defined in the subclass, including defaults.
%
% proplist (abstract property)
%   Property list is a n-by-4 cell array, where each nth property has the
%   entry
%   	{ 'name', default_value, function_handle, '[ {default} | value ]' }
%   The function handle is a validator function for the property, which
%   should return true if the value is valid for the property, and false
%   otherwise. This can be [], which means input will not be validated.
%
% Methods:
%   opt = set(opt, 'pref1', value1, 'pref2', value2, ...)
%     Sets option properties using name/value pairs. Since optset subclasses
%     are value classes, the action of the set must be captured as function
%     output.
%
%   value = get(opt, 'pref')
%     Gets the value of 'pref' from optset object.
%
%   args = varargs(opt)
%     Returns a cell array of name/value pairs of preferences from the object.
%
% An example is given in SZSET.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties(Abstract, Access=protected)
    proplist
end

methods
    function opt = optset(varargin)
        % Assign input name-value pairs.
        
        opt = set(fillDefaults(opt), varargin{:});
    end

    function disp(opt)
        % Display values for option class.
        
        nump = size(opt.proplist, 1);
        maxstrlen = 0;
        for k = 1:nump
            maxstrlen = max(maxstrlen, length(opt.proplist{k,1}));
        end

        fprintf('Values for %s object, with format\n', class(opt))
        sep(1:maxstrlen) = '-';
        fprintf('  % *s : [ description {default} ] : current value\n%s\n', ...
            maxstrlen, 'name', ['  ' sep]);
        for k = 1:nump
            curval = evalc(sprintf('disp(opt.%s)', opt.proplist{k,1}));
            fprintf('  % *s : %s : %s\n', ...
                maxstrlen, opt.proplist{k,[1, 4]}, strtrim(curval));
        end
        fprintf('\n')
    end
    
    function value = get(opt, name)
       value = opt.(name); 
    end
    
    function [args, exargs] = separateArgs(opt, varargin)
        % Separatate 'name'/value pairs which belong to optset object and
        % those that don't.
        
        n = numel(varargin);
        if mod(n, 2)
            error('CMT:InvalidArgument', 'Expected name/value pairs.')
        end
        
        recognized = properties(opt);
        args = {};
        idx = 1:n;
        for k = 1:2:n-1
            match = strcmpi(varargin{k}, recognized);
            if ~any(match)
                continue
            end
            args(numel(args)+(1:2)) = ...
                {recognized{match}, varargin{k+1}}; %#ok<AGROW>
            idx = idx(idx ~= k & idx ~= k+1);
        end
        if isempty(idx)
            exargs = {};
        else
            exargs = varargin(idx);
        end

    end

    function [opt, exargs] = set(opt, varargin)
        % Set properties based on 'name'/value pairs.
        %
        % opt = set(opt, varargin) errors on unknown properties.
        % [opt, exargs] = set(opt, varargin) puts unknown properties in the
        % extargs cell array without generating an error.
        
        namelist = opt.proplist(:,1);
        exargs = {};
        for k = 1:2:numel(varargin)
            j = find(strcmp(varargin{k}, namelist));
            if isempty(j)
                if nargout > 1
                    exargs(numel(exargs)+(1:2)) = varargin(k+(0:1)); %#ok<AGROW>
                    continue
                else
                    error('CMT:InvalidArgument', ...
                        'Preference name "%s" is not valid for class %s.', ...
                        varargin{k}, class(opt))
                end
            end
            pname = namelist{j};
            isvalid = opt.proplist{j,3};
            if ~isempty(isvalid) && isa(isvalid, 'function_handle')
                if ~isvalid(varargin{k+1})
                    error('CMT:InvalidArgument', ...
                        ['Value for property "%s" is invalid. ' ...
                         'It must satisfy the function\n\t%s'], ...
                        pname, strtrim(evalc('disp(isvalid)')))
                end
            end
            opt.(pname) = varargin{k+1};
        end
    end
    
    function args = varargs(opt)
        % Output cell array of preferences suitable to pass as name/value
        % pairs in a function call.
        
        plist = properties(opt);
        args = cell(1, 2*numel(plist));
        for k = 1:numel(plist)
            args(2*(k-1)+(1:2)) = {plist{k}, opt.(plist{k})};
        end
    end
end

methods(Access=protected)
    function opt = fillDefaults(opt)
        % Fill structure with defined default values.
        
        names = opt.proplist(:,1);
        values = opt.proplist(:,2);
        for k = 1:numel(names)
            opt.(names{k}) = values{k};
        end
    end
end
    
end

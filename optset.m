classdef optset
% OPTSET base abstract class for options structures.
%
% Provides a base set of functions for providing options structures.
% Subclasses must provide a protected property called 'proplist', along
% with public properties matching the names in the property list.
%
% opt = optset(varargin)
%   Constructor does name/value pair matching/setting based on proplist in
%   subclass via the builtin inputParser class.
%
% defaults(optset)
%   Uses sublass proplist to show settable properties with default values.
%
% proplist
%   Property list is a n-by-4 cell array, where each nth property has the
%   entry
%   	{ 'name', default_value, function_handle, '[ {default} | value ]' }
%   The function handle is a validator function for the property, which is
%   passed to inputParse.addParameter. This can be [], which means input
%   will not be validated on construction.
%
% See also inputParser.
% An example is given in szset.

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
        
        input = nvpair(opt, varargin{:});
        for k = 1:size(opt.proplist, 1)
            fname = opt.proplist{k,1};
            opt.(fname) = input.(fname);
        end
    end

    function defaults(opt)
        % Display default values for option class.
        
        nump = size(opt.proplist, 1);
        maxstrlen = 0;
        for k = 1:nump
            maxstrlen = max(maxstrlen, length(opt.proplist{k,1}));
        end

        fprintf('Values for %s object, with format\n', class(opt))
        fprintf('  % *s : [ description {default} ] : current value\n\n', ...
            maxstrlen, 'name');
        for k = 1:nump
            curval = evalc(sprintf('disp(opt.%s)', opt.proplist{k,1}));
            fprintf('  % *s : %s : %s\n', ...
                maxstrlen, opt.proplist{k,[1, 4]}, strtrim(curval));
        end
        fprintf('\n')
    end

    function disp(opt)
        defaults(opt)
    end

    function input = nvpair(opt, varargin)
        % Parse name-value pairs. See proplist.
        
        p = inputParser;
        for k = 1:size(opt.proplist, 1)
            if ~isempty(opt.proplist{k,3}) ...
                    && isa(opt.proplist{k,3}, 'function_handle')
                ix = 1:3;
            else
                ix = 1:2;
            end
            addParameter(p, opt.proplist{k,ix});
        end
        parse(p, varargin{:});
        input = p.Results;
    end
end
    
end

classdef cmtbase
% CMTBASE provides base CMT functionality for classes.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

methods
    function prefs = get(this, property)
        % Retrieve preference values for module (class) stored in app data
        % facility.
        
        prefs = getappdata(0, 'cmt_prefs');
        if isempty(prefs)
            % Set defaults.
            set(this); 
            prefs = getappdata(0, 'cmt_prefs');
        end
        
        module = class(this);
        if isfield(prefs, module)
            prefs = prefs.(module);
        else
            error('CMT:NotDefined', ...
                'Unable to retreive app data for class "%s".', ...
                module)
        end
        
        if nargin > 1
           if isfield(prefs, property)
               prefs = prefs.(property);
           else
               error('CMT:NotDefined', ...
                   'Property "%s" not defined for class "%s%".', ...
                   property, module)
           end
        end
    end
    
    function set(this, varargin)
        % Store preferences for module (class) in app data facility.
        
        prefs = getappdata(0, 'cmt_prefs');        
        module = class(this);
        if (nargin < 2 || ~isfield(prefs, module)) ...
                && isprop(this, 'optsClass')
            % Set default.
            prefs.(module) = eval(this.optsClass);
        end
        
        for k = 1:2:numel(varargin)
            prefs.(module).(varargin{k}) = varargin{k+1};
        end
        setappdata(0, 'cmt_prefs', prefs);
    end
end

end

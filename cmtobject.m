classdef cmtobject
% CMTOBJECT provides some base CMT functionality for classes.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

methods
    function prefs = get(this, property)
        % Retrieve preference values for module (class) stored in app data
        % facility.
        
        prefs = getappdata(0, 'cmt_prefs');
        module = class(this);

        if isempty(prefs)
            % Set defaults.
            if isa(property, 'optset')
                set(this, property);
                prefs = getappdata(0, 'cmt_prefs');
            else
                error('CMT:NotDefined', ...
                    'No default values given to set for class %s.', ...
                    module)
            end
        end
        
        if isfield(prefs, module)
            prefs = prefs.(module);
        else
            % This shouldn't actually happen.
            error('CMT:NotDefined', ...
                'Unable to retreive app data for class %s.', ...
                module)
        end
        
        if nargin > 1 && ischar(property)
           if isfield(prefs, property)
               prefs = prefs.(property);
           else
               error('CMT:NotDefined', ...
                   'Property %s not defined for class %s%.', ...
                   property, module)
           end
        end
    end
    
    function set(this, varargin)
        % Store preferences for module (class) in app data facility.
        
        prefs = getappdata(0, 'cmt_prefs');        
        module = class(this);
        
        if numel(varargin) == 1
            if isa(varargin{1}, 'optset')
                % Set default.
                prefs.(module) = varargin{1};
                varargin = {};
            else
                error('CMT:InvalidArgument', ...
                    'Expected an optset object.')
            end
        end
        
        for k = 1:2:numel(varargin)
            prefs.(module).(varargin{k}) = varargin{k+1};
        end
        
        setappdata(0, 'cmt_prefs', prefs);
    end
end

end

classdef closedcurve < curve
    % CLOSEDCURVE abstract base class for simple planar Jordan curves.
    
    % This file is a part of the CMToolkit.
    % It is licensed under the BSD 3-clause license.
    % (See LICENSE.)
    
    % Copyright Toby Driscoll, 2014.
    % (Re)written by Everett Kropf, 2014,
    % adapted from Toby Driscoll's code, originally 20??.
    
    properties
        % nothing additional
    end
    
    methods
        
        function c = closedcurve(varargin)
            if isa(varargin{1},'curve')
                a = varargin{1};
                varargin = {a.position,a.tangent,a.bounds};
            end
            c = c@curve(varargin{:});
        end
        
        function img = apply(src,f)
            img = closedcurve( apply@curve(src,f) );
        end

        function t = mod(C, t)
            % Returns equivalent parameter value in the allowed range.
            t = t - C.bounds(1);
            t = C.bounds(1) + mod(t-C.bounds(1),diff(C.bounds));
        end
        
    end
    
    
end

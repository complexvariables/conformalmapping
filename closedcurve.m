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
        
        function r = between(p,q)
            % BETWEEN(P,Q)
            % Return the region between the two closedcurve objects P and Q.
            % They can be specified in either order (outer/inner or
            % inner/outer). No strict checking of nonintersection or nesting is done.
            
            pointOnQ = point(q,0.123456789);
            pointOnP = point(p,0.123456789);
            if isinside(p,pointOnQ)
                % P is outer, Q is inner
                r = region(p,q);
            elseif isinside(q,pointOnP)
                % Q is outer, P is inner
                r = region(q,p);
            else
                error('Boundaries do not appear to be nested')
            end
        end
        
        function str = char(c)
            str = sprintf('closed curve parameterized over [%g,%g]',c.bounds);
        end
        
        function r = exterior(c)
            r = region(c,'exterior');
        end
        
        function r = interior(c)
            r = region(c,'interior');
        end
        
        function t = isinside(c,z)
            %ISINSIDE True for points inside a closed curve.
            %  ISINSIDE(C,Z) returns an array the same size as Z. Each
            %  element is true if and only if the corresponding Z value
            %  lies inside the CLOSEDCURVE C.
            
            b = c.bounds;
            cdisc = point(c,linspace(b(1),b(2),300));
            t = inpolygon(real(z),imag(z),real(cdisc),imag(cdisc));
        end
        
        function t = mod(C, t)
            % Returns equivalent parameter value in the allowed range.
            t = t - C.bounds(1);
            t = C.bounds(1) + mod(t-C.bounds(1),diff(C.bounds));
        end
        
    end
    
    
end

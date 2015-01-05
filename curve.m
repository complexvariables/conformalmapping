classdef (InferiorClasses = {?double}) curve
    % CLOSEDCURVE abstract base class for simple planar Jordan curves.
    
    % This file is a part of the CMToolkit.
    % It is licensed under the BSD 3-clause license.
    % (See LICENSE.)
    
    % Copyright Toby Driscoll, 2014.
    % (Re)written by Everett Kropf, 2014,
    % adapted from Toby Driscoll's code, originally 20??.
    
    properties
        cornerList       % For future development
        position        % complex-valued function returning point for each parameter
        tangent          % function returning tangent vector at each parameter
        bounds = [0 1]     % Parameter bounds
    end
    
    methods
        
        function c = curve(positionfun,tangentfun,bounds)
            c.position = positionfun;
            c.tangent = tangentfun;
            if nargin > 2
                c.bounds = bounds;
            end
        end
        
        function img = apply(src,f)
            posn = @(t) f(src.position(t));
            tant = [];   % TODO
            img = curve(posn,tant,src.bounds);
        end
        
        function box = boundbox(C)
            % Return bounding box containing the curve.
            t = linspace(C.bounds(1),C.bounds(2),300);
            box = cmt.boundbox(point(C, t));
        end
        
        function str = char(c)
            str = sprintf('curve parameterized over [%g,%g]',c.bounds);
        end
        
        function varargout = corner(C, varargin)
            % I'm honestly not quite sure about the usage of this, it will probably
            % change drastically in the future. Haven't found it being used yet in any
            % of Toby's code examples. -- EK
            % Other than being set in the polygon constructor. -- EK
            v = C.cornerList;
            
            % Classes of input args?
            inclass = cellfun(@class, varargin, 'uniformoutput', false);
            
            % Did we get an index?
            j = find(strcmp('double', inclass));
            if ~isempty(j)
                v = v(varargin{j});
            end
            
            % Field selection?
            j = find(strcmp('char', inclass));
            if ~isempty(j)
                varargout{1} = cat(1, v.(varargin{j}));
            else
                if nargout == 3
                    t = cat(1, v.param);
                    z = cat(1, v.point);
                    alpha = cat(1, v.alpha);
                    varargout = {t, z, alpha};
                else
                    varargout{1} = v;
                end
            end
        end
        
        function disp(c)
            fprintf(char(c))
            fprintf('\n\n')
        end
        
        function cm = uminus(c)
            cm = c;
            cm.position = @(t) -c.position(t);
            cm.tangent = @(t) -c.tangent(t);
        end
        
        function out = plot(C, varargin)
            % Plot curve in the plane.
            washold = ishold;
            newplot
            
            h = plotCurve(C);
            [cargs, pargs] = cmtplot.closedcurveArgs(varargin{:});
            set(h, pargs{:}, cargs{:});
            
            if ~washold
                axis(plotbox(C, 1.1));
                set(gca, 'dataaspectratio', [1 1 1])
                hold off
            end
            
            if nargout > 0
                out = h;
            end
        end
        
        function box = plotbox(C, scale)
            % Return plot box for curve using evenly spaced points.
            %
            % See also plotbox.
            
            if nargin < 2
                scale = [];
            end
            t = linspace(C.bounds(1),C.bounds(2),300);
            box = cmt.plotbox(point(C, t), scale);
        end
        
        function c  = minus(c,z)
            c = plus(c,-z);
        end
        
        function c  = plus(c,z)
            if isa(z,'curve')
                tmp = z;
                z = c;
                c = tmp;
            end
            validateattributes(z,{'double'},{'scalar','finite'},...
                'mtimes','scalar');
            c.position = @(t) c.position(t) + z;
        end
        
        function c = mtimes(c,z)
            if isa(z,'curve')
                tmp = z;
                z = c;
                c = tmp;
            end
            validateattributes(z,{'double'},{'scalar','finite'},...
                'mtimes','scalar');
            c.position = @(t) c.position(t)*z;
            c.tangent = @(t) c.tangent(t)*z;
        end
        
        function z = point(c,t)
            z = nan(size(t));
            mask = (t>=c.bounds(1)) & (t<=c.bounds(2));
            z(mask) = c.position(t(mask));
        end
        
        function out = rsplot(C, varargin)
            % Plot curve on the Riemann sphere.
            washold = ishold;
            
            % Draw Riemann shpere if not there.
            if isempty(findobj(gca, 'tag','CMT:RiemannSphere')) || ~washold
                [xs, ys, zs] = sphere(36);
                mesh(0.995*xs, 0.995*ys, 0.995*zs, 'edgecolor', .85*[1 1 1], ...
                    'tag', 'CMT:RiemannSphere')
                hold on
            end
            
            % Draw on the sphere.
            function x = rspoint(t)
                z = point(C, t);
                [x1, x2, x3] = c2rs(z);
                x = [x1, x2, x3];
            end
            h = adaptplot(@rspoint, C.bounds);
            set(h, varargin{:});
            
            if ~washold
                hold off
            end
            axis equal
            
            if nargout > 0
                out = h;
            end
        end
        
        function z = subsref(C, S)
            % Equate C(t) with point(C, t). Fallback to builtin subsref otherwise.
            if numel(S) == 1 && strcmp(S.type, '()')
                if numel(S.subs) == 1
                    z = point(C, S.subs{1});
                    return
                else
                    error('Object only takes single parenthised subscript.')
                end
            end
            
            z = builtin('subsref', C, S);
        end
        
        function xy = xypoint(C, t)
            z = point(C, t);
            xy = [real(z), imag(z)];
        end
    end
    
    methods(Hidden)
        function h = plotCurve(C, varargin)
            h = adaptplot(@(t) xypoint(C, t), C.bounds);
        end
    end
    
    methods (Access=private)
        function [x1, x2, x3] = c2rs(z)
            % C2RS Cartesian complex coordinate to Riemann sphere projection.
            
            % This file is a part of the CMToolbox.
            % It is licensed under the BSD 3-clause license.
            % (See LICENSE.)
            
            % Copyright Toby Driscoll, 2014.
            
            theta = angle(z);
            absz = abs(z);
            phi = atan2(absz.^2 - 1, 2*abs(z));
            phi(isinf(z)) = pi/2;
            [x1, x2, x3] = sph2cart(theta, phi, ones(size(theta)));
            
        end
    end
    
end

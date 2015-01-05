classdef circle < closedcurve
    % CIRCLE is a generalized circle class.
    %
    % C = circle(center, radius)
    %   Creates a circle with given center and radius.
    %
    % C = circle([z1, z2, z3])
    %   Creates a generalized circle passing through the three given points.
    
    % This file is a part of the CMToolbox.
    % It is licensed under the BSD 3-clause license.
    % (See LICENSE.)
    
    % Copyright Toby Driscoll, 2014.
    % (Re)written by Everett Kropf, 2014,
    % adapted from code by Toby Driscoll, originally 20??.
    
    properties
        center      % will be an interior point for a line
        radius
    end
    
    methods
        function c = circle(varargin)
            if nargin==0
                % Return an abstract circle with unknown geometry.
                return
            end
            
            switch nargin
                case 1
                    z3 = varargin{1};
                    if isa(z3, 'double') && numel(z3) == 3
                        
                        % Deduce center and radius.
                        % Use standardmap directly to avoid vicious circularity, since mobius
                        % constructs a circle when given two 3-tuples.
                        M = mobius(mobius.standardmap(z3)\mobius.standardmap([1, 1i, -1]));
                        zi = pole(M);
                        if abs(abs(zi) - 1) < 10*eps
                            % Infinite radius case.
                            center_ = nan;
                            radius_ = inf;
                            % If Inf was given, make it the last point.
                            if isreal(z3)
                                z3 = complex(z3);
                            end
                            z3 = sort(z3);
                        else
                            % Inverse of zi maps to center.
                            center_ = M(1/conj(zi));
                            radius_ = abs(z3(1) - center_);
                        end
                        
                        % Find a point in the interior of the curve. For a line, this is a
                        % point to the "left" as seen by following the given points.
                        if isinf(radius_)
                            tangent = diff(z3(1:2));
                            center_ = z3(1) + 1i*tangent;
                        end
                    end
                    
                case 2
                    [center_, radius_] = varargin{:};
                    validateattributes(center_,{'double'},...
                        {'scalar','finite','nonnan'},...
                        'circle','center');
                    validateattributes(radius_,{'double'},...
                        {'scalar','finite','nonnan','positive'},...
                        'circle','radius');
            end
            
            function z = positionfun(t)
                if ~isinf(radius_)
                    z = center_ + radius_*exp(1i*t);
                else
                    % Use homogeneous coordinates to define a reasonable interpolant.
%                     tangent = diff(c.points(1:2)); % must be finite
%                     upper = 2*tangent*(t - 1/2);
%                     lower = 4*t.*(1 - t);
%                     z = double(homog(c.points(1)) + homog(upper, lower));
                end
            end
            
            function zt = tangentfun(t)
                zt = 1i*exp(1i*t);
            end
            
            c = c@closedcurve(@positionfun,@tangentfun,[0 2*pi]);
            c.center = center_;
            c.radius = radius_;

        end
        
%         function gc = apply(gc, m)
%             if ~isa(m, 'mobius')
%                 error('CMT:NotDefined', ...
%                     'Expected a mobius transformation.')
%             end
%             
%             gc = circle(m(gc.points));
%         end
%         
        
        function str = char(gc)
            str = sprintf('circle with center %s and radius %s',...
                num2str(gc.center), num2str(gc.radius));
        end
        
        function d = dist(gc, z)
            % Distance between point and circle.
            if ~isinf(gc)
                v = z - gc.center;
                d = abs(abs(v) - gc.radius);
%             else
%                 v = z - gc.points(1);
%                 s = sign(1i*diff(gc.points(1:2)));
%                 d = abs(real(v)*real(s) + imag(v)*imag(s));
            end
        end
        
        function z = intersect(gc1, gc2)
            % Calculate circle intersections.
            
            % Move to "gencircle"?
            
            % Map first circle to the real axis.
            M = mobius(point(gc1, [1/3, 2/3, 1]), [-1, 1, Inf]);
            gc = M(gc2);
            if isinf(gc)
                % Intersect real axis with a line.
                tau = tangent(gc);
                p = gc.points(1);
                if abs(imag(tau)) > 100*eps
                    t = -imag(p)/imag(tau);
                    z = real(p) + t*real(tau);
                else
                    warning(['Circles are close to tangency.\nIntersection', ...
                        ' problem is not well conditioned.'])
                    z = [];
                end
                z = [z, inf];
            else
                % Intersect real axis with a circle.
                rat = -imag(gc.center)/gc.radius;
                if abs(abs(rat) - 1) < 100*eps
                    warning(['Circles are close to tangency.\nIntersection', ...
                        ' problem is not well conditioned.'])
                end
                theta = asin(rat);                    % find one intersection
                theta = theta(isreal(theta));         % may not have one
                theta = unique([theta, pi - theta]);  % may have a second
                z = real(gc.center + gc.radius*exp(1i*theta));
            end
            z = feval(inv(M), z);
        end
        
        function tf = isinf(gc)
            tf = isinf(gc.radius);
        end
        
        function tf = isinside(gc, z)
            if isinf(gc)
                z = (z - gc.center)/tangent(gc, z);  % borked!
                tf = imag(z) > 0;
            else
                tf = abs(z - gc.center) < gc.radius;
            end
        end
             
        function c = uminus(c)
            c = uminus@curve(c);
            c.center = -c.center;
        end
        
        function c = plus(c,z)
            c = plus@curve(c,z);
            c.center = c.center + z;
        end

        function c = mtimes(c,z)
            c = mtimes@curve(c,z);
            c.center = c.center * z;
            c.radius = c.radius * abs(z);
        end
     
   end
    
    
end

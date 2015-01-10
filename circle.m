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
        center = []     
        radius = []
        % In complex terms, a line can be seen as a circle. If the
        % requested circle is a line, accept it and put it here.
        line = []
    end
    
    methods
        function c = circle(varargin)
        %CIRCLE  Circle object.
        %   CIRCLE(CENTER,RADIUS) returns an object representing the
        %   circle with given complex center and real radius. 
        %
        %   CIRCLE(Z) returns the circle through the three complex points 
        %   in the vector Z. Inf is allowed as an entry of Z. If the points
        %   are collinear, a warning is thrown and a ZLINE is returned
        %   instead. 
        %
        %   CIRCLE by itself returns a circle object with unknown
        %   parameters, representing an abstract circle.
        %
        %   See also ZLINE.
        
            if nargin==0
                % Abstract circle with unknown geometry.
                return
            end
            
            switch nargin
                case 1   % vector of three points given
                    z3 = varargin{1};
                    if isa(z3, 'double') && numel(z3) == 3
                        
                        % Find the mobius map from [1,i,-1] to the given
                        % trio. Use a mobius standardmap as an intermediate
                        % step, because the regular mobius constructor
                        % wants to create a circle, causing infinite
                        % recursion.
                        M = mobius(mobius.standardmap(z3)\mobius.standardmap([1, 1i, -1]));
                        
                        % This is the preimage of infinity. Because mobius
                        % maps preserve reflections, its inverse WRT to the
                        % unit circle is the preimage of the center of the
                        % original circle. 
                        zi = pole(M);
                        
                        % If the pole lies on the unit circle, then
                        % infinity is on the original circle--so it's a
                        % line.
                        if abs(abs(zi) - 1) < 10*eps
                            % Infinite radius case.
                            warning('The requested circle is really a line.')
                            z3 = z3(~isinf(z3));
                            center_ = NaN;
                            radius_ = NaN;
                            line_ = zline( z3(1:2) );
                        else
                            % Unit-inverse of zi maps to the center.
                            center_ = M(1/conj(zi));
                            radius_ = abs(z3(1) - center_);
                        end
                    end
                    
                case 2   % center and radius given
                    [center_, radius_] = varargin{:};
                    validateattributes(center_,{'double'},...
                        {'scalar','finite','nonnan'},...
                        'circle','center');
                    validateattributes(radius_,{'double'},...
                        {'scalar','finite','nonnan','positive'},...
                        'circle','radius');
            end
            
            function z = positionfun(t)
                z = center_ + radius_*exp(1i*t);
             end
            
            function zt = tangentfun(t)
                zt = 1i*exp(1i*t);
            end
            
            c = c@closedcurve(@positionfun,@tangentfun,[0 2*pi]);
            c.center = center_;
            c.radius = radius_;
            
            if isnan(center_)
                c.line = line_;
            end
        end
        
%         function c = apply(c, m)
%             if ~isa(m, 'mobius')
%                 error('CMT:NotDefined', ...
%                     'Expected a mobius transformation.')
%             end
%             
%             c = circle(m(c.points));
%         end
%         
        
        function str = char(c)
            if ~isinf(c)
                str = sprintf('circle with center %s and radius %s',...
                    num2str(c.center), num2str(c.radius));
            else
                str = ['(generalized circle) ',char(c.line)];
            end
        end
        
        function d = dist(c,z)
            % Distance between point and circle.
            if ~isinf(c)
                v = z - c.center;
                d = abs(abs(v) - c.radius);
            else
                d = dist(c.line,z);
            end
        end
        
        function z = intersect(c1, c2)
            % Calculate circle intersections.
             
            % Map first circle to the real axis.
            M = mobius(point(c1, [1/3, 2/3, 1]), [-1, 1, Inf]);
            c = M(c2);
            if isinf(c)
                % Intersect real axis with a line.
                tau = tangent(c);
                p = c.points(1);
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
                rat = -imag(c.center)/c.radius;
                if abs(abs(rat) - 1) < 100*eps
                    warning(['Circles are close to tangency.\nIntersection', ...
                        ' problem is not well conditioned.'])
                end
                theta = asin(rat);                    % find one intersection
                theta = theta(isreal(theta));         % may not have one
                theta = unique([theta, pi - theta]);  % may have a second
                z = real(c.center + c.radius*exp(1i*theta));
            end
            z = feval(inv(M), z);
        end
        
        function tf = isinf(c)
            tf = isnan(c.radius);
        end
        
        function tf = isinside(c, z)
            if isinf(c)
                % We'll use the "to the left" definition based on the given
                % tangent direction.
                z0 = point(c.line,0);
                z = (z - z0) / c.line.tangent(0);  
                tf = imag(z) > 0;
            else
                tf = abs(z - c.center) < c.radius;
            end
        end
             
        function c = uminus(c)
            if ~isinf(c)
                c = uminus@curve(c);
                c.center = -c.center;
            else
                c.line = -c.line;
            end
        end
        
        function out = plot(c,varargin)
            if ~isinf(c)
                h = plot@curve(c,varargin{:});
            else
                h = plot(c.line,varargin{:});
            end
            
            if nargout > 0
                out = h;
            end
        end
            
        function c = plus(c,z)
            if ~isinf(c)
                c = plus@curve(c,z);
                c.center = c.center + z;
            else
                c.line = c.line + z;
            end
        end

        function c = mtimes(c,z)
            if ~isinf(c)
                c = mtimes@curve(c,z);
                c.center = c.center * z;
                c.radius = c.radius * abs(z);
            else
                c.line = cline*z;
            end
        end
        
        function out = rsplot(c)
            if isinf(c)
                h = rsplot(c.line);
            else
                h = rsplot@closedcurve(c);
            end
            if nargout > 0
                out = h;
            end
        end
     
        function c = truncate(c)
            if isinf(c)
                c = truncate(c.line);
            end
        end
   end
    
    
end

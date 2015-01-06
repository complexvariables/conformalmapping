classdef ellipse < closedcurve
    % ELLIPSE class represents a parameterized ellipse.
    
    % This file is a part of the CMToolkit.
    % It is licensed under the BSD 3-clause license.
    % (See LICENSE.)
    
    % Copyright Toby Driscoll, 2014.
    % Written by Everett Kropf, 2014.
    
    properties
        yradius = 1
        xradius = 1
        rotation = 0
    end
        
    methods
        function E = ellipse(varargin)
            rot_ = 0;
            if nargin == 1 % eccentricity given
                e = varargin{1};
                if e < 0 || 1 <= e
                    error('CMT:InvalidArgument', ...
                       'Eccentricity must be between 0 and 1.')
                end
                
                m1 = (1-e^2)^(.25);
                m2 = 1/m1;
            elseif nargin >= 2
                if nargin > 3
                    error('CMT:InvalidArgument', 'Too many arguments.')
                end
                [m1,m2] = varargin{1:2};
                if nargin > 2
                    rot_ = varargin{3};
                end
            end
            
            E = E@closedcurve(@point,@tangent,[0 2*pi]);
            
            function z = point(t)
                z = m1*cos(t) + 1i*m2*sin(t);
                if rot_~=0
                    z = z*exp(1i*rot_);
                end
            end
            
            function zt = tangent(t)
                zt = -m1*sin(t) + 1i*m2*cos(t);
                if rot_~=0
                    zt = zt.*exp(1i*rot_);
                end
            end
          
            E.yradius = m1;
            E.xradius = m2;
            E.rotation = rot_;
            
        end
        
        function str = char(e)
            str = sprintf('ellipse with radii %s and %s (eccentricity %.2f)',...
                num2str(e.yradius), num2str(e.xradius), ...
                eccentricity(e) );
        end
        
        function e = eccentricity(E)
            r = radii(E);
            e = sqrt(1-r(1)^2/r(2)^2);
        end
        
        function f = foci(E)
            r = radii(E);
            f = sqrt(r(2)^2-r(1)^2);
            f = exp(1i*E.rotation)*[f,-f];
        end
        
        function r = radii(E)
            % Minor and major radii (sorted in that order).
            r = sort([E.xradius,E.yradius]);
        end
                
        function th = theta_exact(E, t)
            % Boundary correspondence of conformal map to circle.
            % Exact to machine precision.
            % See Henrici, vol. 3, p. 391, equation for theta at bottom of page.
            
            r = radii(E);
            ep = diff(r) / sum(r);
            if ep > 0.95
                warning('The formula is not accurate for this ellipse.')
            end
            
            t = mod(t,2*pi);
            
            % ellipse term magnitude function
            emf = @(m, e) e.^m./(1 + e.^(2*m))./m;
            % term function
            termfun = @(t, m) (-1).^m.*emf(m, ep).*sin(2*m.*t);
            
            % Keep adding terms, 20 at a time, until too tiny to make a change. Good
            % up to about e = 0.95.
            th = t(:);
            for k = 1:60
                m = (k-1)*20 + (1:20);
                th = th + 2*sum(bsxfun(termfun, t, m), 2);
                if all(emf(m(end),ep) < eps(th))
                    break
                end
            end
            th = reshape(th, size(t));
        end
    end
    
end

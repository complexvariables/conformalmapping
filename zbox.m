classdef zbox < polygon
    %ZBOX  Box object.
    %   zbox(BOUNDS) returns an object representing the box (border of a
    %   rectangle) in the complex plane with real part bounded by
    %   BOUNDS(1:2) and imaginary part bounded by BOUNDS(3:4).
    %
    %   zbox(Z1,Z2) returns the box whose opposite corners
    %   are given by the complex points Z1 and Z2.
    %
    %   zbox by itself returns a box object with unknown
    %   parameters, representing an abstract object.
    
    
    % This file is a part of the CMToolbox.
    % It is licensed under the BSD 3-clause license.
    % (See LICENSE.)
    
    % Copyright Toby Driscoll, 2014.
    % (Re)written by Everett Kropf, 2014,
    % adapted from code by Toby Driscoll, originally 20??.
    
    methods
        function c = zbox(varargin)
            
            if nargin==0
                % Abstract with unknown geometry.
                return
            end
            
            switch nargin
                case 1   % vector of bounds given
                    bounds_ = varargin{1};
                case 2
                    zz = cat(1,varargin{:});
                    bounds_ = [min(real(zz)), max(real(zz)),...
                        min(imag(zz)), max(imag(zz)) ];
            end
            
            cr = bounds_([1 2 2 1]);
            ci = bounds_([3 3 4 4]);
            corners = complex(cr,ci);
            
            c = c@polygon(corners);
        end
        
        function b = bounds(r)
            z = r.vertex;
            b = [min(real(z)),max(real(z)),min(imag(z)),max(imag(z))];
        end
        
        function str = char(r)
            b = bounds(r);
            str = sprintf('box with Re(z) in [%g,%g] and Im(z) in [%g,%g]\n',...
                b );
        end
        
        function disp(r)
            disp(char(r))
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
        
    end
    
    
end

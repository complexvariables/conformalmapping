classdef annulus < region
    % ANNULUS is a region bounded by two circles.
    
    % This file is a part of the CMToolbox.
    % It is licensed under the BSD 3-clause license.
    % (See LICENSE.)
    
    % Copyright Toby Driscoll, 2014.
    % (Re)written by Everett Kropf, 2014,
    % adapted from an idea by Toby Driscoll, 20??.
    
    methods
        function A = annulus(varargin)
            %  ANNULUS(CENTER,RAD1,RAD2)
            
            if nargin==0
                icirc = [];  ocirc = [];  % abstract annulus
            else
                
                % TODO: Argument parsing and validity checking.
                [c,r1,r2] = deal(varargin{:});
                if r1 > r2
                    [r2,r1] = deal(r1,r2);  % swap
                end
                icirc = circle(c,r1);
                ocirc = circle(c,r2);
                
            end
            
            A = A@region(icirc,'exterior',ocirc,'interior');
        end
        
        function t = hasgrid(~)
            t = true;
        end
        
        function gd = grid(A,type)
            
            if nargin < 2
                type = 'curves';
            end
            
            rb = [innerRadius(A),outerRadius(A)];
            switch type
                case 'curves'
                    r = linspace(rb(1),rb(2),10);
                    gd = polargrid(A,r(2:end-1),12,'curves',rb);
                    
                case 'mesh'
                    gd = polargrid(A,61,120,'mesh',rb);
                    
                    %             case 'carleson'
                    %                 gd = carlesongrid(A);
                    
                otherwise
                    error('CMT:NotDefined', ...
                        'Grid type "%s" not recognized.', type)
            end
        end
        
        function z = center(A)
            c = outer(A);
            z = c.center;
        end
        
        function r = innerRadius(A)
            c = inner(A);
            r = c.radius;
        end
        
        function r = outerRadius(A)
            c = outer(A);
            r = c.radius;
        end
        
        function m = modulus(A)
            m = innerRadius(A) / outerRadius(A);
        end
        
    end
    
end

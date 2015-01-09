classdef halfplane < region
% HALFPLANE is a region bounded by a line.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% (Re)written by Everett Kropf, 2014,
% adapted from an idea by Toby Driscoll, 20??.

properties
    boundingLine
    interiorPoint
end


methods
    function H = halfplane(bline,interiorpoint)
        %TODO: Argument checking. 
        % We use a polygon as the boundary, as it makes some things easier.
        
        zp = point(bline,[-.5 .5]);
        tau = diff(zp);  % tangent 
        % Is the interior point to the left? If not, reverse orientation.
        s = imag( (interiorpoint-zp(1)) / tau );
        if s < 0
            zp = zp([2 1]);
            tau = -tau;
        end
        p = polygon([zp(1),zp(2),homog(tau,0),homog(-tau,0)]);
        H = H@region(p,'interior');
        H.boundingLine = bline;
        H.interiorPoint = interiorpoint;
    end

    function t = hasgrid(~)
        t = true;
    end
    
    function gd = grid(H,type)
        
        %FIXME: Only for grid aligned with coordinate axes!
        if nargin < 2
            type = 'curves';
        end
        
        switch type
            case 'curves'
                gd = cartesiangrid(H,[-4 4 0 4],12,12,'curves');
                
            case 'mesh'
                gd = cartesiangrid(H,[-4 4 0 4],100,100,'mesh');
               
            otherwise
                error('CMT:NotDefined', ...
                    'Grid type "%s" not recognized.', type)
        end
    end
       
end

end

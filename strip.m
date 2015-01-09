classdef strip < region
%STRIP  Strip region in the complex plane.
    
% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% (Re)written by Everett Kropf, 2014,
% adapted from an idea by Toby Driscoll, 20??.

properties
    ylim
end

methods
    function S = strip(ybounds)
        
        if nargin < 1
            ybounds = [0 1];
        end
        
        p = polygon( [1i*ybounds(1) infcorner(0,0) 1i*ybounds(2) infcorner(pi,pi) ] );
        
        S = S@region(p,'interior');
        S.ylim = ybounds;
    end
       
    function g = grid(S,type)
        
        if nargin < 2
            type = 'curves';
        end
        
        bnd = 0.8./(1-0.8^2);
        bnd = [ -bnd bnd S.ylim ];
        switch type
            case 'curves'
                x = linspace(-0.8,0.8,15);
                x = x./(1-x.^2);
                y = linspace(S.ylim(1),S.ylim(2),11);
                g = cartesiangrid(S,bnd,x,y,'curves');
                
            case 'mesh'
                g = cartesiangrid(S,bnd,100,100,'mesh');
        end
    end
    
    function t = hasgrid(~)
        t = true;
    end

end

end

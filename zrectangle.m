classdef zrectangle < region
%ZRECTANGLE  Rectangular region in the complex plane.
    
% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% (Re)written by Everett Kropf, 2014,
% adapted from an idea by Toby Driscoll, 20??.

methods
    function B = zrectangle(varargin)
        
        if nargin>1 || ~isa(varargin{1},'zbox')
            r = zbox(varargin{:});
        else
            r = varargin{1};
        end
        
        B = B@region(r);
    end
    
    function b = bounds(B)
        b = bounds(B.boundary);
    end
    
    function g = grid(B,type)
        
        if nargin < 2
            type = 'curves';
        end
        
        switch type
            case 'curves'
                g = cartesiangrid(B,bounds(B),11,11,'curves');
                
            case 'mesh'
                g = cartesiangrid(B,bounds(B),100,100,'mesh');
        end
    end
    
    function t = hasgrid(~)
        t = true;
    end

end

end

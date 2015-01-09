classdef zbox < region
% DISK is a region bounded by a circle.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% (Re)written by Everett Kropf, 2014,
% adapted from an idea by Toby Driscoll, 20??.

methods
    function B = zbox(varargin)
        
        if nargin>1 || isa(varargin{1},'zrectangle')
            r = zrectangle(varargin{:});
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
            type = 'mesh';
        end
        
        switch type
            case 'curves'
                g = cartesiangrid(B,bounds(B),11,11,'curves');
                
            case 'mesh'
                g = cartesiangrid(B,bounds(B),100,100,'mesh');
        end
    end

end

end

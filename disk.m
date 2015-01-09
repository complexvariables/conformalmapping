classdef disk < region
% DISK is a region bounded by a circle.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% (Re)written by Everett Kropf, 2014,
% adapted from an idea by Toby Driscoll, 20??.

methods
    function D = disk(center, radius)
        badargs = false;
        switch nargin
            case 0
                C = [];

            case 1
                if isa(center, 'disk')
                    C = center.outerboundary;
                elseif isa(center, 'double') && numel(center) == 3
                    C = circle(center);
                elseif isa(center, 'circle') && ~isinf(center)
                    C = center;
                else
                    badargs = true;
                end

            case 2
                if isa(center, 'double') && isa(radius, 'double') ...
                        && numel(center) == 1 && numel(radius) == 1
                    C = circle(center, radius);
                else
                    badargs = true;
                end

            otherwise
                badargs = true;
        end
        if badargs
            error('CMT:InvalidArgument', ...
                'Expected 3 points or a center and radius.')
        end

        if isempty(C)
            supargs = {};
        else
            supargs = {C};
        end
        
        D = D@region(supargs{:});
        get(D, gridset);
    end

    function gd = grid(D,type)
        
        if nargin < 2
            type = cmtgetpref('grid','type');
        end
        
        switch type
            case 'curves'
                gd = polargrid(D,8,12,'curves');
                
            case 'mesh'
                gd = polargrid(D,60,100,'mesh');

            case 'carleson'
                gd = carlesonGrid(D, opts);
                
            otherwise
                error('CMT:NotDefined', ...
                    'Grid type "%s" not recognized.', type)
        end
    end
    
    function z = center(D)
        c = boundary(D);
        z = c.center;
    end
    
    function r = radius(D)
        c = boundary(D);
        r = c.radius;
    end


    function tf = hasgrid(~)
        tf = true;
    end
    
    function gd = polarGrid(D, opts)
        nrad = opts.numRadialLines;
        ncirc = opts.numCircularLines;
        
        npt = 200;
        c = center(outer(D));
        r = radius(outer(D));

        curves = cell(nrad + ncirc, 1);
        zg = (1:npt)'/(npt+1);
        for k = 1:nrad
            curves{k} = c + r*exp(2i*pi*(k-1)/nrad)*zg;
        end
        zg = exp(2i*pi*(0:npt-1)'/(npt-1));
        for k = 1:ncirc
            curves{nrad + k} = c + r*k/(ncirc+1)*zg;
        end

        gd = gridcurves(curves);
    end
end

end

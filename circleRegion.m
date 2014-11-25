classdef circleRegion < region
%circleRegion represents a region bounded by circles
%
% C = circleRegion(clist)
% C = circleRegion(circle1, circle2, ...)
% Uses circles in cell array clist, or alternatively a list of circles as
% individual arguments, to create a region with zero or more circles as
% boundaries. If the first circle bounds the following circles, the region
% is considered bounded. It's an unbounded region otherwise.
%
% This class represents one of two types of regions:
%   1. A region bounded by a circle which contains one or more circular
%   holes.
%   2. An unbounded region containing one or more circular holes.
% In either case none of the circular boundaries intersect, including
% tangentialy.
% The single circle case is ambiguous, thus the class constructor defaults
% to the unbounded case for a single circle. In this case, and this case
% only, one may use the following to convert to a bounded region:
%
%   C = circleRegion(circle(center, radius));
%   C.bounded = true;
%
% See also circle, region.

% This file is a part of the CMToolit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% E. Kropf, 2014

% UNDOCUMENTED:
% C = circleRegion(..., 'nocheck')
% Skips the bounded/unbounded/intersect check. Use with extreme caution, as
% this breaks the intended class representation.

properties(SetAccess=protected)
    centers
    radii
end

properties(Dependent)
    bounded
end

methods
    function R = circleRegion(varargin)
        args = {};
        if nargin
            if ischar(varargin{nargin})
                clist = varargin(1:nargin-1);
                varargin = varargin(nargin);
            else
                clist = varargin;
                varargin = {};
            end
            if numel(clist) == 1 && isa(clist{1}, 'cell')
                clist = clist{1};
            end
            
            if ~isempty(clist) 
                if any(cellfun(@(x) ~isa(x, 'circle'), clist))
                    error('CMT:InvalidArgument', ...
                        'Expected one or more circles.')
                end
            
                m = numel(clist);
                cv(m,1) = 1i;
                rv(m,1) = 0;
                for j = 1:m
                    cv(j) = clist{j}.center;
                    rv(j) = clist{j}.radius;
                end
                
                if numel(varargin) == 0 || ~strcmp(varargin{1}, 'noCheck')
                    % Check bounded/unbounded/intersect.
                    isinside = circleRegion.circleCheck(cv, rv);
                    if numel(clist) > 1 && all(isinside)
                        args = {clist{1}, clist(2:end)};
                    elseif numel(clist) == 1 || ~any(isinside)
                        args = {clist, 'exteriorto'};
                    else
                        error('CMT:InvalidArgument', ...
                            ['The circles given do not define an expected region.' ...
                            ' See help.'])
                    end
                end
            end
        end
        
        R = R@region(args{:});
        if ~nargin || isempty(clist)
            return
        end
        
        R.centers = cv;
        R.radii = rv;
    end

    function R = apply(R, M)
        % Apply Mobius map to circle region.
        
        if ~isa(M, 'mobius')
            error('CMT:NotDefined', ...
                'Applying class %s to circleRegion is not defined.', ...
                class(M))
        end
        
        C = boundary(R);
        if ~iscell(C)
            C = {C};
        end
        for j = 1:numel(C)
            C{j} = M(C{j});
        end
        R = circleRegion(C);
    end
    
    function R = inv(R)
        % Invert circle region; inner and outer boundaries change roles.
        C = boundary(R);
        if hasouter(R)
            R = region(C(2:end), C{1});
        else
            R = region(C, 'interiorto');
        end
    end
    
    function invfill(R, varargin)
        newplot
        washold = ishold;
        hold on
        
        if isbounded(R)
            fill(region(R.outerboundary{1}, 'exteriorto'))
            fill(region(R.innerboundary, 'interiorto'))
        else
            fill(inv(R))
        end
        
        if ~washold
            hold off
            axis(plotbox(R))
            aspectequal
        end
    end
    
    function tf = isbounded(R)
        tf = R.bounded;
    end
    
    function str = replicate(R)
        s = sprintf('%s = circleRegion({...\n', inputname(1));
        m = numel(R.radii);
        for j = 1:m
            s = sprintf('%s    circle(%s, %s)', s, ...
                num2str(R.centers(j), '%.6g'), num2str(R.radii(j), '%.6g'));
            if j < m
                s = sprintf('%s, ...\n', s);
            else
                s = sprintf('%s});\n', s);
            end
        end
        if nargout
            str = s;
        else
            fprintf('%s', s)
        end
    end
    
    %%%%% get/set %%%%%
    function b = get.bounded(R)
        b = hasouter(R);
    end
    
    function R = set.bounded(R, b)
        if R.m ~= 1
            error('CMT:InvalidOperation', ...
                'Bounded status may only be specified for the single circle case.')
        end
        if ~(islogical(b) && numel(b) == 1)
            error('CMT:InvalidArgument', ...
                'Expected a single logical value.')
        end
        
        bstatus = R.bounded;
        if ~bstatus && b
            % Not bounded, but make it so.
            R.outerboundary = R.innerboundary;
            R.innerboundary = {};
        elseif bstatus && ~b
            % Bounded, but make it not so.
            R.innerboundary = R.outerboundary;
            R.outerboundary = {};
        end
    end
end

methods(Static)
    function isinside = circleCheck(cv, rv)
        m = numel(rv);
        isinside = false(m-1, 1);
        intersect = false;
        for i = 1:m
            for j = i+1:m
                sep = abs(cv(i) - cv(j));
                if sep < rv(i)
                    if rv(i) - sep <= rv(j)
                        intersect = true;
                        break
                    elseif i == 1
                        isinside(j-1) = true;
                    end
                elseif sep > rv(i)
                    if sep <= rv(i) + rv(j)
                        intersect = true;
                        break
                    end
                else
                    intersect = true;
                    break
                end
            end
            if intersect
                error('CMT:InvalidArgument', ...
                    'Circle intersection detected.')
            end
        end
    end
end

end

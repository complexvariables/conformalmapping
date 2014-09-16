classdef circleRegion < region
%circleRegion represents a region bounded by circles
%
% C = circleRegion(clist)
% Uses circles in cell array clist to create region. If clist{1} bounds the
% following circles, the region is considered bounded. It's an unbounded
% region otherwise.
%
% This class represents one of two types of regions:
%   1. A region bounded by a circle which contains one or more circular
%   holes.
%   2. An unbounded region containing one or more circular holes.
% In either case none of the circular boundaries intersect, including
% tangentialy.
% A region interior to a single circle is not currently handled by this
% class. For this use the unitdisk class.
%
% See also circle, region, unitdisk.

% This file is a part of the CMToolit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% E. Kropf, 2014

properties
    centers
    radii
end

methods
    function R = circleRegion(clist)
        args = {};
        if nargin
            cliscell = iscell(clist);
            if (cliscell && ~all(cellfun(@(x) isa(x, 'circle'), clist))) ...
                    || (~cliscell && ~isa(clist, 'circle'))
                error('CMT:InvalidArgument', ...
                    'Expected a single circle or cell array of circles.')
            end
            
            if ~iscell(clist)
                clist = {clist};
            end
            m = numel(clist);
            cv(m,1) = 1i;
            rv(m,1) = 0;
            for j = 1:m
                cv(j) = clist{j}.center;
                rv(j) = clist{j}.radius;
            end
            
            % Check bounded/unbounded/intersect.
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
        
        R = R@region(args{:});
        if ~nargin
            return
        end
        
        R.centers = cv;
        R.radii = rv;
    end
    
    function R = inv(R)
        % Invert circle region; inner and outer boundaries change roles.
        C = boundary(R);
        if hasouter(R)
            R = region(C{2:end}, C{1});
        else
            R = region(C, 'interiorto');
        end
    end
end

end

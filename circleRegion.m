classdef circleRegion < region
%circleRegion represents a region bounded by circles
%
% C = circleRegion(clist)
% C = circleRegion(circle1, circle2, ...)
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

% UNDOCUMENTED:
% C = circleRegion(..., 'nocheck')
% Skips the bounded/unbounded/intersect check. Use with extreme caution, as
% this breaks the intended class representation.

properties(SetAccess=protected)
    centers
    radii
end

methods
    function R = circleRegion(varargin)
        args = {};
        if nargin % && ~isempty(clist)
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
            
            if isempty(clist) || any(cellfun(@(x) ~isa(x, 'circle'), clist))
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
        
        R = R@region(args{:});
        if ~nargin || isempty(clist)
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

classdef cmt
% CMT is a helper function wrapper class.
%
% See methods(cmt) for a list of functions.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

methods(Static)
    function z = bb2z(box)
        % BB2Z axis bounding box to vertices.

        z = complex(box([1, 2, 2, 1]), box([3, 3, 4, 4]));
        z = z(:);
    end

    function box = boundbox(points)
        % BOUNDBOX calculates the bounding box around points.
        %
        % box = boundbox(points) calculates a bounding box around a set of points
        % in the complex plane, and returns coordinates in AXIS format.
        %
        % See also axis, plotbox.

        box([1 3]) = min([real(points(:)) imag(points(:))], [], 1);
        box([2 4]) = max([real(points(:)) imag(points(:))], [], 1);
    end
    
    function clearPrefs
        % Remove appdata prefs storage.
        %
        % Needed since "clear classes" complains if any objects are stored
        % in appdata.
        
        try
            rmappdata(0, 'cmt_prefs')
        catch err
            if strcmp(err.identifier, ...
                    'MATLAB:HandleGraphics:Appdata:InvalidNameAppdata')
                % Nothing there, that's ok.
                return
            else
                rethrow(err)
            end
        end
    end
    
    function dispPrefs
        % Preference dump to console.
        
        prefs = getappdata(0, 'cmt_prefs');
        if isempty(prefs)
            fprintf('No CMT preferences found.\n')
            return
        end
        
        fn = fieldnames(prefs);
        for k = 1:numel(fn)
            fprintf(['\n==========================\n' ...
                       ' For class %s:\n' ...
                       '--------------------------\n'], fn{k})
            disp(prefs.(fn{k}))
        end
    end

    function box = plotbox(points, scale)
        % PLOTBOX returns padded axis coordinates around points.
        %
        % box = plotbox(points) calculates a 1-by-4 array to pass to AXIS which
        % sets a padded square box around the given points.
        %
        % box = plotbox(points, scale) allows setting the padding scale, which
        % defaults to 1.2 times the largest axis of the bounding box.
        %
        % See also axis, boundbox.

        if nargin < 2 || isempty(scale)
            scale = 1.2;
        end

        box = cmt.boundbox(points);

        dbox = scale/2*max(diff(box(1:2)), diff(box(3:4)))*[-1 1];

        box(1:2) = mean(box(1:2)) + dbox;
        box(3:4) = mean(box(3:4)) + dbox;
    end
    
    function result = runTests()
        import matlab.unittest.TestSuite;

        suite_all_tests = TestSuite.fromFolder('test');
        status = run(suite_all_tests);

        if nargout
            result = status;
        end
    end
end

end

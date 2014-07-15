classdef comap
    % This is a dummy class used to provide static methods whose scope is
    % not naturally limited to one of the other classes.
    
    % This file is a part of the CMToolbox.
    % It is licensed under the BSD 3-clause license.
    % (See LICENSE.)
    
    % Copyright Toby Driscoll, 2014.
    % Written by Everett Kropf and Toby Driscoll, 2014.
    
    methods (Static)
        function result = runTests()
            import matlab.unittest.TestSuite;
            
            suite_all_tests = TestSuite.fromFolder('test');
            status = run(suite_all_tests);
            
            if nargout
                result = status;
            end
        end
        
        function box = plotbox(points, scale)
            % Padded axis coordinates around points.
            %
            % box = plotbox(points) calculates a 1-by-4 array to pass to AXIS which
            % sets a padded square box around the given points.
            %
            % box = plotbox(points, scale) allows setting the padding scale, which
            % defaults to 1.2 times the largest axis of the bounding box.
            %
            % See also axis, boundbox.
            
            % Copyright Toby Driscoll, 2014.
            % Written by Everett Kropf, 2014.
            
            if nargin < 2 || isempty(scale)
                scale = 1.2;
            end
            
            box = boundbox(points);
            
            dbox = scale/2*max(diff(box(1:2)), diff(box(3:4)))*[-1 1];
            
            box(1:2) = mean(box(1:2)) + dbox;
            box(3:4) = mean(box(3:4)) + dbox;
        end
        
    end
    
end
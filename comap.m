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
    end
    
end
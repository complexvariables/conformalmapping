function result = runTests()
%runTests runs CMT unit test suite.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

import matlab.unittest.TestSuite;

suite_all_tests = TestSuite.fromFolder('test');
status = run(suite_all_tests);

if nargout
    result = status;
end

end

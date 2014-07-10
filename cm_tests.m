function result = run_tests()
% Run tests in test directory.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

import matlab.unittest.TestSuite;

suite_all_tests = TestSuite.fromFolder('test');
status = run(suite_all_tests);

if nargout
  result = status;
end
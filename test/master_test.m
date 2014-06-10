classdef master_test < matlab.unittest.TestCase
% Common setup/teardown for all tests.
%   * Set path for running tests.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

methods(TestClassSetup)
  function add_the_path(test)
    test.addTeardown(@path, addpath('..'));
  end
end

end
classdef test_mobius < master_test
% Test class for mobius.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

methods(Test)
  function basic_mobius_check(test)
    M = mobius([1, 1i, -1], [0, 5, 7i]);
    R = inv(M);
    cond = isequal(size(M.matrix), [2, 2]) ...
           && isequal(size(R.matrix), [2, 2]);
    test.verifyTrue(cond);
  end
end

end

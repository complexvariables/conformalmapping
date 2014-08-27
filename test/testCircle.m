classdef testCircle < masterTest
% Test class for circle.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

methods(Test)
  function threePointCheck(test)
    gc = circle([0, 5, 7i]);
    test.verifyEqual(radius(gc), hex2num('40113463fa37014d'), ...
        'AbsTol', 1e-15);
  end
  
  function zlineCheck(test)
    z1 = circle([0, 1, inf]);
    cond = isinside(z1, 1i) && ~isinside(z1, -1i);
    test.verifyTrue(cond);
  end
end

end

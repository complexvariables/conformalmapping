classdef test_homog < master_test
% Unit tests for HOMOG class.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

methods(Test)
  function create_homog_numbers(test)
    h = [0; homog(1,0); homog(1i,0); 1i];
    test.verifyEqual(double(h), [0; Inf; Inf; 1i]);
  end
  
  function check_homog_angle(test)
    h = [0; homog(1,0); homog(1i,0); 1i];
    test.verifyEqual(angle(h), pi/2*[0; 0; 1; 1]);
  end
end

end
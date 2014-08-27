classdef testRegion < masterTest

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

methods(Test)
    function basicRegionCheck(test)
        c0 = circle(0, 1);
        c1 = circle(0.41+0.16i, 0.24);
        c2 = circle(-0.53+0.09i, 0.12);
        c3 = circle(-0.10-0.34i, 0.25);
        
        r = region(c0, {c1, c2, c3});
        
        cond = hasouter(r) & r.numouter == 1 & hasinner(r) & r.numinner == 3;
        test.verifyTrue(cond);
    end
end

end

classdef test_mobius < master_test
% Test class for mobius.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties
    maps
    data
end

methods(TestMethodSetup)
    function createMaps(test)
        z3 = [1, 1i, -1];
        w3 = [-1, -2i, 0];
        m1 = mobius(z3, w3);
        s3 = [0, 3, 1i];
        m2 = mobius(w3, s3);
        test.maps = struct('m1', m1, 'm2', m2);
        test.data = struct('z3', z3, 'w3', w3, 's3', s3);
    end
end

methods(Test)
    function basicMobiusCheck(test)
        M = mobius([1, 1i, -1], [0, 5, 7i]);
        R = inv(M);
        cond = isequal(size(M.matrix), [2, 2]) ...
            && isequal(size(R.matrix), [2, 2]);
        test.verifyTrue(cond);
    end
  
    function compositionCheck1(test)
        z3 = test.data.z3;
        m1 = test.maps.m1;
        m2 = test.maps.m2;
        m = m2*m1;
        test.verifyLessThan(norm(m2(m1(z3)).' - m(z3).'), 1e-15)
    end
    
    function compositionCheck2(test)
        z3 = test.data.z3;
        m1 = test.maps.m1;
        m2 = test.maps.m2;
        m = conformalmap(m1, m2);
        test.verifyLessThan(norm(m2(m1(z3)).' - m(z3).'), 1e-15)
    end
    
    function plotCheck(test)
        m1 = test.maps.m1;
        m2 = test.maps.m2;
        h = figure;
        plot(m1)
        clf
        plot(m2*m1)
        close(h)
    end
end

end

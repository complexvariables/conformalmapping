classdef testPolygon < masterTest
% Test class for polygon class

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

methods(Test)
    function basicPolygonCheck(test)
        p1 = polygon([0, 1, 1+1i, 1i]);
        % Took out because polygon temporarily changed to match SCT.
        %     p2 = polygon([0, homog(1, 0), homog(1i, 0)]);
        cond = sum(angle(p1)) == 2; % && sum(angle(p2)) == 0;
        test.verifyTrue(cond);
    end
end

end

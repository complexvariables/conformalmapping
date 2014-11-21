%% New diskmap test.
clear


%%
% 'L' polygon.

p = polygon([1i; -1+1i; -1-1i; 1-1i; 1; 0]);
f = diskmap(p)

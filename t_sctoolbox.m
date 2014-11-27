%% New diskmap test.
clear


%%
% 'L' polygon.

p = polygon([1i; -1+1i; -1-1i; 1-1i; 1; 0]);


%%

f = hplmap(p);
g = diskmap(f);


%%

g = center(g, 0.4-0.6i);
plot(exp(1i*pi/4)*(g-center(g)))

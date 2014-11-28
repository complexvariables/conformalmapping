%% New diskmap test.
clear


%%
% 'L' polygon.

p = polygon([1i; -1+1i; -1-1i; 1-1i; 1; 0]);


%%

% f = hplmap(p);
% plot(f)


%%

% f = diskmap(f);
% f = center(f, 0.4-0.6i);
% plot(exp(1i*pi/4)*(f-center(f)))


%%

% f = extermap(p);


%%

f = rectmap(p, [1 2 4 5]);

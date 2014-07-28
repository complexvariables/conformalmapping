%% Try to get a Carleson grid.
clear


%%

gd = carleson(unitdisk);


%%
G = sample_splines(3);
opts = szset('numCollPts', 1024, 'numFourierPts', 1024);
f = szmap(G, 0, opts);
g = f';


%%
clf
hold on
plot(f(gd))
plot(g(carleson(diskex(0, 1))))
plot(G)
axis(plotbox(G))
aspectequal
axis off

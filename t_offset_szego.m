%% Test offset exterior Szego map.
clear


%%
G1 = splinep([ ...
    1.9532 + 1.2164i; 0.9123 + 2.2222i; 0.0234 + 1.2398i
    1.1462 + 0.0819i ...
]);
a1 = 0.998057406068837 + 1.15689434530228i;

% G1 = G1 - a1;
% a1 = 0;

%%
% Map to region containing origin by complex inversion. Get Szego kernel
% for this interior map, create exterior map by double inversion.

Gi = cinvcurve(G1, a1);

S = szego(Gi, 'confCenter', 0);

nF = get(szego, 'numFourierPts');
t = (0:nF-1)'/nF*2*pi;
s = invtheta(S, t);


%%
% % Inverse theta function fails. Check that tangents are being done correctly.
% 
% ntt = 50;
% tt = (0:ntt-1)'/ntt;
% 
% figure(1), clf
% subplot(1,2,1)
% z = G1(tt);
% zt = tangent(G1, tt);
% quiver(real(z), imag(z), real(zt), imag(zt));
% hold on
% plot(G1)
% aspectequal
% axis(plotbox(G1))
% 
% subplot(1,2,2)
% z = Gi(tt);
% zt = tangent(Gi, tt);
% quiver(real(z), imag(z), real(zt), imag(zt));
% hold on
% plot(Gi)
% aspectequal
% axis(plotbox(Gi))

%% Specification by example -- Szego kernel map.
% This file best viewed via the command
%
%   open(publish('spec_mobius_grids')).
%
clear

%%
% A spline to work with. Although technically the Szego kernel represents
% an operator on L^2, the code only works well for now on a closed curve
% that thinks it is C^2.
G = splinep([ ...
    0.2398 + 0.6023i; 0.3567 + 1.0819i; 0.2632 + 1.5965i
    -0.5205 + 1.7485i; -1.0585 + 1.1170i; -1.0702 + 0.5088i
    -0.5906 + 0.0994i; -0.7778 - 0.4269i; -1.2924 - 0.6140i
    -1.4561 - 1.2456i; -0.5439 - 1.3509i; 0.2515 - 1.0702i
    0.3099 - 0.6023i; 0.7427 - 0.5906i; 1.1053 - 0.1813i
    1.2807 + 0.3567i ...
]);

%%
% Make the map and draw it.
f = szmap(G, 0);
plot(f)

%%
% So how is this done? Start by computing the Szego kernel for |G|.
S = szego(G, 0);
disp(S)


%%
% Take some (20 here) evenly spaced points around G, and consider the boundary
% correspondence.
t = (0:19)'/20;

clf
subplot(1,2,1)
plot(G), hold on, plot(G(t), 'rd')
axis off
subplot(1,2,2)
plot(circle(0, 1)), hold on, plot(exp(1i*theta(S, t)), 'rd')
axis off

%%
% Now consider the inverse boundary correspondence.

clf
subplot(1,2,1)
plot(circle(0, 1)), hold on, plot(exp(2i*pi*t), 'rd')
axis off
subplot(1,2,2)
plot(G), hold on, plot(G(invtheta(S, 2*pi*t)), 'rd')
axis off

%%
% Using this inverse boundary correspondence and an FFT, create a conformal
% map.
N = 512;
th = 2*pi*(0:N-1)'/N;
t = invtheta(S, th);
w = G(t);
c = fft(w)/N;
f = @(z) polyval(flipud(c), z);

%%
% Aaaaaand draw.
gd = grid(unitdisk);
gdi = cell(1, numel(gd));
for k = 1:numel(gd)
  gdi{k} = f(gd(k));
end
gdi = gridcurves(gdi);

subplot(1,2,1)
plot(gd)
subplot(1,2,2)
plot(gdi)


%%
% Now we'll do a shape fingerprint. Let's use a new shape.
clear
G = splinep([ ...
    0.8129 + 0.0643i; 1.1871 + 0.5673i; 1.0468 + 1.3275i
    0.2865 + 1.6316i; -0.0643 + 1.2339i; -0.4035 + 0.6257i
    -1.0117 + 0.8480i; -1.4211 + 0.5789i; -1.4795 - 0.2982i
    -1.1871 - 1.1520i; -0.2865 - 0.5205i; 0.4503 - 0.7076i ...
]);


%%
% Create a map for the interior and exterior of the curve. Note the use of
% the conjugate inverse for the exterior map.
fin = szmap(G, 0);
g = szmap(G', 0);
fout = @(z) 1./conj(g(1./conj(z)));

%%
% Compute the outer grid.
gde = grid(diskex(0, 1));
gdi = cell(1, numel(gde));
for k = 1:numel(gde)
  gdi{k} = fout(gde{k});
end
gdi = gridcurves(gdi);


%%
% Compute the shape fingerprint, which is just the graph of the boundary
% correspondence angles between the shape and the unit disk with respect to the
% interior and exterior maps.
Sin = fin.kernel_;
Sout = g.kernel_;
s = (0:199)'/200;                 % Points in [0, 1) for the 
                                  % parameterized  boundary G.
tin = unwrap(theta(Sin, s));      % Interior boundary image angles.
tout = unwrap(theta(Sout, s));    % Exterior boundary image angles.


%%
clf
subplot(1,2,1), hold on
plot(gdi)
plot(fin)
plot(G(0), 'r.')
axis(plotbox(G))
aspectequal
axis off

subplot(1,2,2)
plot(tin/pi, tout/pi)
xlabel('\theta_{inner}/\pi')
ylabel('\theta_{outer}/\pi')
axis([0 2 0 2])
grid on
aspectequal
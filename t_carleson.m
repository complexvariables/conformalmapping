%% Try to get a Carleson grid.
clear


%%
% Zipper Carleson example.

gdz = zipper.cgrid;


%%
% Try the |mu| levels of circles. First level has |nu| radial lines. Second
% level adds |nu| radial lines for a total of |2*nu|. Third level adds |2*nu|
% radial lines for a total of |4*nu|. At level |mu|, there should be a
% total of |2^(mu-1)*nu| radial lines.

mu = 6;
nu = 32;

gc = cell(mu + 2^(mu-1)*nu, 1);
r = zeros(mu, 1);
r(1) = 0.6;

% Level 1 circles
ncp = 200;
gc{1} = r(1)*exp(2i*pi*(0:ncp-1)'/(ncp-1));

% Level 1 radial lines.
ppul = 200; % Number of radial line points per unit length.
ncr = ceil(ppul*(1 - r(1)));
rln = linspace(r(1), 1 - 1e-8, ncr);
for k = 1:nu
  gc{1 + k} = rln*exp(2i*pi*(k-1)/nu);
end

idx = 2 + nu;
for j = 2:mu
  ncp = j*ncp;
  dr = 1 - r(j-1);
  r(j) = 0.5*dr + r(j-1);
  gc{idx} = r(j)*exp(2i*pi*(0:ncp-1)'/(ncp-1));
  
  ncr = ceil(j*ppul*(1 - r(j)));
  rln = linspace(r(j), 1 - 1e-6, ncr);
  nuj = 2^(j-2)*nu; % number of radial lines being added
  for k = 1:nuj
    gc{idx + k} = rln*exp(1i*pi/nuj*(2*k-1));
  end
  
  idx = idx + 1 + nuj;
end

gd = gridcurves(gc);


%%
G = sample_splines(3);
opts = szset('numCollPts', 1024, 'numFourierPts', 1024);
f = szmap(G, 0, opts);


%%
clf
hold on
plot(f(gd))
plot(G)
aspectequal
axis off

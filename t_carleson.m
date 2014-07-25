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

mu = 5;
nu = 32;

gc = cell(1 + mu + 2^(mu-1)*nu, 1);
r = 0.6;

% Level 0 circle.
ncp = 200;
gc{1} = r*exp(2i*pi*(0:ncp-1)'/(ncp-1));

% Base number of radial line points per unit length.
ppul = 200;

idx = 1;
for j = 1:mu
  if j > 1
    nuj = 2^(j-2)*nu;
  else
    nuj = nu;
  end
  ncr = ceil(j*ppul*(1 - r));
  rln = linspace(r, 1 - 1e-8, ncr)';
  dt = 2*pi/nuj;
  off = (j > 1)*dt/2;
  for k = 1:nuj
    gc{idx + k} = rln*exp(1i*(off + (k-1)*dt));
  end
  
  idx = idx + nuj + 1;
  r = (1 + r)/2;
  np = (j+1)*ncp;
  gc{idx} = r*exp(2i*pi*(0:np-1)'/(np-1));
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

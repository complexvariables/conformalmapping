function handle = plot(p,varargin)
%PLOT Plot a polygon.

z = p.vertex;
alpha = p.angle;
zh = p.hvertex;
zhind = p.hindex;
incoming = p.incoming;

if isempty(z), return, end

n = length(z);
atinf = isinf(z);

newplot
washold = ishold;
autoscale = strcmp(get(gca,'xlimmode'),'auto') & ...
  strcmp(get(gca,'ylimmode'),'auto') | ~ishold;

zf = z(~atinf);
if autoscale
  lim = [min(real(zf)),max(real(zf)),min(imag(zf)),max(imag(zf))];
  maxdiff = max(diff(lim(1:2)),diff(lim(3:4)));
  if maxdiff < 100*eps, maxdiff = 1; end
  fac = .6 + .15*(any(atinf));
  lim(1:2) = mean(lim(1:2)) + fac*maxdiff*[-1,1];
  lim(3:4) = mean(lim(3:4)) + fac*maxdiff*[-1,1];
else
  lim = axis;
end
R = 10*max([1,lim(2)-lim(1),lim(4)-lim(3)]);

% An unbounded polygon will be truncated to make it safe for plotting
zplot = vertex( truncate(p) );

% % Accumulate the vertices for ploting.
% % An unbounded polygon will be truncated to make it safe for plotting and
% % filling.
% j = 0;
% zplot = [];
% for k = 1:n
%   j = j+1;
%   if ~isinf(z(k))
%     zplot(j) = z(k);
%   else
%     kp = mod(k-2,n) + 1;
%     za = z(kp) + R*sign(incoming(k));
%     kn = mod(k,n) + 1;
%     zb = z(kn) - R*sign(incoming(kn));
%     % Insert a NaN to prevent drawing between these two "truncated
%     % infinities"
%     zplot(j:j+2) = [za NaN zb];
%     j = j+2;
%   end
% end

zplot = zplot([1:end 1]);
h = plot(real(zplot),imag(zplot),'linewidth',1.25);
set(h,varargin{:});
if autoscale, axis(lim), axis square, end

if nargout > 0
  handle = h;
end


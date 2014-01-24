function box = clippingbox(data,scale)

x = real( data(~isinf(data)) );
y = imag( data(~isinf(data)) );
x = [ min(x) max(x) ];
y = [ min(y) max(y) ];

center = [ sum(x) sum(y) ]/2;
size = [ diff(x) diff(y) ]/2;
box(1:2) = center(1) + size(1)*scale*[-1 1];
box(3:4) = center(2) + size(2)*scale*[-1 1];

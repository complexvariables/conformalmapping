%% Specification by example for Mobius.
clear

z3 = [1, 1i, -1];
w3 = [-1, -2i, 0];
m1 = mobius(z3, w3);
figure,plot(m1)

s3 = [0, 3, 1i];
m2 = mobius(w3, s3);
figure,plot(m2)

% composition!
m = m2*m1;
fprintf('composition comparison: %.g\n', norm(m2(m1(z3)).' - m(z3).'))

mf = conformalmap(m1, m2);
fprintf('other comparison: %.g\n', norm(m2(m1(z3)).' - mf(z3).'))

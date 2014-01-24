function display(zeta)

n = size(zeta.numer);
sizestr = [ sprintf('%i-by-',n(1:end-1)), sprintf('%i',n(end)) ];
fprintf(['\n\t' sizestr ' array of homogeneous coordinates:\n\n'])
fprintf('numerator = \n\n')
disp(zeta.numer)
fprintf('\ndenominator = \n\n')
disp(zeta.denom)
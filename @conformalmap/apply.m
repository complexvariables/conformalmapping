function gnew = apply(f,g)

if isa(g,'gridcurves')
  c = curve(g);
  for m = 1:length(c)
    %paramfun = @(t) feval( f, c{m}(t) );
    [T,Z{m}] = adaptplot( @paramfun, [100*eps 1-100*eps] );
  end
  gnew = gridcurves(f.image,Z{:});
else
  error('Use APPLY to apply a map to a GRIDCURVES object')
end

  function z = paramfun(t)
    % The weird syntax seems to be needed because of a matlab overloading bug.
    S = struct('type','()','subs',{ {c{m}(t)} });
    z = subsref(f,S);
  end

end
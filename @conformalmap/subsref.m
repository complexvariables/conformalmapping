function w = subsref(f,S)

if length(S)==1 & strcmp(S.type,'()') & length(S.subs)==1
  if iscell(f.function)
    % Composite map; call recursively
    for m = 1:length(f.function)
      S.subs{1} = subsref(f.function{m},S);
    end
    w = S.subs{1};
  else
    z = S.subs{1};
    if isa(z,'double')
      w = f.function(z);
    end
  end
else
  error('Unsupported subscript reference')
end

    
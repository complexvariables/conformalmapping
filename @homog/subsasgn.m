function a = subsasgn(a,s,b)


if strcmp(s.type,'()')
  if length(s.subs)==1
    a = homog(a);
    b = homog(b);
    index = s.subs{1};
    a.numer(index) = b.numer;
    a.denom(index) = b.denom;
  else
    error('HOMOG objects support linear indexing only')
  end
else
  error('Unsupported assignment syntax')
end
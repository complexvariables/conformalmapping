function z = subsref(this,S)

if length(S) == 1 & strcmp(S.type,'()')
  z = point(this,S.subs{1});
else
  error('Only a single parenthesized subscript is allowed')
end

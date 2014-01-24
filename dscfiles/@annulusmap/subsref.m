function w = subsref(map,S)

if strcmp(S.type,'()') 
  w = feval(map,S.subs{1});
end
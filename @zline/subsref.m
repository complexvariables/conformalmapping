function z = subsref(this,S)
%SUBSREF 

%   Copyright (c) 1998, 2006 by Toby Driscoll.

if length(S) == 1 & strcmp(S.type,'()')
  z = point(this,S.subs{1});
else
  error('Unsupported referencing syntax')
end

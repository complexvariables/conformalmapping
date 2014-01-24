h36053
s 00002/00000/00046
d R 1.2 97/04/24 09:55:57 tad 2 1
c Added isdiff and isinv fields.
c 
e
s 00046/00000/00000
d D 1.1 97/04/24 09:52:54 tad 1 0
c 
e
u
U
f b 
f e 0
t
T
I 1
function M = scmap(poly,opt)
%SCMAP Construct generic Schwarz-Christoffel map object.
%   SCMAP(P) creates a generic parent scmap object whose target polygon
%   is given by P. SCMAP(P,OPTIONS) accepts an options structure
%   produced by SCMAPOPT.
%   
%   SCMAP(M), where M is already an scmap object, returns M. SCMAP by
%   itself returns an object with empty polygon and default options.
%   
%   You do not need to create an scmap object directly. Use one of the
%   specific child classes instead.
%   
%   See also SCMAPOPT, and classes DISKMAP, HPLMAP, EXTERMAP, STRIPMAP,
%   RECTMAP, CRMAP, CRRECTMAP.

%   Copyright 1997 by Toby Driscoll. Last updated %G%.

superiorto('double');

if nargin == 0
  % Leave everything empty
  poly = [];
  opt = [];
else

  % Branch based on class of first argument
  switch class(poly)
  case 'scmap'
    % Self-return
    M = poly;
    return
  case 'polygon'
    if nargin == 1
      opt = [];
    end
  otherwise
    msg = 'Expected ''%s'' to be of class polygon or scmap.';
    error(sprintf(msg,inputname(1)))
  end

end
  
M.polygon = poly;
M.options = scmapopt(opt);

M = class(M,'scmap');
E 1

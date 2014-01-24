function M = scmap(domain,image,opt)
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
%   RECTMAP, CRDISKMAP, CRRECTMAP.

%   Copyright 1998 by Toby Driscoll.
%   $Id: scmap.m,v 2.1 1998/05/10 04:24:18 tad Exp $

superiorto('double');

M.domain = [];
M.image = [];
M.opt = [];

switch nargin
  case 0
    M = class(M,'scmap',conformalmap);
    return
  case 1 
    if isa(domain,'scmap')
      M = domain;
      return
    end
  case 3
    M.domain = domain;
    M.image = image;
    M.opt = opt;
end

M = class(M,'scmap',conformalmap([],domain,image));

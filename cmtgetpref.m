function value = cmtgetpref(group,prop)
%CMTGETPREF  Get preference(s) for the CM Toolkit.
%  CMTGETPREF(GROUP,PARAM) returns the value of the Conformal Mapping
%  Toolkit preference parameter PARAM in the group GROUP.
%
%  CMTGETPREF(GROUP) returns a structure with all the parameter values in
%  the group.

%  Copyright 2015 by Toby Driscoll.

cmtgroup = ['CMT',group];

if ~isappdata(0,cmtgroup)   % first call; create prefs group
    cmtsetpref(group,'factory');
end

data = getappdata(0,cmtgroup);

if nargin < 2    % return all values in group
    value = data;
   
else    % return a particular value
    if ~isfield(data,prop)
        error('Property %s of group %s is not recognized.',prop,group)
    end
    
    value = data.(prop);
    
end

end
function value = cmtgetpref(group,prop)

cmtgroup = ['CMT',group];
if ~isappdata(0,cmtgroup)
    error('group %s is not recognized.',group)
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
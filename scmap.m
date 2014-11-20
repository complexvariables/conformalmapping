classdef scmap < conformalmap
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

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.

properties(SetAccess=protected)
    polygon         % Polygon for the map.
    options         % Option structure for the map.
end

methods
    function map = scmap(poly, opt)
        if ~nargin
            poly = [];
            opt = [];
        else
            % Branch based on class of first argument.
            switch class(poly)
                case 'polygon'
                    if nargin == 1
                        opt = [];
                    end
                otherwise
                    msg = 'Expected ''%s'' to be of class polygon or scmap.';
                    error(msg,inputname(1))
            end
        end
       
        map.polygon = poly;
        map.options = sctool.scmapopt(opt);
    end
    
    function s = char(~)
        s = '  scmap object (generic)';
    end
    
    function Md = diff(M)
        %DIFF   Differentiated SC map object.
        %   DIFF(M) returns an object formally representing the derivative of
        %   the map M.
        
        Md = scmapdiff(M);
    end
    
    function M = diskmap(M)
        %DISKMAP Convert generic Schwarz-Christoffel map object to disk map.
        %   DISKMAP(M) creates a diskmap object based on the polygon and
        %   options contained in M.
        %
        %   See the DISKMAP class documentation.

        M = diskmap(M.polygon,M.options);
    end
    
    function disp(f)
        s = char(f);
        if isstr(s)
            disp(s)
        elseif iscell(s)
            fprintf('\n  SC %s:\n\n',class(f));
            for n = 1:length(s)
                disp(s{n})
            end
        end
    end
    
    function display(f)
        fprintf('\n%s =\n',inputname(1))
        disp(f)
    end
    
    function M = extermap(M)
        %EXTERMAP Convert generic Schwarz-Christoffel map object to exterior map.
        %   EXTERMAP(M) creates a extermap object based on the polygon and
        %   options contained in M.
        %
        %   See the EXTERMAP class documentation.
        
        M = extermap(M.polygon,M.options);
    end
    
    function varargout = get(map,varargin)
        %GET    Get map parameters.
        %   Each SC map stores data needed to compute with the map. GET is a
        %   generic access to these data.
        %
        %   [VAL1,VAL2,...] = GET(F,'PROP1','PROP2',...) returns the values of
        %   the map F corresponding to the requested properties. Valid properties
        %   vary by map type. Use GET(F) to see a list for the type associated
        %   with F. Note that field name abbreviations are NOT allowed.
        
        if nargin==1   % typeout available fields
            
            % Get names.
            names = fieldnames(map);
            % Strip out the generic scmap (this class!).
            names( strmatch('scmap',names) ) = [];
            % Add in the always-present 'polygon' and 'options'.
            names = { 'polygon', 'options', names{:} };
            
            fprintf('\nAvailable fields in class %s:\n',class(map))
            fprintf('  %s\n',names{:})
            fprintf('\n')
            
        else   % get values
            
            for j = 1:length(varargin)
                switch lower(varargin{j}(1:min(3,length(varargin{j}))))
                    case 'pol'
                        varargout{j} = polygon(map);
                    case 'opt'
                        varargout{j} = options(map);
                    otherwise
                        try
                            % This is an affront to OO programming! However, it's very
                            % convenient.
                            m = struct(map);
                            varargout{j} = m.(lower(varargin{j}));
                        catch
                            warning(sprintf('Field ''%s'' not recognized.\n',varargin{j}))
                            varargout{j} = [];
                        end
                end
            end
            
        end
    end
    
    function M = hplmap(M)
        %HPLMAP Convert generic Schwarz-Christoffel map object to half-plane map.
        %   HPLMAP(M) creates a hplmap object based on the polygon and
        %   options contained in M.
        %
        %   See the HPLMAP class documentation.
        
        M = hplmap(M.polygon,M.options);
    end
    
    function Mi = inv(M)
        %INV    Inverse SC map object.
        %   INV(M) returns an object formally representing the inverse of the
        %   map M.
        
        Mi = scmapinv(M);
    end
    
    function M = minus(M,a)
        %   Subtract a contant from the image of an SC map.
        
        M = M + (-a);
    end
    
    function M = mrdivide(M,a)
        %   Divide the image of an SC map by a constant.
        
        if isa(a,'double')
            M = M * (1/a);
        else
            error('Cannot divide by an SC map.')
        end
    end
    
    function M = mtimes(M,c)
        %   Scale the image of a map by a complex constant.
        
        % Usually this will be invoked by a child object, which needs to adjust
        % its own constant as well.
        
        % May need to swap arguments
        if isa(M,'double') & isa(c,'scmap')
            tmp = M;
            M = c;
            c = tmp;
        end
        
        M.polygon = c*M.polygon;
    end
    
    function M = plus(M,a)
        %   Add a constant to the image of an SC map (i.e., translate image).
        
        % May need to swap arguments
        if isa(M,'double') & isa(a,'scmap')
            tmp = M;
            M = a;
            a = tmp;
        end
        
        if length(a)==1 & isa(a,'double')
            M.polygon = M.polygon + a;
        else
            error('Addition is not defined for these operands.')
        end
    end
    
    function z = prevertex(M)
        %PREVERTEX Extract a vector of the prevertices of an S-C map.
        
        tmp = parameters(M);
        if strmatch('prevertex',fieldnames(tmp))
            z = tmp.prevertex;
        else
            msg = sprintf('Prevertices not defined for map of class %s\n',class(M));
            error(msg)
        end
    end
    
    function opt = scmapopt(M,varargin)
        %SCMAPOPT Options structure for a Schwarz--Christoffel map object.
        %   Same as the regular SCMAPOPT, but the first argument is the options
        %   structure contained in the map M.
        
        opt = sctool.scmapopt(M.options,varargin{:});
    end
    
    function M = stripmap(M)
        %STRIPMAP Convert generic Schwarz-Christoffel map object to strip map.
        %   STRIPMAP(M) creates a stripmap object based on the polygon and
        %   options contained in M.
        %
        %   See the STRIPMAP class documentation.
        
        M = stripmap(M.polygon,M.options);
    end
    
    function wp = subsref(M,S)
        %SUBSREF Evaluate map by subscript notation.
        %   M(ZP), where M is a SC map and ZP is a vector of points in the
        %   canonical domain of the map, returns the image of the points in ZP.
        %
        %   This just a synonym for EVAL(M,ZP).
        %
        %   See also EVAL, SCMAP.
        
        if length(S) == 1 & strcmp(S.type,'()')
            wp = eval(M,S.subs{1});
        else
            error('Only syntax for SCMAP is a single parenthesized subscript.')
        end
    end
    
    function M = uminus(M)
        %   Negate the image of an SC map.
        
        M = (-1)*M;
    end
    
    function M = uplus(M)
    end
end

end

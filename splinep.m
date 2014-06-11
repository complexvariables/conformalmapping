classdef splinep < closedcurve
% SPLINEP class represents a periodic spline in the plane.

% This file is a part of the CMToolkit.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014.

properties
  xk_                               % knot x-coordinates
  yk_                               % knot y-coordinates
  
  pp_                               % piecewise polynomial array
  chordal_arclength_ = 0            % total chordal arclength
end

methods
  function S = splinep(varargin)
    if ~nargin
      return
    end
    
    needpts = false;
    tofignum = [];
    
    switch nargin
      case 1
        tmp = varargin{1};
        if isa(varargin{1}, 'splinep')
          xk = tmp.xk_;
          yk = tmp.yk_;
        elseif isempty(tmp) || (numel(tmp) == 1 && ...
            ishghandle(tmp(1)) && strcmp(get(tmp, 'type'), 'figure'))
          % Need to prompt for points.
          needpts = true;
          tofignum = tmp;
        else
          % Assume it's a complex vector.
          xk = real(tmp(:));
          yk = imag(tmp(:));
        end
        
      case 2
        % Just assume it's two real vectors.
        xk = varargin{1}(:);
        yk = varargin{2}(:);
        if numel(xk) ~= numel(yk)
          error('CMT:InvalidArgument', ...
                'Input vectors must have the same number of elements.')
        end
        
      otherwise
        error('CMT:InvalidArgument', ...
              'Expected a complex vector or two real vectors.')
    end
    
    if needpts
      [xk, yk] = splinep.get_pts_(tofignum);
    end
    
    if xk(1) ~= xk(end)
      xk = [xk; xk(1)];
      yk = [yk; yk(1)];
    end
    
    % Superclass constructor here.
    
    % Spline data.
    S.xk_ = xk;
    S.yk_ = yk;
    [S.pp_, S.chordal_arclength_] = splinep.spline_(xk, yk);
  end
  
  function out = apply(S, op)
    % Apply operator to spline knots.
    switch class(op)
      case 'mobius'
        M = matrix(op);
        z = zpts(S);
        out = splinep((M(1)*z + M(3))./(M(2)*z + M(4)));
        
      otherwise
        error('CMT:NotDefined', 'Application of %s to %s is not defined.', ...
              class(op), class(S))
    end
  end
  
  function L = arclength(S)
    L = S.chordal_arclength_;
  end
  
  function disp(S)
    fprintf('splinep object:\n\n')
    fprintf('  defined with %d spline knots,\n', numel(S.xk_))
    lstr = strtrim(evalc('disp(arclength(S))'));
    fprintf('  total chordal arc length %s\n\n', lstr)
  end
  
  function S = minus(S, z)
    S = plus(S, -z);
  end
  
  function S = mrdivide(S, z)
    if isa(z, 'splinep')
      [z, S] = deal(S, z);
    end
    if ~(isa(z, 'double') && numel(z) == 1)
      error('CMT:NotDefined', ...
            'Only scalar division allowed.')
    end
    S = mtimes(S, 1/z);
  end
  
  function S = mtimes(S, z)
    % Scalar multiplication.
    if ~isa(S, 'splinep')
      [z, S] = deal(S, z);
    end
    if isa(z, 'double') && numel(z) == 1
      S = splinep(z*zpts(S));
    else
      error('CMT:NotDefined', ...
            'Only scalar multiplication defined.')
    end
  end
  
  function S = plus(S, z)
    % Translate a spline.
    if isa(z, 'splinep')
      [z, S] = deal(S, z);
    end
    if ~(isa(z, 'double') && numel(z) == 1)
      error('CMT:NotDefined', ...
            'Only translation by a scalar allowed.')
    end
    S = splinep(zpts(S) + z);
  end
  
  function z = point(S, t)
    t = modparam(S, t)*S.chordal_arclength_;
    z = complex(ppval(S.pp_{1,1}, t), ...
                ppval(S.pp_{2,1}, t));
  end
  
  function replicate(S)
    % Print in format for pasting into scripts, etc.
    fprintf('%s = splinep([ ...\n    ', inputname(1));
    z = zpts(S);
    n = numel(z) - 1; % don't repeat first point
    for k = 1:n
      fprintf('%.4f', real(z(k)));
      if imag(z(k)) < 0
        fprintf(' - ');
      else
        fprintf(' + ');
      end
      fprintf('%.4fi', abs(imag(z(k))));
      
      if k ~= n
        if mod(k, 3)
          fprintf('; ');
        else
          fprintf('\n    ');
        end
      end
    end
    
    fprintf(' ...\n]);\n');
  end
  
  function z2 = second(S, t)
    t = modparam(S, t)*S.chordal_arclength_;
    z2 = complex(ppval(S.pp_{1,3}, t), ...
                 ppval(S.pp_{2,3}, t));
  end
  
  function zt = tangent(S, t)
    t = modparam(S, t)*S.chordal_arclength_;
    zt = complex(ppval(S.pp_{1,2}, t), ...
                 ppval(S.pp_{2,2}, t));
  end
  
  function S = uminus(S)
    S = splinep(-zpts(S));
  end
  
  function x = xpts(S)
    x = S.xk_;
  end
  
  function y = ypts(S)
    y = S.yk_;
  end
  
  function z = zpts(S)
    z = complex(S.xk_, S.yk_);
  end
end

methods(Access=protected)
  function h = dumb_plot_(S, varargin)
    t = (0:200)'/200;
    h = plot(point(S, t), varargin{:});
  end
end

methods(Access=protected, Static)
  function [x, y] = get_pts_(tofignum)
    fprintf(['\n' ...
             '    Left mouse button picks points.\n' ...
             '    Right mouse button picks last point,\n' ...
             '      or <Enter> ends selection.\n' ...
             '    Point selection order determines curve\n' ...
             '      orientation.\n\n']);
    
    if isempty(tofignum)
      figure;
      axis(2*[-1 1 -1 1]);
      set(gca, 'dataaspectratio', [1 1 1]);
    else
      figure(tofignum);
    end
    hold on
    grid on
    box on
    
    x = [];
    y = [];
    sh = [];
    np = 0;
    button = 1;
    while button == 1
      [xi, yi, button] = ginput(1);
      if isempty(button)
        % return pressed
        break
      end
      plot(xi, yi, 'bo');
      np = np + 1;
      x(np) = xi; %#ok<AGROW>
      y(np) = yi; %#ok<AGROW>
      
      if ~isempty(sh)
        delete(sh);
      end
      if numel(x) > 1
        sh = dumb_plot_(splinep(x, y), 'k');
      end
      text(xi, yi, ['  ' int2str(np)], 'EraseMode', 'background');
    end
    
    x = [x(:); x(1)];
    y = [y(:); y(1)];
  end
  
  function [pp, tl] = spline_(x, y)
    % This algorithm is from " PERIODIC CUBIC SPLINE INTERPOLATION USING    
    % PARAMETRIC SPLINES" by W.D. Hoskins and P.R. King, Algorithm 73, The
    % Computer Journal, 15, 3(1972) P282-283. Fits a parametric periodic
    % cubic spline through n1 points (x(i), y(i)) (i = 1, ... ,n1) with 
    % x(1) = x(n1) and y(1) = y(n1). This function returns the first three
    % derivatives of x and y, the chordal distances h(i) of (x(i),y(i)) and
    % (x(i + 1), y(i + 1)) (i = 1, ..., n1 - 1) with h(n1) = h(1) and the
    % total distance.
    %
    % See also INTERP_1, INTERP_2, MATLAB function SPLINE.

    % Thomas K. DeLillo, Lianju Wang 07-05-99.
    % modified a bit by E. Kropf, 2013, 2014.
    
    if abs(x(1) - x(end)) > 100*eps || abs(y(1) - y(end)) > 100*eps
      x(end+1) = x(1);
      y(end+1) = y(1);
    end
    nk = length(x);
    n = nk - 1;
    dx = diff(x);
    dy = diff(y);
    h = sqrt(dx.^2 + dy.^2);
    tl = sum(h);
    h(nk) = h(1);
    p = h(1:n);
    q = h(2:nk);
    a = q./(p + q);
    b = 1 - a;
    c = spdiags(...
      [ [b(n);ones(n-1,1)] [a(2:n);0] 2*ones(n,1) [0;b(1:n-1)] ...
        [ones(n-1,1);a(1)] ], ...
      [-n+1 -1 0 1 n-1], n, n);
    d1 = 3*(a.*dx./p + b.*[dx(2:n); x(2) - x(nk)]./q);
    mmdflag = spparms('autommd');
    spparms('autommd', 0);
    x1 = c\d1;
    spparms('autommd', mmdflag);
    x1(2:nk) = x1;
    x1(1) = x1(nk);
    d = 3*(a.*dy./p + b.*[dy(2:n); y(2) - y(nk)]./q);
    mmdflag = spparms('autommd');
    spparms('autommd', 0);
    y1 = c\d;
    spparms('autommd', mmdflag);
    y1(2:nk) = y1;
    y1(1) = y1(nk);
    x2(2:nk) = 2*(x1(1:n) + 2*x1(2:nk) - 3*dx./p)./p;
    y2(2:nk) = 2*(y1(1:n) + 2*y1(2:nk) - 3*dy./p)./p;
    x2(1) = x2(nk);
    y2(1) = y2(nk);
    x2 = x2';
    y2 = y2';
    x3  = diff(x2)./p;
    x3(nk) = x3(1);
    y3 = diff(y2)./p;
    y3(nk) = y3(1);
    
    % Make pp for later evaluation.
    pp = cell(2, 3);
    t = [0; cumsum(h)];
    
    pp{1,1} = mkpp(t, [x3/6, x2/2, x1, x]);
    pp{2,1} = mkpp(t, [y3/6, y2/2, y1, y]);
    for j = 1:2
      coef = pp{j,1}.coefs;
      pp{j,2} = mkpp(t, [3*coef(:,1), 2*coef(:,2), coef(:,3)]);
      pp{j,3} = mkpp(t, [6*coef(:,1), 2*coef(:,2)]);
    end
  end
end

end

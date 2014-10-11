classdef mobius < conformalmap
% MOBIUS transformation class.
%   MOBIUS(Z,W) creates the Mobius transformation that maps the
%   3-vector Z to W. One infinity is allowed in each of Z and W.
%   
%   MOBIUS(a,b,c,d) creates the transformation
%   
%         a*z  +  b
%         ---------
%         c*z  +  d
%         
%   MOBIUS([a b; c d]) is also allowed. In either of these cases, a,b,c,d
%   should be finite complex numbers.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% (Re)written by Everett Kropf, 2014,
% adapted from code by Toby Driscoll, originally 20??.

properties
  matrix_
end

methods
  function M = mobius(varargin)
    domain = [];
    range = [];
    matrix = [];
    
    switch nargin        
      case 1
        A = varargin{1};
        if isa(A, 'double') && isequal(size(A), [2, 2])
          matrix = A;
        else
          error('CMT:InvalidArgument', ...
                'Single argument should be a 2-by-2 matrix.')
        end
        
      case 2
        [z, w] = deal(varargin{:});
        if isa(z, 'circle') && isa(w, 'circle')
          z = point(z, pi*[0.5, 1, 1.5]);
          w = point(w, pi*[0.5, 1, 1.5]);
        end
        if (isa(z, 'double') && length(z) == 3) && (isa(w, 'double') && ...
            length(w) == 3)
          A1 = mobius.standardmap(z);
          circ = circle(z);
          if ~isinf(circ)
            domain = disk(circ);
          end
          A2 = mobius.standardmap(w);
          circ = circle(w);
          if ~isinf(circ)
            range = disk(circle(w));
          end
          matrix = A2\A1;
        else
          error('CMT:InvalidArgument', ...
                'Invalid arguments; see help for mobius.')
        end
        
      case 4
        matrix = reshape(cat(1, varargin{:}), [2, 2]).';
    end
    
    if ~isempty(matrix) && rcond(matrix) < eps
      warning('Mobius map appears to be singular.')
    end
    
    if ~nargin
      supargs = {};
    else
      supargs = {domain, range};
    end
    M = M@conformalmap(supargs{:});
    M.matrix_ = matrix;
  end % ctor
    
  function out = char(map)
    % CHAR Pretty-print a Mobius map.
    
    %   Copyright (c) 1998-2006 by Toby Driscoll.

    % Numerator
    num = '';
    a = map.matrix_(1,1);
    if a~=0
      if a~=1
        if isreal(a)
          num = [num num2str(a,4) '*'];
        elseif isreal(1i*a)
          num = [num num2str(imag(a),4) 'i*'];
        else
          num = [num '(' num2str(a,4) ')*'];
        end
      end
      num = [num 'z'];
    end
    a = map.matrix_(1,2);
    if a~=0
      if ~isempty(num)
        s = sign(real(a));
        if s==0, s = sign(imag(a)); end
        if s > 0
          num = [num ' + '];
        else
          num = [num ' - '];
          a = -a;
        end
      end
      if isreal(a)
        num = [num num2str(a,4)];
      elseif isreal(1i*a)
        num = [num num2str(imag(a),4) 'i'];
      else
        num = [num '(' num2str(a,4) ')'];
      end
    end
    
    % Denominator
    den = '';
    a = map.matrix_(2,1);
    if a~=0
      if a~=1
        if isreal(a)
          den = [den num2str(a,4) '*'];
        elseif isreal(1i*a)
          den = [den num2str(imag(a),4) 'i*'];
        else
          den = [den '(' num2str(a,4) ')*'];
        end
      end
      den = [den 'z'];
    end
    a = map.matrix_(2,2);
    if a~=0
      if ~isempty(den)
        s = sign(real(a));
        if s==0, s = sign(imag(a)); end
        if s > 0
          den = [den ' + '];
        else
          den = [den ' - '];
          a = -a;
        end
      end
      if isreal(a)
        den = [den num2str(a,4)];
      elseif isreal(1i*a)
        den = [den num2str(imag(a),4) 'i'];
      else
        den = [den '(' num2str(a,4) ')'];
      end
    end
    
    L = [length(num),length(den)];
    D = (max(L)-L)/2;
    num = [blanks(floor(D(1))) num blanks(ceil(D(1)))];
    den = [blanks(floor(D(2))) den blanks(ceil(D(2)))];
    fline = repmat('-',1,max(L));
    
    out = sprintf('\n  %s\n  %s\n  %s\n',num,fline,den);
  end % char
  
  function disp(f)
    if isempty(f.matrix_)
      fprintf('\n\tempty transformation matrix\n\n')
    else
      disp(char(f))
    end
  end
  
  function w = feval(M, z)
    warning('mobius.feval() is depricated, use mobius.apply() instead.')
    w = apply_map(M, z);
  end
  
  function Minv = inv(M)
    % Inverse transformation.
    
    % Original code turned off builtin singular matrix warning and supplied
    % a repeat of warning supplied in constructor. Skipping here.
    Minv = mobius(inv(M.matrix_));
  end
  
  function A = matrix(M)
    A = M.matrix_;
  end
  
  function M = mrdivide(M1, M2)
    % Divide Mobius map by a scalar, or reciprocate it.
    %   1/M, for Mobius map M, swaps the numerator and denominator of M.
    %   M/c, for scalar c, multiplies the denominator of M by c.
    
    %   Copyright (c) 1998 by Toby Driscoll.
    %   $Id: mrdivide.m,v 1.1 1998/07/01 20:14:22 tad Exp $
    
    if isequal(M1, 1)
      % Exchange numerator and denominator
      A = M2.matrix_;
      M = mobius(A([2, 1],:));
    elseif isa(M2, 'double') && length(M2) == 1
      A = M1.matrix_;
      M = mobius([1, 0; 0, M2]*A);
    else
      error('Division not defined for these operands.')
    end
  end
  
  function M = mtimes(M1, M2)
    % Multiply Moebius transformation by a scalar.
    
    if isa(M1, 'double')
      % Make the first one mobius.
      tmp = M1;
      M1 = M2;
      M2 = tmp;
    elseif isa(M2, 'mobius')
      M = mobius(M1.matrix_*M2.matrix_);
      M.theDomain = domain(M2);
      M.theRange = range(M1);
      return
    elseif isa(M2, 'conformalmap')
      M = mtimes@conformalmap(M2, M1);
      return
    end
    
    A = M1.matrix_;
    if isa(M2, 'double') && length(M2) == 1
      A(1,:) = A(1,:)*M2;
      M = mobius(A);
    else
      error('CMT:NotDefined', 'Operation not defined.')
    end
  end
  
  function z = pole(M)
    % Return the pole of the Moebius map.
    A = M.matrix_;
    z = double(homog(-A(2,2)/A(2,1)));
  end
  
  function M = pretty(M, c)
    % Normalize a Moebius transformation for 'nicer' numbers.
    %   PRETTY(M), where M is a moebius map, divides the coefficients of M
    %   by the constant term in the demoninator, or (if that is zero) the
    %   coefficient of the linear term in the denominator. Then any real or
    %   imaginary parts that are close to machine precision are rounded to
    %   exactly zero.
    %
    %   PRETTY(M,C) instead normalizes all coefficients by C, then cleans up
    %   small numbers.
    
    
    %   Copyright (c) 2004, 2006 by Toby Driscoll.
    %   $Id$
    
    A = M.matrix_;
    
    % Normalization
    if nargin==1
      d = A(2,2);
      if d==0
        d = A(2,1);
      end
      A = A/d;
    else
      A = c*A;
    end
    
    % Rounding
    index = abs(imag(A)) < 100*eps;
    A(index) = real(A(index));
    index = abs(real(A)) < 100*eps;
    A(index) = 1i*imag(A(index));
    
    M = mobius(A);
  end % pretty
    
  function M = uminus(M)
    M = -1*M;
  end

  function z = zero(M)
    % Return the zero of the Mobius map.
    A = M.matrix_;
    z = double(homog(-A(1,2)/A(1,1)));
  end
end

methods(Access=protected)
  function w = apply_map(M, z)
    % Evaluate Mobius transformation.
    
    %   Copyright (c) 2006 by Toby Driscoll.

    switch(class(z))
      case {'circle', 'zline'}
        zp = pole(M);
        if dist(z, zp) < 10*eps(zp)
          % Result appears to be a line.
          zp = apply_map(M, point(z, [0.5, 1.5]*pi));
          w = zline(zp);
        else
          % Find new circle using three points.
          zp = apply_map(M, point(z, [0.5, 1, 1.5]*pi));
          w = circle(zp);
        end
        
      case 'double'
        w = NaN(size(z));
        
        % Convert inputs to homogeneous coordinates, and reshape
        z = homog(z);
        Z = [numer(z(:)).'; denom(z(:)).'];
        
        % Apply map
        W = M.matrix_*Z;
        
        % Convert to complex without DBZ warnings
        w(:) = double(homog(W(1,:), W(2,:)));
        
      otherwise
        error('Mobius maps can be applied to floats or circles/zlines only')
    end
  end
end

methods(Static)
  function A = standardmap(z)
    % Return mobius matrix for map from z(1:3) to [0, 1, Inf].
    if isinf(z(1))
      A = [0, z(2) - z(3); 1, -z(3)];
    elseif isinf(z(2))
      A = [1, -z(1); 1, -z(3)];
    elseif isinf(z(3))
      A = [1, -z(1); 0, z(2) - z(1)];
    else
      rms = z(2) - z(3);
      rmq = z(2) - z(1);
      A = [rms, -z(1)*rms; rmq, -z(3)*rmq];
    end
  end % standardmap
end

end
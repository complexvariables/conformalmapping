classdef homog
% HOMOG is the homogenous coordinate class.
%
% Replace complex number z by pair (z1, z2), such that z = z1/z2. Interacts
% quietly with builtin double data type.
%
% zeta = homog(z1, z2) -- If z2 is not supplied, we assume z2=1.

% This file is a part of the CMToolbox.
% It is licensed under the BSD 3-clause license.
% (See LICENSE.)

% Copyright Toby Driscoll, 2014.
% Written by Everett Kropf, 2014,
% adapted to new classdef from Toby Driscoll's code, originally 20??.

properties
  numer_
  denom_
end

methods
  function zeta = homog(z1, z2)
    % Constructor
    if nargin > 0
      if nargin < 2
        if isa(z1, 'homog')
          zeta = z1;
          return
        end
        z2 = [];
      end
      
      if ~isequal(size(z1), size(z2))
        if isempty(z2)
          % Assume 1 for denominator
          z2 = ones(size(z1));
        elseif numel(z2) == 1
          % Scalar expansion.
          z2 = repmat(z2, size(z1));
        else
          error('Input arguments must be scalar or of the same size.')
        end
      end
      
      % Transform infinities to finite representation. There is no unique
      % choice, so we arbtrarily pick [\pm 1,0] based on the sign.
      idx = isinf(z1);
      z1(idx) = sign(real(z1(idx))) + 1i*sign(imag(z1(idx)));
      z2(idx) = 0;
      
      zeta.numer_ = z1;
      zeta.denom_ = z2;
    end
  end % ctor
  
  function r = abs(zeta)
    % Return absolute value.
    r = abs(double(zeta));
  end
  
  function theta = angle(zeta)
    % Return phase angle, standardised to [-pi, pi).
    theta = mod(angle(zeta.numer_) - angle(zeta.denom_) + pi, 2*pi) - pi;
  end
  
  function zeta = cat(dim, varargin)
    % Override double cat().
    numers = cell(nargin - 1, 1);
    denoms = numers;
    for n = 1:nargin - 1
      h = homog(varargin{n});
      numers{n} = h.numer_;
      denoms{n} = h.denom_;
    end
    
    try
      zeta = homog(cat(dim, numers{:}), cat(dim, denoms{:}));
    catch
      error('Argument dimensions are not consistent.')
    end
  end
  
  function zetbar = conj(zeta)
    % Complex conjugate.
    zetbar = homog(conj(zeta.numer_), conj(zeta.denom_));
  end
  
  function eta = ctranspose(zeta)
    % Complex transpose.
    eta = homog(ctranspose(zeta.numer_), ctranspose(zeta.denom_));
  end
  
  function z2 = denom(zeta)
    % Return denominator.
    z2 = zeta.denom_;
  end
  
  function display(zeta)
    % Format for viewing pleasure.
    n = size(zeta.numer_);
    fprintf('\n\t%s array of homogenous coordinates:\n\n', ...
            [sprintf('%i-by-', n(1:end-1)), sprintf('%i', n(end))])
    fprintf('numerator = \n\n')
    disp(zeta.numer_)
    fprintf('\ndenominator = \n\n')
    disp(zeta.denom_)
  end
  
  function z = double(zeta)
    % Convert to double.
    
    % Driscoll's original turned of divide by zero warning. Do we still need
    % this? Newer versions of MATLAB don't give this warning.
    
    z = zeta.numer_./zeta.denom_;
    % Ensure imag(z(isinf(z))) = 0 reliably.
    z(isinf(z)) = Inf;
  end
  
  function e = end(zeta, k, n)
    % Return array end indexes.
    if n == 1
      e = length(zeta.numer_);
    else
      e = size(zeta.numer, k);
    end
  end
  
  function zeta = horzcat(varargin)
    % Provide horizontal contatenation.
    zeta = cat(2, varargin{:});
  end
  
  function z = inv(zeta)
    % Return 1/zeta.
    z = homog(zeta.denom_, zeta.numer_);
  end
    
  function tf = isinf(zeta)
    tf = zeta.denom_ == 0 & zeta.numer_ ~= 0; 
  end
  
  function n = length(zeta)
    % Length of zeta.
    n = length(zeta.numer_);
  end
  
  function c = minus(a, b)
    % Provide subtraction.
    c = plus(a, -b);
  end
  
  function c = mldivide(a, b)
    % Provide matrix left divide.
    c = mtimes(inv(a), b);
  end
  
  function c = mrdivide(a, b)
    % Provide matrix right divide.
    c = mtimes(a, inv(b));
  end
  
  function c = mtimes(a, b)
    % Provide multiplication.
    if isfloat(a)
      a = homog(a);
    end
    if isfloat(b)
      b = homog(b);
    end
    c = homog(a.numer_*b.numer_, a.denom_*b.denom_);
  end
  
  function n = numel(zeta, varargin)
    n = numel(zeta.numer_, varargin{:});
  end
  
  function z1 = numer(zeta)
    % Return numerator.
    z1 = zeta.numer_;
  end
  
  function c = plus(a, b)
    % Provide addition.
    if isfloat(a)
      a = homog(a);
    end
    if isfloat(b)
      b = homog(b);
    end
    c = homog(a.numer_*b.denom_ + a.denom_*b.numer_, a.denom_*b.denom_);
  end
  
  function zeta = subsref(zeta, s)
    % Provide double-like indexing.
    switch s.type
      case '()'
        zeta = homog(subsref(zeta.numer_, s), subsref(zeta.denom_, s));
      otherwise
        error('This type of indexing is not supported by homog objects.')
    end
  end
  
  function zeta = subsasgn(zeta, s, val)
    % Provide double-like assignment.
    switch s.type
      case '()'
        if length(s.subs) == 1
          zeta = homog(zeta);
          val = homog(val);
          index = s.subs{1};
          zeta.numer_(index) = val.numer_;
          zeta.denom_(index) = val.denom_;
        else
          error('HOMOG objects support linear indexing only.')
        end
      otherwise
        error('Unspported assignment syntax.')
    end
  end
  
  function eta = transpose(zeta)
    % Provide basic transpose.
    eta = homog(transpose(zeta.numer_), transpose(zeta.denom_));
  end
  
  function zeta = vertcat(varargin)
    % Provide vertical contatenation.
    zeta = cat(1, varargin{:});
  end
  
  function b = uminus(a)
    % Unitary minus.
    b = homog(-a.numer_, a.denom_);
  end
end

end
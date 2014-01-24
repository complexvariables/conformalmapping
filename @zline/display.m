function display(this)
% Pretty-print a line.

%   Copyright (c) 2006 by Toby Driscoll.

fprintf(['\nLine passing through ' num2str(this.base) ...
  ' and parallel to ' num2str(this.tangent) '\n\n'])

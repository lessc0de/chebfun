function I = eye(d)
%EYE       Identity operator.
%   This function is deprecated. Use OPERATORBLOCK.EYE.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org/ for Chebfun information.

I = linop( operatorBlock.eye(d) );

end

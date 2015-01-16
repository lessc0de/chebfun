function f = fracCumSum(f, m)
%FRACCUMSUM   Indefinite integral of a CHEBFUN for a fractional order.
%  FRACCUMSUM(F, M) is called by DIFF(F, M) and CUMSUM(F, M) when M is not 
%  integer and computes the fractional indefinite integral of order M of F, as 
%  defined by the Riemann–Liouville integral.
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Trivial case:
if ( isempty(f) )
    return
end

% Deal with the integer part of M:
for k = 1:m
    f = cumsum(f);
end

% The fractional part:
alph = m - floor(m);

% Grab some fields from f:
funs = f.funs;
numFuns = numel(funs);
dom = f.domain;

% Loop over each FUN by calling FRACCumSum@BNDFUN:
for k = 1:numFuns
    f.funs{k} = fracCumSum(f.funs{k}, alph, dom([1,end]));
end

% Simplify:
f = simplify(f);

% Grab new function values at the breakpoints:
pointValues = chebfun.getValuesAtBreakpoints(f.funs);

% Update the function value at the breakpoints:
f.pointValues = pointValues;

end
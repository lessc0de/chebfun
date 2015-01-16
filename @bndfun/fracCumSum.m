function f = fracCumSum(f, alph, dom)
%FRACUMSUM   Indefinite integral of a BNDFUN for a fractional order.
%  FRACUMSUM(FUNS, ALPH) returns the fractional indefinite integral of order 
%  FRACM defined by the piecewise smooth CLASSICFUNs stored in the cell FUNS and 
%  its domain is same as the CLASSICFUN given in the last entry of FUNS. Note 
%  that ALPH needs to be non-integer (fractional).
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

scl = (diff(f.domain)/2)^alph;

if ( nargin == 3 && f.domain(1) == dom(1) )
    % Cast to a SINGFUN:
    f.onefun = singfun(f.onefun);
end

f.onefun = scl * fracCumSum(f.onefun, f.mapping, alph);

end
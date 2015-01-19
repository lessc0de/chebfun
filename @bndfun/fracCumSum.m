function f = fracCumSum(f, alph)
%FRACUMSUM   Indefinite integral of a BNDFUN for a fractional order.
%  FRACUMSUM(FUNS, ALPH) returns the fractional indefinite integral of order 
%  ALPH for scalar ALPH in (0, 1).
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Scaling for domain:
scl = (diff(f.domain)/2)^alph;

if ( ~isa(f.onefun, 'singfun') )
    % Cast ONEFUN to a SINGFUN:
    f.onefun = singfun(f.onefun);
    % TODO: Do this here, or at ONEFUN level?
end

f.onefun = scl * fracCumSum(f.onefun, alph);

end
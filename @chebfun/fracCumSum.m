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


% Constant scaling:
scl = 1/gamma(alph);

% Grab some fields from f:
funs = f.funs;
numFuns = numel(funs);
dom = f.domain;

% For constructing the global contribution:
const = zeros(1, numFuns);
K = @(t, tau) (t-tau).^(alph-1);
data.exponents = [1-alph, 0];

% This shouldn't be needed, but I can't make the global contribution happy with
% out it. Very annyoing. NH 01/2015.
p = chebfunpref;
p.fixedLength = 65;

% Loop over each FUN by calling BNDFUN/FRACCUMSUM:
for k = 1:numFuns
    
    % Global contribution:
    I = 0;
    for j = 1:k-1 % Loop over all intervals before the current one:
        
        % Function which comuputes the fractional integral for a given t:
        op = @(t) scl*innerProduct(KFun(K, t, alph, dom(j:j+1)), funs{j}); 
        % Construct a bndfun for this in t.        
        data.domain = dom(k:k+1);            
        const(j) = vec(op, dom(k));  % Subtract const term so that Ij is of the 
        opj = @(t) op(t) - const(j); % form (1-x)^(1-alph)*<smooth>.
        Ij = bndfun(@(t) vec(opj, t), data, p);
        I = I + Ij;
    end
    
    % Local contribution:
    J = fracCumSum(funs{k}, alph);
    
    % Sum global and local contributions:
    %  (TODO: Except for the first piece -- in which case I = 0 -- (I + J) 
    %   *should* simplify to have zero exponents but often doesn't. Grrr.)
    warnState = warning('off', 'CHEBFUN:SINGFUN:plus:exponentDiff');
    f.funs{k} = (I + J) + sum(const); % Replace the constant terms.
    warning(warnState);

end

% Simplify:
f = simplify(f);

% Update the function value at the breakpoints:
f.pointValues = chebfun.getValuesAtBreakpoints(f.funs);

end

function y = KFun(K, t, alph, dom)
%KFUN  Construct a BNDFUN of the kernel for a given t, alph and domain.
    data.domain = dom;
    data.exponents = [0 0];
    if ( t == dom(2) )
        data.exponents = [0, alph-1];
    end
    y = bndfun(@(tau) K(t, tau), data);
end

function y = vec(op, t)
%VEC  Vectorize a given function handle.
    y = zeros(size(t));
    for l = 1:numel(t)
        y(l,:) = op(t(l));
    end
end

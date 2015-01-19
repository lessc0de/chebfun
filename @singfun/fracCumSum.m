function f = fracCumSum(f, alph)
%FRACCUMSUM   Fractional indefinite integral of a SINGFUN.
%  FRACCSUMSUM(ONEFUNS, ALPH) returns the fractional indefinite integral of
%  order ALPH for scalar ALPH in (0, 1).
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Exponents:
oldExps = f.exponents;
newExps = oldExps;
newExps(1) = newExps(1) + alph;

% Call the SINGFUN constructor to construct the Riemann–Liouville fractional 
% integral:
data.exponents = newExps;
f = (1/gamma(alph))*singfun(@(x) op(f, x, alph), data);

end

%% The operator of the Riemann–Liouville fractional integral:
function y = op(g, x, alph)

% Preallocation and vectorization:
sx = size(x);
x = x(:);
lx = length(x);
y = zeros(lx, 1);

% Loop over each sample point:
for k = 1:lx
    
    if ( x(k) == -1 )
        % When sampled at -1, namely the left endpoint of the entire CHEBFUN or 
        % the left edge of a piece, set the value of the fractional integral 
        % at -1 as NaN and this value will be extrapolated during construction:
        y(k) = NaN;
        continue
    end

    % Kernel mapped to [-1 1]:
    data.exponents = [0, alph-1];
    scl = (.5*(x(k)+1)).^alph;
    kernel = scl*singfun(@(s) (1-s).^(alph-1), data);
    
    if ( x(k) == 1 )
        % At x(k) = 1 (i.e., rightmost pt) we do not need to restrict.
        gk = g;
    else
        % Restrict g to [-1, x(k)]:
        gk = restrict(g, [-1, x(k)]);
    end
    
    % Integrate with innerProduct():
    y(k) = innerProduct(kernel, gk);
    
end

% Reshape the result:
y = reshape(y, sx);

end

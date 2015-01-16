function f = fracCumSum(f, map, alph)
%FRACCUMSUM   Fractional indefinite integral of a ONEFUN.
%  FRACCSUMSUM(ONEFUNS, MAPS, FRACM) returns the fractional indefinite integral
%  of order FRACM defined by the ONEFUNs stored in the cell ONEFUNS and their 
%  corresponding map structure which maps between the real domain and [-1 1]. 
%  Note that FRACM needs to be non-integer (fractional).
%
% See also CUMSUM, DIFF.

% Copyright 2014 by The University of Oxford and The Chebfun Developers.
% See http://www.chebfun.org for Chebfun information.

% Grab the exponent of the last piece:
exps = f.exponents;

% Exponents for the result:
newExps = exps;

newExps(1) = newExps(1) + alph;

% Call the SINGFUN constructor to construct the Riemann–Liouville fractional 
% integral:

data.exponents = newExps;
f = (1/gamma(alph))*singfun(@(x) op(f, map, x, alph), data);

end

%% The operator of the Riemann–Liouville fractional integral:
function y = op(g, map, x, alph)

% Get the smooth part and the exponents of the last piece:
exps = g.exponents;
sp = g.smoothPart;

% Preallocation and vectorization:
sx = size(x);
x = x(:);
lx = length(x);
y = zeros(lx, 1);

warnState = warning('off', 'MATLAB:integral:MinStepSize');

% Loop over each sample point:
for k = 1:lx
    t = x(k);
    
    


%     t = map.Inv(map.For(x(k)))
%     data.exponents = [0 0];
%     kernel = singfun(@(s) (t-s).^(alph-1), data);
% 
%     I = sum(kernel.*g);
    
    if ( x(k) == -1 || x(k) == 0 )
        % When sampled at -1, namely the left endpoint of the entire CHEBFUN or 
        % the left edge of a piece, set the value of the fractional integral 
        % at -1 as NaN and this value will be extrapolated during construction:
        y(k) = NaN;
    else

        y(k) = integral(@(tau) (t-tau).^(alph-1).*feval(g, tau), -1, t, ...
        'abstol', eps, 'reltol', eps);
        
%     elseif ( x(k) == 1 )
%         % Call SINGFUN constructor:
%         data.exponents = [exps(1) exps(2)+alph-1];
%         h = singfun(sp, data);
%         y(k) = I + sum(h);
%     else
%         % Call SINGFUN constructor:
%         rsp = restrict(sp, [-1 x(k)]);
%         data.exponents = [exps(1) alph-1];
%         h = singfun(rsp, data);
%         y(k) = I + ((x(k)+1)/2)^(alph + exps(1))*sum(h);
    end
    
end

warning(warnState);

% Reshape the result:
y = reshape(y, sx);

end
% Test file for fractional calculas.

% function pass = test_fracCalc(pref)

% Obtain preferences.
% if ( nargin == 0 )
    pref = chebfunpref();
% end

% Generate a few random points in [-1 1] to use as test values.
seedRNG(7681);
xr = sort(2 * rand(1000, 1) - 1);

%% Differentiation:
f = chebfun(@(x) sin(x), [-1 1]);
g = diff(f, 0.3);
h = diff(g, 0.7);
v = diff(f);
pass(1) = norm(h-v, inf) < 1e3*epslevel(v)*hscale(v);

%% Differentiation of a piecewise CHEBFUN:
f = chebfun({@(x) sin(x) @(x) sin(x)}, [-1 0 1]);
g = diff(f, 0.3);
h = diff(g, 0.7);
v = diff(f);
pass(2) = norm(h-v, inf) < 1e4*epslevel(v)*hscale(v);

%% Integration:
f = chebfun(@(x) sin(x), [-1 1]);
g = cumsum(f, 0.3);
h = cumsum(g, 0.7);
v = cumsum(f);
pass(3) = norm(h-v, inf) < 1e1*epslevel(v)*hscale(v);

%% Integration of a piecewise smooth CHEBFUN:
f = chebfun({@(x) sin(x) @(x) sin(x)}, [-1 0 1]);
g = cumsum(f, 0.3);
h = cumsum(g, 0.7);
v = cumsum(f);
pass(4) = norm(h-v, inf) < 1e3*epslevel(v)*hscale(v);

%% Differentiation of a singular function:
f = chebfun(@(x) sin(x).*((1+x).^0.3), [-1 1], 'exps', [0.3 0]);
g = diff(f, 0.3);
h = diff(g, 0.7);
v = diff(f);
err = feval(h-v, xr);
pass(5) = norm(err, inf) < 1e5*epslevel(v)*hscale(v);

%% Differentiation of a singular function represented by a piecewise smooth CHEBFUN:
f = chebfun({@(x) sin(x).*((1+x).^0.3) @(x) sin(x).*((1+x).^0.3)}, [-1 0 1], 'exps', [0.3 0 0 0]);
g = diff(f, 0.3);
h = diff(g, 0.7);
v = diff(f);
err = feval(h-v, xr);
pass(6) = norm(err, inf) < 1e4*epslevel(v)*hscale(v);

%% Integration of of a singular function:
f = chebfun(@(x) sin(x).*((1+x).^0.3), [-1 1], 'exps', [0.3 0]);
g = cumsum(f, 0.3);
h = cumsum(g, 0.7);
v = cumsum(f);
err = feval(h-v, xr);
pass(7) = norm(err, inf) < 1e1*epslevel(v)*hscale(v);

%% Integration of a singular function represented by a piecewise smooth CHEBFUN:
f = chebfun({@(x) sin(x).*((1+x).^0.3) @(x) sin(x).*((1+x).^0.3)}, [-1 0 1], 'exps', [0.3 0 0 0]);
g = cumsum(f, 0.3);
h = cumsum(g, 0.7);
v = cumsum(f);
err = feval(h-v, xr);
pass(8) = norm(err, inf) < 1e1*epslevel(v)*hscale(v);

%% Caputo definition:
f = chebfun(@(x)exp(x), [-2 7]);
g = diff(f, 0.3, 'caputo');
h = diff(g, 0.7, 'caputo');
v = diff(f);
pass(9) = norm(h-v, inf) < 2e3*epslevel(v)*hscale(v);

%% V4 tests.

tol = 1e-10;
dom = [0, 1];

% Polynomials
x = chebfun('x', dom);
q = sqrt(2)/2;
k = 9;
for n = [1 4]
    k = k + 1;
    xn = x.^n;
    xnpq = diff(xn, q, 'Caputo');
    tru = gamma(n+1)./gamma(n+1-q)*chebfun(@(x) x.^(n-q), dom, 'exps', [n-q, 0]); 
    pass(k) = norm(tru - xnpq, inf) < tol;
end

% Exponential
u = chebfun('exp(x)', dom);
trueRL = chebfun('erf(sqrt(x)).*exp(x) + 1./sqrt(pi*x)', dom, 'exps', [-.5 0]);
trueC = chebfun('erf(sqrt(x)).*exp(x)', dom, 'exps', [.5, 0]);

% RL
up05 = diff(u, .5, 'RL');
uint05 = cumsum(u, .5);
xr = xr(xr > 0);
pass(12) = norm(feval(trueRL, xr) - feval(up05, xr), inf) < 50*tol;
% Answer for int is actually same as C.
pass(13) = norm(trueC - uint05, inf) < tol; 

% Caputo
up05 = diff(u, .5, 'Caputo');
trueC = chebfun('erf(sqrt(x)).*exp(x)', dom, 'exps', [.5, 0]);
pass(14) = norm(trueC - up05, inf) < tol;

% end
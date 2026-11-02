clear; clc;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), ...
        fullfile(root,'analysis'),   fullfile(root,'test'));

p = params();

% --- TEST 2a: f = 1, control values ---
[xe, flag] = equilibrium(p, 1);
assert(flag > 0, 'fsolve did not converge')
assert(abs(xe(4) - 2.647)    < 5e-4, 'cda off')
assert(abs(xe(5) - 79.02)    < 5e-3, 'vda off')
assert(abs(xe(9) - 0.001975) < 5e-7, 'eda off')
assert(abs(g_auto(xe(9), p.eds) - 1) < 5e-5, 'g should be 1')
fprintf('TEST 2a (f=1) passed\n')

% --- TEST 2b: sweep in f, the seed is the previous solution ---
flist = [1 0.5 0.2 0.1 0.05 0.01];
x0 = xe;
for k = 1:numel(flist)
    fk = flist(k);
    [xk, flag] = equilibrium(p, fk, x0);
    res = max(abs(rhs(0, xk, p, fk)));
    assert(flag > 0,        'no convergence at f=%g', fk)
    assert(res < 1e-8,      'residual too big at f=%g: %.1e', fk, res)
    assert(all(xk > -1e-9), 'negative state at f=%g', fk)
    fprintf('f=%5.2f  cda=%7.4f  vda=%8.3f  eda=%.6f  |F|max=%.1e\n', ...
            fk, xk(4), xk(5), xk(9), res);
    x0 = xk;
end
fprintf('TEST 2b (sweep) passed\n')
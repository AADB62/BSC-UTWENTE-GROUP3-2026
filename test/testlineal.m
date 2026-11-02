% test_lineal.m
% Checks lineal.m: (1) A vs central finite differences, (2) A real and finite,
% (3) local stability (eigenvalues of A).
clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), ...
        fullfile(root,'analysis'),   fullfile(root,'test'));

p     = params();
fList = [1 0.5 0.2 0.1 0.05 0.01];

for f = fList
    [A, ~, ~, ~, xe] = lineal(p, f);

    % (1) Reference: central finite differences (step relative to each state)
    AFD = zeros(9);
    for k = 1:9
        hk = 1e-6*max(abs(xe(k)), 1e-3);
        e  = zeros(9,1);  e(k) = hk;
        AFD(:,k) = (rhs(0, xe+e, p, f) - rhs(0, xe-e, p, f)) / (2*hk);
    end
    err = max(abs(A(:) - AFD(:))) / max(abs(AFD(:)));

    % (2) The result must be real and finite
    ok = isreal(A) && all(isfinite(A(:)));

    % (3) Stability: eigenvalues of A (h^-1)
    lam = eig(A);

    fprintf('f = %5.2f | max rel err vs FD = %.2e | real/finite = %d | max Re(lambda) = %.3e h^-1\n', ...
            f, err, ok, max(real(lam)));
end
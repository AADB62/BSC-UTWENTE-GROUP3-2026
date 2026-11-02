% run_prcc.m  -  Global sensitivity: Latin Hypercube Sampling + PRCC (no toolbox needed)
%   PRCC = partial rank correlation coefficient between each parameter and each output,
%   controlling for all the other parameters (rank-transformed data).
clear; clc; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), ...
        fullfile(root,'analysis'),   fullfile(root,'test'));
p = params();

names    = setdiff(fieldnames(p), {'c'}, 'stable');     % c = unit conversion
k        = numel(names);
N        = 500;                                          % samples per scenario
D        = 100;                                          % mg
fList    = [1 0.5 0.2 0.1 0.05 0.01];
outNames = {'eda_eq','peak','t_peak','AUC','t_half'};
rng(1);                                                  % reproducible

% ---- ranges: [nominal/2, nominal*2], log-uniform; f5 is a fraction -> [0.35, 1] ----
nom = cellfun(@(n) p.(n), names);                        % k x 1
lo  = nom/2;   hi = nom*2;
i5  = strcmp(names, 'f5');   lo(i5) = 0.35;   hi(i5) = 1;

% ---- Latin Hypercube in [0,1]: one value in each of N strata, shuffled per parameter ----
U = zeros(N, k);
for j = 1:k
    U(:,j) = (randperm(N)' - rand(N,1)) / N;
end
X = lo.' .* (hi.'./lo.').^U;                             % N x k, log-uniform in [lo, hi]

% ---- run the model on every sample and scenario ----
Y = nan(N, numel(outNames), numel(fList));
for jf = 1:numel(fList)
    tic;
    for s = 1:N
        ps = p;
        for q = 1:k, ps.(names{q}) = X(s,q); end
        Y(s,:,jf) = outputs(ps, fList(jf), D);
        if mod(s,100) == 0
            fprintf('f = %.2f: %d/%d samples (%.0f s)\n', fList(jf), s, N, toc);
        end
    end
    fprintf('f = %.2f: failed samples = %d of %d\n', fList(jf), nnz(isnan(Y(:,2,jf))), N);
end

% ---- PRCC ----
R   = nan(k, numel(outNames), numel(fList));
SIG = false(k, numel(outNames), numel(fList));
for jf = 1:numel(fList)
    for m = 1:numel(outNames)
        [R(:,m,jf), SIG(:,m,jf)] = prcc(X, Y(:,m,jf));
    end
end

% ---- print: ten largest |PRCC| for peak and eda_eq in each scenario (* = significant) ----
for jf = 1:numel(fList)
    for m = [2 1]
        [~, ord] = sort(abs(R(:,m,jf)), 'descend');
        fprintf('\nPRCC of %s, D = %d mg, f = %.2f\n', outNames{m}, D, fList(jf));
        for q = ord(1:10)'
            fprintf('  %-7s %7.3f %s\n', names{q}, R(q,m,jf), repmat('*',1,SIG(q,m,jf)));
        end
    end
end

% ---- heat map of PRCC_peak (parameter x f) ----
cmap = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ...
        ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];
fig = figure('Name','prcc_peak');
imagesc(squeeze(R(:,2,:)));  colormap(cmap);  caxis([-1 1]);  colorbar
set(gca, 'XTick',1:numel(fList), 'XTickLabel', arrayfun(@(f) sprintf('%.2f',f), fList, 'UniformOutput',false), ...
         'YTick',1:k, 'YTickLabel', names, 'FontSize', 7);
xlabel('f');  title(sprintf('PRCC of the peak, D = %d mg', D));

% ---- monotonicity check: scatter of the peak vs the 4 top parameters at f = 0.01 ----
[~, ord] = sort(abs(R(:,2,end)), 'descend');
fig = figure('Name','prcc_scatter');
tiledlayout(1, 4, 'TileSpacing','compact');
for q = 1:4
    nexttile;
    loglog(X(:,ord(q)), Y(:,2,end)*1e3, '.');
    xlabel(names{ord(q)});  ylabel('peak \delta eda (nM)');  grid on
end

save(fullfile(root,'prcc.mat'), 'X','Y','R','SIG','names','fList','D','N','outNames','lo','hi');
exportgraphics(gcf, fullfile(root,'figures','sens_prcc_heatmap.pdf'), 'ContentType','vector');

% ================= local functions =================
function y = outputs(p, f, D)
    y = nan(1,5);
    try
        xe     = equilibrium(p, f);
        [t, x] = simulate(p, f, D, (0:0.02:24)');
        m      = metrics_eda(t, x, xe, Inf);
        y      = [xe(9), m.peak, m.t_peak, m.auc, m.t_half];
    catch
    end
end

function [r, sig] = prcc(X, y)
    ok = all(isfinite(X),2) & isfinite(y);               % drop failed samples
    X = X(ok,:);  y = y(ok);
    kk = size(X,2);  n = size(X,1);
    C  = corrcoef([rankcol(X) rankcol(y)]);
    P  = inv(C);                                         % precision matrix
    d  = diag(P);
    r  = -P(1:kk, end) ./ sqrt(d(1:kk)*d(end));          % partial correlation with y, given the rest
    dof = n - kk - 1;
    t   = r .* sqrt(dof ./ (1 - r.^2));
    sig = abs(t) > 1.96;                                 % ~5 % two-sided, normal approximation
end

function Rk = rankcol(A)
    [n, m] = size(A);  Rk = zeros(n, m);
    for j = 1:m
        [~, idx] = sort(A(:,j));
        Rk(idx, j) = (1:n)';
    end
end
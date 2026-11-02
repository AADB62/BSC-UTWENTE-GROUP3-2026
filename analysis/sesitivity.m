% run_sensitivity.m  -  OAT sensitivity over a grid of doses and surviving fractions
clear; clc; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), ...
        fullfile(root,'analysis'),   fullfile(root,'test'));
p = params();

DList    = [10 100 250];                 % mg  (10 mg: near-linear regime; 100, 250: clinical range)
fList    = [1 0.5 0.2 0.1 0.05 0.01];
delta    = 0.10;                         % +-10 %
outNames = {'eda_eq','peak','t_peak','AUC','t_half'};
names    = setdiff(fieldnames(p), {'c'}, 'stable');    % must match sensitivity_oat.m
nP = numel(names);  nO = numel(outNames);  nD = numel(DList);  nF = numel(fList);

S4 = nan(nP, nO, nD, nF);                % S4(param, output, dose, f)
Y0 = nan(nO, nD, nF);                    % nominal outputs
nTot = nD*nF;  n = 0;
for i = 1:nD
    for j = 1:nF
        tic;
        [S, ~, y0] = sensitivity_oat(p, fList(j), DList(i), delta);
        S4(:,:,i,j) = S;  Y0(:,i,j) = y0;
        n = n + 1;
        fprintf('scenario %2d/%d  (D = %3d mg, f = %.2f)  %.0f s  | remaining ~ %.1f min\n', ...
                n, nTot, DList(i), fList(j), toc, (nTot-n)*toc/60);
    end
end
fprintf('NaN entries (failed runs): %d\n', nnz(isnan(S4)));

% ---- 1. Ranking: max |S| over all scenarios, for each output ----
for m = 1:nO
    worst = squeeze(max(max(abs(S4(:,m,:,:)), [], 3), [], 4));   % nP x 1
    [~, ord] = sort(worst, 'descend');
    fprintf('\nOutput %s: ten most influential parameters (max |S| over all scenarios)\n', outNames{m});
    for k = ord(1:10)'
        fprintf('  %-7s %7.3f\n', names{k}, worst(k));
    end
end

% ---- 2. Non-influential parameters: |S| < 0.05 for every output and scenario ----
worstAll = squeeze(max(max(max(abs(S4), [], 2, 'omitnan'), [], 3), [], 4));
fprintf('\nParameters with max |S| < 0.05 for ALL outputs and scenarios (candidates to fix):\n  ');
fprintf('%s ', names{worstAll < 0.05});  fprintf('\n');

% ---- 3. S_peak vs f for each dose (top 10 parameters by max |S_peak|) ----
worstPk = squeeze(max(max(abs(S4(:,2,:,:)), [], 3), [], 4));
[~, ordPk] = sort(worstPk, 'descend');  top = ordPk(1:10);
for i = 1:nD
    fprintf('\nS_peak, D = %d mg\n', DList(i));
    fprintf('%-7s', 'param');  fprintf(' f=%-6.2f', fList);  fprintf('\n');
    for k = top'
        fprintf('%-7s', names{k});
        fprintf(' %8.3f', squeeze(S4(k,2,i,:)));  fprintf('\n');
    end
end

% ---- 4. Heat maps of S_peak (parameter x f), one panel per dose ----
cmap = [linspace(0,1,128)' linspace(0,1,128)' ones(128,1); ...
        ones(128,1) linspace(1,0,128)' linspace(1,0,128)'];     % blue-white-red
lim  = max(abs(S4(:,2,:,:)), [], 'all', 'omitnan');
fig = figure('Name','sensitivity_heatmap');
tiledlayout(1, nD, 'TileSpacing','compact');
for i = 1:nD
    nexttile;
    imagesc(squeeze(S4(:,2,i,:)));
    colormap(cmap);  caxis([-lim lim]);
    set(gca, 'XTick',1:nF, 'XTickLabel', arrayfun(@(f) sprintf('%.2f',f), fList, 'UniformOutput',false), ...
             'YTick',1:nP, 'YTickLabel', names, 'FontSize', 7);
    xlabel('f');  title(sprintf('S_{peak}, D = %d mg', DList(i)));
    if i == nD, colorbar; end
end

save(fullfile(root,'sensitivity.mat'), 'S4','Y0','names','DList','fList','delta','outNames');
exportgraphics(gcf, fullfile(root,'figures','sens_local_heatmap.pdf'), 'ContentType','vector');

function [S, names, y0] = sensitivity_oat(p, f, D, delta)
% SENSITIVITY_OAT  Local one-at-a-time normalised sensitivity of the outputs to each parameter.
%   S(k,m) = [Y_m(p_k(1+delta)) - Y_m(p_k(1-delta))] / (2*delta*Y_m(p))
%   Outputs m: 1 eda_eq, 2 peak, 3 t_peak, 4 AUC, 5 t_half
%   p_k is perturbed alone; the equilibrium is recomputed for each perturbed value.

    names = setdiff(fieldnames(p), {'c'}, 'stable');   % c = unit conversion, not a parameter
    y0    = outputs(p, f, D);
    S     = nan(numel(names), numel(y0));

    for k = 1:numel(names)
        pp = p;  pm = p;
        pp.(names{k}) = p.(names{k})*(1 + delta);
        pm.(names{k}) = p.(names{k})*(1 - delta);
        yp = outputs(pp, f, D);
        ym = outputs(pm, f, D);
        S(k,:) = (yp - ym) ./ (2*delta*y0);
    end
end

function y = outputs(p, f, D)
    y = nan(1,5);
    try
        xe     = equilibrium(p, f);
        tspan  = (0:0.01:24)';                  % fine grid: t_peak resolution 0.01 h
        [t, x] = simulate(p, f, D, tspan);
        m      = metrics_eda(t, x, xe, Inf);    % theta not needed here
        y      = [xe(9), m.peak, m.t_peak, m.auc, m.t_half];
    catch
        % equilibrium or integration failed: leave NaN
    end
end
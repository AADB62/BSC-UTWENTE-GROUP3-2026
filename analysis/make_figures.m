% make_figures.m
% Report figures from the nonlinear model (reads results.mat from main_simulation.m).
clear; clc; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), ...
        fullfile(root,'analysis'),   fullfile(root,'test'));
load(fullfile(root,'results.mat'), 'res', 'DList', 'fList', 'theta');
p = params();

outDir = fullfile(root,'figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
set(groot, 'defaultAxesFontSize', 11, 'defaultLineLineWidth', 1.5);
col  = parula(numel(fList)+1);  col = col(1:end-1,:);        % one colour per f (no yellow)
fSub = [1 0.5 0.2 0.1 0.05 0.01];                            % curves that appear in the legends
savefig_pdf = @(fig,name) exportgraphics(fig, fullfile(outDir,[name '.pdf']), 'ContentType','vector');
iD  = @(D) find(DList == D, 1);                              % index of a dose

% ---- Fig 1: absolute eda(t), 100 mg ----
fig = figure('Name','fig_eda_absolute');  hold on
i = iD(100);
for j = 1:numel(fList)
    hv = 'off';  if ismember(fList(j), fSub), hv = 'on'; end
    plot(res(i,j).t, res(i,j).x(:,9)*1e3, 'Color', col(j,:), ...
         'DisplayName', sprintf('f = %.2f', fList(j)), 'HandleVisibility', hv);
end
yline(theta*1e3, 'k--', 'HandleVisibility', 'off');          % threshold (Reed rule)
text(23.8, theta*1e3, sprintf('\\theta = %.2f nM', theta*1e3), ...
     'HorizontalAlignment','right', 'VerticalAlignment','bottom');
xlabel('t (h)');  ylabel('eda (nM)');  xlim([0 24]);  grid on
legend('Location', 'northeast');
title('Extracellular dopamine, D = 100 mg');
savefig_pdf(fig, 'fig1_eda_absolute');

% ---- Fig 2: delta eda(t), 100 and 250 mg ----
fig = figure('Name','fig_delta_eda');
tl = tiledlayout(1, 2, 'TileSpacing', 'compact');
for D = [100 250]
    nexttile;  hold on
    i = iD(D);
    for j = 1:numel(fList)
        hv = 'off';  if ismember(fList(j), fSub), hv = 'on'; end
        plot(res(i,j).t, (res(i,j).x(:,9) - res(i,j).xe(9))*1e3, 'Color', col(j,:), ...
             'DisplayName', sprintf('f = %.2f', fList(j)), 'HandleVisibility', hv);
    end
    xlabel('t (h)');  ylabel('\delta eda (nM)');  xlim([0 24]);  grid on
    title(sprintf('D = %d mg', D));
end
legend('Location', 'northeast');                             % last tile only
savefig_pdf(fig, 'fig2_delta_eda');

% ---- Fig 3: metrics vs f ----
fig = figure('Name','fig_metrics');
tl = tiledlayout(2, 2, 'TileSpacing', 'compact');
names  = {'peak','auc','t_half','rel_rise'};
ylabs  = {'peak \delta eda (nM)','AUC_{0-24 h} (nM\cdoth)','half-peak duration (h)','peak / eda_{eq} (-)'};
scales = [1e3 1e3 1 1];
for k = 1:4
    nexttile;  hold on
    for D = [100 250]
        i = iD(D);
        v = arrayfun(@(j) res(i,j).m.(names{k}), 1:numel(fList)) * scales(k);
        plot(fList, v, '.-', 'MarkerSize', 12, 'DisplayName', sprintf('%d mg', D));
    end
    set(gca, 'XDir', 'reverse', 'XScale', 'log');
    xlabel('f (surviving fraction)');  ylabel(ylabs{k});  grid on
    if k == 1, legend('Location','northwest'); end
end
savefig_pdf(fig, 'fig3_metrics_vs_f');

% ---- Fig 4: baseline eda at equilibrium vs f ----
fine = logspace(-2, 0, 40);
eq   = zeros(size(fine));
for k = 1:numel(fine)
    xe = equilibrium(p, fine(k));
    eq(k) = xe(9)*1e3;
end
fig = figure('Name','fig_baseline');  hold on
plot(fine, eq, 'k-');
plot(fList, arrayfun(@(j) res(1,j).xe(9)*1e3, 1:numel(fList)), 'ro', 'MarkerFaceColor','r');
yline(theta*1e3, '--', sprintf('\\theta = %.2f nM', theta*1e3), ...
      'LabelHorizontalAlignment','right', 'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
set(gca, 'XDir', 'reverse', 'XScale', 'log');
xlabel('f (surviving fraction)');  ylabel('eda_{eq} (nM)');  grid on
legend('equilibrium','simulated f','Location','southwest');
savefig_pdf(fig, 'fig4_baseline_vs_f');

% ---- Fig 5: dose-response (peak vs D), log-log, selected f ----
doses = DList(DList > 0);
fig = figure('Name','fig_dose_response');  hold on
for j = find(ismember(fList, fSub))
    pk = arrayfun(@(D) res(iD(D),j).m.peak, doses)*1e3;
    loglog(doses, pk, 'o-', 'Color', col(j,:), 'DisplayName', sprintf('f = %.2f', fList(j)));
end
i1  = iD(doses(1));                                          % smallest dose: linear reference
pk1 = res(i1,1).m.peak*1e3;                                  % f = 1 is fList(1)
loglog(doses, pk1*doses/doses(1), 'k--', 'DisplayName', 'linear (f = 1)');
set(gca, 'XScale','log', 'YScale','log');
xlabel('dose D (mg)');  ylabel('peak \delta eda (nM)');  grid on
legend('Location','northwest');
title('Dose-response (1800 mg: non-clinical, diagnostic)');
savefig_pdf(fig, 'fig5_dose_response');

% ---- Fig 6: time constants vs f (all modes, local stability) ----
fig = figure('Name','fig_tau');  hold on
for j = 1:numel(fList)
    lam = eig(lineal(p, fList(j)));
    tau = -60 ./ real(lam);                                  % min  (lambda in 1/h)
    semilogy(fList(j)*ones(size(tau)), tau, 'k.', 'MarkerSize', 14);
end
set(gca, 'XDir','reverse', 'XScale','log', 'YScale','log');
xlabel('f (surviving fraction)');  ylabel('\tau = -1/Re(\lambda) (min)');  grid on
title('Time constants of the 9 modes (all Re \lambda < 0)');
savefig_pdf(fig, 'fig6_time_constants');
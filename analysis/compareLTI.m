% compare_lti.m
% Compares the nonlinear response with the local LTI approximation:
%   nonlinear : res(i,j) from main_simulation.m (ode15s)
%   linear    : d(eda)(t) = D * C * expm(A t) * B   (linear in D, computed once per mg)
% Output: error table, waveform figures, peak error vs f, and error vs dose (tangent test).
clear; clc; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), ...
        fullfile(root,'analysis'),   fullfile(root,'test'));
load(fullfile(root,'results.mat'), 'res', 'DList', 'fList', 'theta');
p = params();

outDir = fullfile(root,'figures');
if ~exist(outDir,'dir'), mkdir(outDir); end
set(groot, 'defaultAxesFontSize', 11, 'defaultLineLineWidth', 1.5);
savefig_pdf = @(fig,name) exportgraphics(fig, fullfile(outDir,[name '.pdf']), 'ContentType','vector');

fSub  = [1 0.5 0.2 0.1 0.05 0.01];          % f shown in tables and legends
dt    = 0.01;                               % h, uniform grid for the linear response
tl    = (0:dt:24)';
N     = numel(tl);
doses = find(DList > 0);
nD = numel(doses);  nF = numel(fList);
col = parula(nF+1);  col = col(1:end-1,:);

peak_err = nan(nD,nF);  auc_err = nan(nD,nF);  nrmse = nan(nD,nF);
tpk_nl   = nan(nD,nF);  tpk_li  = nan(nD,nF);
dl_all = cell(nD,nF);   dn_all = cell(nD,nF);
pk_nl  = nan(nD,nF);    pk_li   = nan(nD,nF);

for j = 1:nF
    f = fList(j);
    [A, B, C, ~, xe] = lineal(p, f);        % LTI at x_e(f); independent of the dose

    % ---- unit response (1 mg): r(t) = C expm(A t) B on the uniform grid ----
    Ad = expm(A*dt);                        % one-step propagator
    xk = B;
    r  = zeros(N,1);
    for k = 1:N
        r(k) = C*xk;                        % output at t_k
        xk   = Ad*xk;                       % advance one step
    end

    for ii = 1:nD
        i = doses(ii);
        D = DList(i);
        dl = D*r;                           % linear response to D mg

        % ---- reuse metrics_eda on both responses ----
        xl = zeros(N, 9);  xl(:,9) = xe(9) + dl;
        ml = metrics_eda(tl, xl, xe, theta);
        mn = metrics_eda(res(i,j).t, res(i,j).x, res(i,j).xe, theta);

        % ---- waveform error on the common grid ----
        dn = interp1(res(i,j).t, res(i,j).x(:,9) - res(i,j).xe(9), tl);

        peak_err(ii,j) = 100*(ml.peak - mn.peak)/mn.peak;
        auc_err(ii,j)  = 100*(ml.auc  - mn.auc )/mn.auc;
        nrmse(ii,j)    = 100*norm(dl - dn)/norm(dn);
        tpk_nl(ii,j)   = mn.t_peak;   tpk_li(ii,j) = ml.t_peak;
        pk_nl(ii,j)    = mn.peak;     pk_li(ii,j)  = ml.peak;
        dl_all{ii,j}   = dl;          dn_all{ii,j} = dn;
    end
end

% ---- table (selected f only) ----
fprintf('\n%5s %5s | %10s %10s %7s | %6s %6s | %7s | %8s\n', ...
        'D','f','pk_nl(nM)','pk_lin(nM)','err%','tpk_nl','tpk_li','AUCerr%','nrmse%');
for ii = 1:nD
    for j = find(ismember(fList, fSub))
        fprintf('%5d %5.2f | %10.3f %10.3f %7.1f | %6.2f %6.2f | %7.1f | %8.1f\n', ...
            DList(doses(ii)), fList(j), pk_nl(ii,j)*1e3, pk_li(ii,j)*1e3, peak_err(ii,j), ...
            tpk_nl(ii,j), tpk_li(ii,j), auc_err(ii,j), nrmse(ii,j));
    end
end

% ---- Fig A: waveforms, one figure per dose (not for 1800 mg), six f each ----
jSub = find(ismember(fList, fSub));
for ii = 1:nD
    D = DList(doses(ii));
    if D > 400, continue, end
    fig = figure('Name', sprintf('LTI vs nonlinear, D = %d mg', D));
    tiledlayout(2, 3, 'TileSpacing','compact');
    for q = 1:numel(jSub)
        j = jSub(q);
        nexttile;
        plot(tl, dn_all{ii,j}*1e3, 'b-'); hold on
        plot(tl, dl_all{ii,j}*1e3, 'r--');
        xlabel('t (h)');  ylabel('\delta eda (nM)');
        title(sprintf('f = %.2f', fList(j)));  grid on
        if q == 1, legend('nonlinear','LTI','Location','best'); end
    end
    sgtitle(sprintf('Dose %d mg', D));
    savefig_pdf(fig, sprintf('fig7_lti_waveforms_D%d', D));
end

% ---- Fig B: peak error vs f, one curve per dose (<= 400 mg) ----
fig = figure('Name','Peak error of the LTI approximation');  hold on
for ii = 1:nD
    D = DList(doses(ii));
    if D > 400, continue, end
    plot(fList, peak_err(ii,:), '.-', 'MarkerSize', 12, 'DisplayName', sprintf('%d mg', D));
end
yline(0, 'k:', 'HandleVisibility', 'off');
set(gca, 'XDir','reverse', 'XScale','log');
xlabel('f (surviving fraction)');  ylabel('peak error LTI vs nonlinear (%)');
legend('Location','best');  grid on
savefig_pdf(fig, 'fig8_lti_peak_error');

% ---- Fig C: waveform error vs dose (tangent test) ----
fig = figure('Name','LTI error vs dose (tangent test)');  hold on
for j = jSub
    loglog(DList(doses), nrmse(:,j), 'o-', 'Color', col(j,:), ...
           'DisplayName', sprintf('f = %.2f', fList(j)));
end
yline(10, 'k--', '10 %', 'HandleVisibility', 'off');
set(gca, 'XScale','log', 'YScale','log');
xlabel('dose D (mg)');  ylabel('waveform error nrmse (%)');
legend('Location','northwest');  grid on
title('The error tends to 0 as D tends to 0 (1800 mg: non-clinical)');
savefig_pdf(fig, 'fig9_lti_error_vs_dose');

save(fullfile(root,'lti_comparison.mat'), 'peak_err','auc_err','nrmse','tpk_nl','tpk_li', ...
     'pk_nl','pk_li','tl','DList','fList','doses');
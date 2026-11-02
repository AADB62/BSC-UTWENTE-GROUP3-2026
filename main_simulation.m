% main_simulation.m
% Runs the nonlinear model for every (D, f) pair, computes the metrics and stores everything.
clear; clc;
root = fileparts(mfilename('fullpath'));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), ...
        fullfile(root,'analysis'),   fullfile(root,'test'));

p     = params();
DList = [0 10 100 250 400 1800];                       % mg: 0 control | 10 linear regime | 100-400 clinical | 1800 Reed-equivalent (non-clinical)
fList = round([1:-0.05:0.05, 0.02, 0.01], 2);          % surviving fraction of SNc cells (round avoids 0.2000000001-type values)
fSub  = [1 0.5 0.2 0.1 0.05 0.01];                     % f shown in the printed tables
tspan = [0 24];                                        % h

% Preallocate: res(i,j) <-> DList(i), fList(j)
res = repmat(struct('D',[],'f',[],'t',[],'x',[],'xe',[],'m',[]), ...
             numel(DList), numel(fList));

% ---- Simulations ----
dev0 = 0;                                              % largest drift with D = 0 (numerical check)
for i = 1:numel(DList)
    for j = 1:numel(fList)
        D = DList(i);
        f = fList(j);

        xe = equilibrium(p, f);                        % equilibrium of this scenario
        [t, x] = simulate(p, f, D, tspan);

        res(i,j).D  = D;
        res(i,j).f  = f;
        res(i,j).t  = t;
        res(i,j).x  = x;
        res(i,j).xe = xe;

        if D == 0
            dev0 = max(dev0, max(abs(x(:,9) - xe(9))));
        end
    end
    fprintf('dose %4d mg done (%d values of f)\n', DList(i), numel(fList));
end
fprintf('D = 0: largest deviation of eda from its equilibrium = %.2e uM\n', dev0);

% ---- Metrics ----
xe1   = equilibrium(p, 1);
xe02  = equilibrium(p, 0.2);
theta = 0.5*(xe1(9) + xe02(9));                        % Reed rule (decision D1)
fprintf('\ntheta = %.3f nM\n', theta*1e3);
fprintf('%5s %6s | %9s %7s %10s %9s %8s %8s %8s\n', ...
        'D','f','peak(nM)','t_pk(h)','AUC(nM*h)','t_abv(h)','t_half','rel','base>th');
for i = 1:numel(DList)
    for j = 1:numel(fList)
        if res(i,j).D == 0, continue, end              % no dose: nothing to measure
        m = metrics_eda(res(i,j).t, res(i,j).x, res(i,j).xe, theta);
        res(i,j).m = m;                                % keep them with the run (all f)
        if ismember(fList(j), fSub)                    % print only the reference f
            fprintf('%5d %6.2f | %9.3f %7.2f %10.3f %9.2f %8.2f %8.2f %8d\n', ...
                res(i,j).D, res(i,j).f, m.peak*1e3, m.t_peak, m.auc*1e3, ...
                m.t_above, m.t_half, m.rel_rise, m.baseline_above);
        end
    end
end

% ---- Sensitivity of theta to the reference f (decision D1) ----
% Same rule as before, but changing which diseased state defines the threshold.
fRefList = [0.2 0.1 0.05];
thetaK   = zeros(size(fRefList));
for k = 1:numel(fRefList)
    xeRef = equilibrium(p, fRefList(k));
    thetaK(k) = 0.5*(xe1(9) + xeRef(9));
end

jPrint = find(ismember(fList, fSub));
for i = find(DList > 0 & DList <= 400)                 % clinical doses only
    fprintf('\nD = %d mg: time above theta (h); * = baseline already above theta\n', DList(i));
    fprintf('%6s', 'f');
    for k = 1:numel(fRefList)
        fprintf(' | ref f=%-4.2f (%.3f nM)', fRefList(k), thetaK(k)*1e3);
    end
    fprintf('\n');
    for j = jPrint
        fprintf('%6.2f', fList(j));
        for k = 1:numel(fRefList)
            m = metrics_eda(res(i,j).t, res(i,j).x, res(i,j).xe, thetaK(k));
            mark = ' '; if m.baseline_above, mark = '*'; end
            fprintf(' | %19.2f%s           ', m.t_above, mark);
        end
        fprintf('\n');
    end
end

save(fullfile(root,'results.mat'), 'res', 'DList', 'fList', 'tspan', 'theta', 'thetaK', 'fRefList');
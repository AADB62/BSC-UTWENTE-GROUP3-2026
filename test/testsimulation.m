clear;clc;

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), ...
        fullfile(root,'analysis'),   fullfile(root,'test'));

makePlots = true;                  % set to false to run only the tests
p     = params();
scale = [100; 5; 5; 3; 80; 5; 3; 10; 0.002];
flist = [1 0.5 0.2 0.1 0.05 0.01];

%% ---------- TEST 3: with D = 0 the system stays at the equilibrium ----------
for f = flist
    xe = equilibrium(p, f);
    [t, x] = simulate(p, f, 0, [0 24]);
    drift = max(abs(x  - xe.') ./ scale.', [], 'all');
    assert(drift < 1e-6, 'drifts away from equilibrium at f=%g: %.1e', f, drift)
end
fprintf('TEST 3 passed: D=0 stays at the equilibrium (all f)\n')

%% ---------- TEST 4: no state becomes negative (with a dose) ----------
for D = [100 250]
    for f = flist
        [t, x] = simulate(p, f, D, [0 24]);
        worst = min(x, [], 'all');
        assert(worst > -1e-7, 'negative state at D=%g f=%g: %.1e', D, f, worst)
    end
end
fprintf('TEST 4 passed: positivity (D = 100, 250 mg; all f)\n')

%% ---------- Quick look: 100 mg, f = 1 (not pass/fail) ----------
[t, x] = simulate(p, 1, 100, [0 24]);
xe     = equilibrium(p, 1);
[pk, i] = max(x(:,9) - xe(9));
fprintf('100 mg, f=1: %d output points, peak of delta eda = %.3e uM at t = %.2f h\n', numel(t), pk, t(i));

%% ---------- Comparison A: tolerance (convergence study, seed) ----------
% Same case solved with different RelTol; the finest one is the reference.
% Output on a fixed grid tg so the curves can be subtracted.
f = 0.1;  D = 100;
tg   = linspace(0, 24, 2000).';
xeA  = equilibrium(p, f);
tols = {[], 1e-6, 1e-8, 1e-10};                    % [] = MATLAB defaults
names = {'default', '1e-6', '1e-8', '1e-10'};
DEDA = zeros(numel(tg), numel(tols));  NST = zeros(1, numel(tols));
for k = 1:numel(tols)
    [xg, NST(k)] = run15s(p, f, D, tg, tols{k}, scale);
    DEDA(:,k) = xg(:,9) - xeA(9);
end
ref = DEDA(:,end);
fprintf('\nTolerance study (f=%g, D=%g mg), reference = RelTol 1e-10\n', f, D)
fprintf('%-8s %8s %14s %10s %14s %12s\n', 'RelTol','steps','peak (uM)','t_peak(h)','AUC (uM*h)','dPeak (rel)');
for k = 1:numel(tols)
    [pk, i] = max(DEDA(:,k));
    auc = trapz(tg, DEDA(:,k));
    [pkr, ~] = max(ref);
    fprintf('%-8s %8d %14.6e %10.3f %14.6e %12.1e\n', names{k}, NST(k), pk, tg(i), auc, abs(pk-pkr)/pkr);
end

%% ---------- Comparison B: solvers (stiffness), f = 0.01, first hour ----------
fB = 0.01;  DB = 100;  tf = 10;
x0B = equilibrium(p, fB);  x0B(1) = x0B(1) + DB;
hB  = @(t,x) rhs(t,x,p,fB);
oB  = odeset('RelTol',1e-6,'AbsTol',1e-6*scale);
solvers = {@ode15s, @ode23s, @ode45, @ode23};
snames  = {'ode15s','ode23s','ode45','ode23'};
nst = zeros(1,4);  Tm = zeros(1,4);  treach = zeros(1,4);
fprintf('\nSolver comparison (f=%g, D=%g mg, 0-%g h, RelTol 1e-6; explicit ones stopped after 20 s)\n', fB, DB, tf)
for k = 1:4
    [nst(k), Tm(k), treach(k)] = runsolver(solvers{k}, hB, x0B, tf, oB);
    fprintf('%-7s %8d steps  %7.2f s  reached t = %.3f h\n', snames{k}, nst(k), Tm(k), treach(k));
end

%% ---------- Figures ----------
if makePlots
    tp = linspace(0, 24, 2000).';

    % Fig 1: delta eda(t) for every f (D = 100 mg)
    figure('Name','delta eda vs f'); hold on
    for f = flist
        xe_ = equilibrium(p, f);
        [xg, ~] = run15s(p, f, 100, tp, 1e-8, scale);
        plot(tp, xg(:,9) - xe_(9), 'LineWidth', 1.4, 'DisplayName', sprintf('f = %g', f));
    end
    xlabel('time (h)'); ylabel('\delta eda (\muM)'); grid on; legend('Location','northeast')
    title('Extracellular dopamine response to 100 mg, by surviving fraction f')

    % Fig 2: all states, f = 1, D = 100 mg
    xe1 = equilibrium(p, 1);
    [xg, ~] = run15s(p, 1, 100, tp, 1e-8, scale);
    figure('Name','states f=1 D=100');
    subplot(2,2,1); plot(tp, xg(:,1:3), 'LineWidth', 1.2); grid on
    legend('A_g (mg)','C_p (\muM)','C_{per} (\muM)'); xlabel('time (h)'); title('Body')
    subplot(2,2,2); plot(tp, xg(:,6:8), 'LineWidth', 1.2); grid on
    legend('ld5','cda5','vda5'); xlabel('time (h)'); ylabel('\muM'); title('5-HT terminal')
    subplot(2,2,3); plot(tp, 100*(xg(:,[4 5 9]) ./ xe1([4 5 9]).' - 1), 'LineWidth', 1.2); grid on
    legend('cda','vda','eda'); xlabel('time (h)'); ylabel('change vs equilibrium (%)'); title('DA terminal')
    subplot(2,2,4); plot(tp, xg(:,9) - xe1(9), 'LineWidth', 1.4); grid on
    xlabel('time (h)'); ylabel('\delta eda (\muM)'); title('Output')

    % Fig 3: tolerance study
    figure('Name','tolerance study');
    subplot(1,2,1); hold on
    plot(tg, DEDA(:,1), 'LineWidth', 1.2, 'DisplayName', 'default');
    plot(tg, ref, '--', 'LineWidth', 1.2, 'DisplayName', 'RelTol 1e-10');
    xlabel('time (h)'); ylabel('\delta eda (\muM)'); grid on; legend; title('Default vs reference')
    subplot(1,2,2); hold on
    for k = 1:numel(tols)-1
        semilogy(tg, abs(DEDA(:,k) - ref) + eps, 'LineWidth', 1.2, 'DisplayName', names{k});
    end
    set(gca,'YScale','log'); xlabel('time (h)'); ylabel('|\delta eda - reference| (\muM)'); grid on; legend
    title('Difference with the reference')

    % Fig 4: solver comparison
    figure('Name','solvers');
    bar(categorical(snames, snames), nst); set(gca,'YScale','log'); grid on
    ylabel('steps (0-1 h)'); title(sprintf('Steps per solver, f=%g (explicit ones may have been stopped)', fB))
end

%% ---------- local functions (must be at the end of the file) ----------
function [xg, nsteps] = run15s(p, f, D, tg, reltol, scale)
% ode15s with output on the grid tg. reltol = [] uses MATLAB defaults.
    xe = equilibrium(p, f);  x0 = xe;  x0(1) = x0(1) + D;
    h  = @(t,x) rhs(t,x,p,f);
    if isempty(reltol)
        opts = odeset();
    else
        opts = odeset('RelTol', reltol, 'AbsTol', reltol*scale);
    end
    sol    = ode15s(h, [tg(1) tg(end)], x0, opts);
    xg     = deval(sol, tg).';          % rows = times, columns = states
    nsteps = numel(sol.x) - 1;
end

function [nsteps, T, treached] = runsolver(solver, h, x0, tf, opts)
% Runs one solver and stops it after 20 s so ode45/ode23 cannot hang the test.
    opts = odeset(opts, 'OutputFcn', @timeguard);
    tic
    sol      = solver(h, [0 tf], x0, opts);
    T        = toc;
    nsteps   = numel(sol.x) - 1;
    treached = sol.x(end);
end

function status = timeguard(~, ~, flag)
    persistent t0
    status = false;
    switch flag
        case 'init'
            t0 = tic;
        case ''
            if toc(t0) > 20, status = true; end
    end
end

%%
% test/plot_eda_vs_f.m
% Extracellular DA (absolute, nM) vs time for several f.
% Top: D = 0 (baselines). Bottom: dose D (mg), comparable to Reed Fig. 7A.

root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), fullfile(root,'analysis'));

p      = params();
fList  = [1 0.2 0.1 0.05 0.01];
D_list = [0 100];                 % mg; change the second value to try other doses
tEnd   = 16;                      % h after the dose
thr    = 1.65;                    % therapeutic threshold in Reed Fig. 7A (nM, by eye)

scale = [100;5;5;3;80;5;3;10;0.002];
opts  = odeset('RelTol',1e-8,'AbsTol',1e-8*scale);
cols  = lines(numel(fList));

figure('Name','eda vs f','Position',[100 100 700 650]);
for d = 1:numel(D_list)
    subplot(2,1,d); hold on; grid on;
    for k = 1:numel(fList)
        f  = fList(k);
        xe = equilibrium(p, f);
        x0 = xe; x0(1) = x0(1) + D_list(d);       % dose enters Ag
        [t, X] = ode15s(@(t,x) rhs(t,x,p,f), [0 tEnd], x0, opts);
        plot(t, X(:,9)*1e3, 'Color', cols(k,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('f = %g', f));
    end
    yline(thr, 'k--', 'DisplayName', 'threshold (Reed, by eye)');
    xlabel('time after dose (h)'); ylabel('extracellular DA (nM)');
    title(sprintf('D = %g mg', D_list(d)));
    legend('Location','northeast');
end
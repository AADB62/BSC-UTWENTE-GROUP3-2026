% saturation_test.m  -  which step makes the dose-response sublinear?
% Counterfactual: make one candidate step 10x less saturable and recompute the linearity index
%   L = [peak(250 mg)/peak(100 mg)] / (250/100)   (L = 1: linear response)
clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), fullfile(root,'analysis'));
p0 = params();

cand   = {'VmM','KM','VmL','Kap','VmD','KD','VmM5','KM5','kMAO'};   % kMAO: control (linear step)
factor = 10;
fTest  = [1 0.01];

for f = fTest
    fprintf('\nf = %.2f\n', f);
    fprintf('%-12s %8s\n', 'case', 'L');
    fprintf('%-12s %8.3f\n', 'nominal', linearity(p0, f));
    for k = 1:numel(cand)
        p = p0;
        p.(cand{k}) = p0.(cand{k})*factor;        % Vmax x10 (more capacity) or K x10 (less saturation)
        fprintf('%-8s x%-3d %8.3f\n', cand{k}, factor, linearity(p, f));
    end
end

function L = linearity(p, f)
    L = NaN;
    try
        doses = [100 250];
        xe = equilibrium(p, f);
        pk = zeros(1,2);
        for i = 1:2
            [t, x] = simulate(p, f, doses(i), (0:0.02:24)');
            pk(i) = max(x(:,9) - xe(9));
        end
        L = (pk(2)/pk(1)) / (doses(2)/doses(1));
    catch
    end
end
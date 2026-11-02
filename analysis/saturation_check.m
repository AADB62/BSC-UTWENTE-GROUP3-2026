% saturation_check.m  -  which step saturates the response to the dose?
% For each dose: occupancy of each Michaelis-Menten step and fate of the entering L-Dopa.
clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), fullfile(root,'analysis'));
load(fullfile(root,'results.mat'), 'res', 'DList', 'fList');
p = params();

for f = [1 0.01]
    j = find(fList == f, 1);
    fprintf('\nf = %.2f\n', f);
    fprintf('%5s | %8s %8s %8s %8s | %10s %10s %8s\n', 'D', 'max Cp', 'VL occ', 'VMAT occ', 'DAT occ', ...
            'int VL/D', 'int dMAO/D', 'MAO/VL');
    for i = find(DList > 0)
        t = res(i,j).t;  x = res(i,j).x;  xe = res(i,j).xe;  D = DList(i);
        Cp = x(:,2);  cda = x(:,4);  eda = x(:,9);

        occVL   = max(Cp ./ (p.Kap + Cp));          % fraction of Vmax of L-Dopa uptake
        occVMAT = max(cda ./ (p.KM + cda));         % fraction of Vmax of VMAT2
        occDAT  = max(eda ./ (p.KD + eda));         % fraction of Vmax of DAT

        VL   = p.VmL*Cp ./ (p.Kap + Cp);            % uM/h, L-Dopa entering the terminal (0 at D = 0)
        dMAO = p.kMAO*(cda - xe(4));                % uM/h, extra MAO flux caused by the dose
        intVL   = trapz(t, VL);                     % uM
        intdMAO = trapz(t, dMAO);                   % uM

        fprintf('%5d | %8.2f %8.3f %8.3f %8.4f | %10.3f %10.3f %8.2f\n', ...
            D, max(Cp), occVL, occVMAT, occDAT, intVL/D, intdMAO/D, intdMAO/intVL);
    end
end
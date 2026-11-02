clear;  clc;
% TEST 2c: equilibrium() vs the V-3 reference (PART 6: fsolve with exact Jacobian;
% PART 7: independent fzero on one scalar equation in eda).
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), ...
        fullfile(root,'analysis'),   fullfile(root,'test'));

[Fnum, Jnum] = build_Fnum_V3();
p     = params();
flist = [1 0.5 0.2 0.1 0.05 0.01];
seed  = [0;0;0; 2.65; 79.1; 0;0;0; 0.00198];

% --- reference 1: PART 6 of the V-3 script ---
optsV3 = optimoptions('fsolve','Display','off','SpecifyObjectiveGradient',true, ...
         'FunctionTolerance',1e-12,'StepTolerance',1e-14,'OptimalityTolerance',1e-12);

% --- reference 2: PART 7 (typed by hand, eds = 0.0019746) ---
x5f = @(z,f) 8000*z./(0.2+z) + 30*z./(3+z) + (400./f).*z;
x4f = @(x5) 3*(1.4*x5/236)./(1 - 1.4*x5/236);
gD2 = @(z) 4.5./(8*(z/0.0019746).^4+1) + 0.5;
VTHs = @(z,x4) 0.56/(1+126/160) * gD2(z) * 125*126*319 ./ ...
               (126*319 + 46*319 + 46*60*(1+x4/110));
Phi  = @(z,f) VTHs(z, x4f(x5f(z,f))) - 10*x4f(x5f(z,f)) - 30*z./(3+z) - (400./f).*z;
x5max = 0.99*236/1.4;

xv = seed;  xm = seed;
for k = 1:numel(flist)
    fk = flist(k);

    FJ = @(z) deal(Fnum(z,fk), Jnum(z,fk));
    [xv, ~, flagV] = fsolve(FJ, xv, optsV3);          % V-3 PART 6
    [xm, flagM]    = equilibrium(p, fk, xm);          % ours

    zmax = fzero(@(z) x5f(z,fk) - x5max, [0 0.01]);   % V-3 PART 7
    zf   = fzero(@(zz) Phi(zz,fk), [1e-7 zmax]);

    scale = [1;1;1; xv(4); xv(5); 1;1;1; xv(9)];      % body and 5-HT are 0
    err6  = max(abs(xm - xv) ./ scale);
    err7  = abs(xm(9) - zf);

    fprintf('f=%5.2f  ours: cda=%7.4f vda=%8.3f eda=%.6f | V-3: cda=%7.4f vda=%8.3f eda=%.6f | err6=%.1e err7(eda)=%.1e\n', ...
            fk, xm(4), xm(5), xm(9), xv(4), xv(5), xv(9), err6, err7);
    fprintf('   flagV=%d (V-3 fsolve)  flagM=%d (ours)  |F|max V-3=%.1e  ours=%.1e\n', ...
            flagV, flagM, max(abs(Fnum(xv,fk))), max(abs(rhs(0,xm,p,fk))));
    assert(flagM > 0, 'our equilibrium() did not converge at f=%g', fk)
    assert(flagV > 0, 'V-3 reference fsolve did not converge at f=%g', fk)
    assert(err6 < 1e-8, 'differs from V-3 PART 6 at f=%g', fk)
    assert(err7 < 1e-9, 'differs from V-3 PART 7 (fzero) at f=%g', fk)
end
fprintf('TEST 2c passed: equilibrium() matches V-3 (fsolve and fzero) for all f\n')
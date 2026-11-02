% Add the project folders to the path (this file lives in <root>/test)
root = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(root,'parameters'), fullfile(root,'model'), ...
        fullfile(root,'analysis'),   fullfile(root,'test'));


% TEST 1: rhs vs Fnum (V-3) at random positive points, several f.
% Self-contained: Fnum is rebuilt below from PART 1-2 of MATLABscriptbcs_V3.m
% (Symbolic Math Toolbox needed). Fnum is only a reference for this test.

Fnum  = build_Fnum_V3();
p     = params();
flist = [1 0.5 0.2 0.1 0.05 0.01];
sc    = [100; 5; 5; 3; 80; 5; 3; 10; 0.01];   % typical size of each state
rng(1)                                         % reproducible random numbers

worst = 0;
for k = 1:numel(flist)
    fk = flist(k);
    for n = 1:10
        x   = sc .* (0.1 + rand(9,1));         % positive, never 0
        a   = rhs(0, x, p, fk);
        b   = Fnum(x, fk);
        err = max(abs(a - b) ./ abs(b));       % relative error per component
        worst = max(worst, err);
        assert(err < 1e-10, 'FAILED at f=%g, point %d: err=%.2e', fk, n, err)
    end
end
fprintf('TEST 1 passed: 60 points, worst relative error = %.2e\n', worst);



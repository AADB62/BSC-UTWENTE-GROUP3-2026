function [xe, flag] = equilibrium(p, f, x0,verbose)
% Equilibrium of the 9-state model (no dose) for a given f.
% Optional x0: starting guess (default: Reed healthy state).
    if nargin < 4
        verbose=false;
    end
    if nargin < 3 || isempty(x0)
        x0 = [0;0;0; 2.65; 79.1; 0;0;0; 0.00198];
    end
    opts = optimoptions('fsolve', 'Display','off', ...
            'FunctionTolerance',1e-12, 'StepTolerance',1e-14, ...
            'OptimalityTolerance',1e-8, ...
            'TypicalX',[1;1;1; 2.6; 79; 1;1;1; 0.002]);
    if verbose
        opts = optimoptions(opts, 'Display','iter', ...
                'PlotFcn',@optimplotfirstorderopt);
    end
    h = @(x) rhs(0, x, p, f);
    [xe, ~, flag] = fsolve(h, x0, opts);
end
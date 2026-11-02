function [t, x] = simulate(p, f, D, tspan)
    xe = equilibrium(p, f); % equilibrium for this f
    x0 = xe;
    x0(1) = xe(1) + D;      % dose in the gut
    scale  = [100; 5; 5; 3; 80; 5; 3; 10; 0.002];
    RelTol = 1e-8;%Use  the defoult  not to argue
    opts   = odeset('RelTol', RelTol, 'AbsTol', RelTol*scale);
    [t, x] = ode15s(@(t,x) rhs(t,x,p,f), tspan, x0, opts);
    %compare with ode45 ode23
end
function [A, B, C, Dss, xe] = lineal(p, f)
% LINEAL  Local LTI approximation of the nonlinear model around x_e(f).
%   d(dx)/dt = A dx + B u ,  d(eda) = C dx + Dss u
%   u = dose impulse (mg): x(0+) = x_e + D*e1  (dose enters as Ag(0) = D)
%
%   A   = dF/dx at x_e(f)      (9x9, complex-step derivative)
%   B   = e1                   (the dose enters only the gut state Ag)
%   C   = e9'                  (output = extracellular dopamine eda)
%   Dss = 0                    (no direct path from the impulse to eda)

    xe = equilibrium(p, f);          % operating point: rhs(xe) = 0
    n  = numel(xe);                  % 9 states
    h  = 1e-30;                      % complex step: no subtractive cancellation

    A = zeros(n, n);
    for k = 1:n
        xc    = xe;
        xc(k) = xc(k) + 1i*h;        % perturb state k along the imaginary axis
        A(:,k) = imag( rhs(0, xc, p, f) ) / h;   % column k = dF/dx_k
    end

    B   = zeros(n, 1);  B(1) = 1;    % dose enters Ag
    C   = zeros(1, n);  C(9) = 1;    % output is eda
    Dss = 0;                         % no direct feedthrough
end
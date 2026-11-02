function out = eigen_table(p, fList)
% EIGEN_TABLE  Eigenvalues and time constants of A(f) for each f.
%   out.lam(:,j) : eigenvalues (h^-1), sorted from slowest to fastest
%   out.tau(:,j) : time constants tau = -1/Re(lambda), in seconds

    n   = 9;
    lam = zeros(n, numel(fList));
    tau = zeros(n, numel(fList));

    for j = 1:numel(fList)
        A = lineal(p, fList(j));
        l = eig(A);
        [~, idx] = sort(abs(real(l)));     % slowest first
        l = l(idx);
        lam(:,j) = l;
        tau(:,j) = -3600 ./ real(l);       % s  (lambda in h^-1)
    end

    % ---- Print: eigenvalues (h^-1) ----
    fprintf('\nEigenvalues of A(f) (h^-1)\n');
    fprintf('%4s', 'mode');
    fprintf(' | f=%-9.2f', fList);  fprintf('\n');
    for i = 1:n
        fprintf('%4d', i);
        for j = 1:numel(fList), fprintf(' | %11.4g', real(lam(i,j))); end
        fprintf('\n');
    end

    % ---- Print: time constants ----
    fprintf('\nTime constants tau = -1/Re(lambda), in minutes unless marked s\n');
    fprintf('%4s', 'mode');
    fprintf(' | f=%-9.2f', fList);  fprintf('\n');
    for i = 1:n
        fprintf('%4d', i);
        for j = 1:numel(fList)
            t = tau(i,j);
            if t < 60, fprintf(' | %8.3g s ', t);
            else,      fprintf(' | %8.3g min', t/60);
            end
        end
        fprintf('\n');
    end

    % ---- Checks ----
    fprintf('\nmax |Im(lambda)| = %.2e h^-1 (0 -> no oscillations)\n', max(abs(imag(lam(:)))));
    fprintf('max Re(lambda)   = %.3e h^-1 (<0 -> stable)\n', max(real(lam(:))));

    out.lam = lam;  out.tau = tau;  out.fList = fList;
end
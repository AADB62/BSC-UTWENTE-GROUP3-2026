% print_prcc_all.m  -  PRCC of ALL outputs from prcc.mat (top 5 per output and f)
clear; clc;
root = fileparts(fileparts(mfilename('fullpath')));
load(fullfile(root,'prcc.mat'));          % R, SIG, names, fList, outNames
for m = 1:numel(outNames)
    for jf = 1:numel(fList)
        [~,ord] = sort(abs(R(:,m,jf)),'descend');
        fprintf('\n%s, f = %.2f\n', outNames{m}, fList(jf));
        for q = ord(1:5)'
            fprintf('  %-7s %7.3f %s\n', names{q}, R(q,m,jf), repmat('*',1,SIG(q,m,jf)));
        end
    end
end
% parameters with negligible PRCC for the peak at every f (max |PRCC| < 0.1)
mx = max(abs(squeeze(R(:,2,:))),[],2);
fprintf('\nmax |PRCC| over f for the peak (ascending):\n');
[v,i] = sort(mx); for q = 1:10, fprintf('  %-7s %6.3f\n', names{i(q)}, v(q)); end
fprintf('\nka, Vp, Vc, CL: PRCC of peak and t_peak by f\n');
for nm = {'ka','Vp','Vc','CL'}
    k = strcmp(names, nm{1});
    fprintf('%-4s peak  :', nm{1}); fprintf(' %6.2f', squeeze(R(k,2,:))); fprintf('\n');
    fprintf('%-4s t_peak:', nm{1}); fprintf(' %6.2f', squeeze(R(k,3,:))); fprintf('\n');
end
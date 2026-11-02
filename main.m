% run_all.m  -  regenerates every result and figure in order
clear; clc;
root = fileparts(mfilename('fullpath'));
run(fullfile(root,'main_simulation.m'));        % simulations and metrics -> results.mat
run(fullfile(root,'analysis','make_figures.m'));
run(fullfile(root,'analysis','compareLTI.m'));
run(fullfile(root,'analysis','sensitivity.m'));
run(fullfile(root,'analysis','run_prcc.m'));
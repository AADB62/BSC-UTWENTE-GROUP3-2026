clear;clc;

root = pwd;   % project root
addpath(fullfile(root,'parameters'), fullfile(root,'model'), fullfile(root,'analysis'));
p   = params();
out = eigen_table(p, [1 0.5 0.2 0.1 0.05 0.01]);
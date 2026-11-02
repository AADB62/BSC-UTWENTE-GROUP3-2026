clear;clc;
p  = params();
serotonin_fluxes(0,0,0,0,p)

% Test: DA fluxes at the healthy equilibrium (f = 1, no dose)
p  = params();
f  = 1;
Cp = 0;  cda = 2.647;  vda = 79.02;  eda = 0.001975;

fl = dopamine_fluxes(Cp, eda, cda, vda, f, p)

% Balance of cda (VL = 0 because Cp = 0): should be ~0
res_cda = fl.VL + fl.VTH - fl.VMAT + fl.VDAT - fl.MAO
% Balance of vda: should be ~0
res_vda = fl.VMAT - fl.Rel
% Balance of eda (5-HT contribution is 0 without dose): should be ~0
res_eda = fl.Rel - fl.VDAT - fl.VCOM - fl.Kclr
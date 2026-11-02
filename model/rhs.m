function dxdt = rhs(t, x, p, f)
% Right-hand side of the 9-state model. x = [Ag;Cp;Cper;cda;vda;ld5;cda5;vda5;eda]
    Ag = x(1); Cp = x(2); Cper = x(3);
    cda = x(4); vda = x(5); eda = x(9);
    ld5 = x(6); cda5 = x(7); vda5 = x(8);

    fb = body_fluxes(Ag,Cp,Cper,p);
    fs = serotonin_fluxes(Cp,ld5,cda5,vda5,p);
    fd = dopamine_fluxes(Cp,eda,cda,vda,f,p);

    % derivatives as sums of fluxes
    dxdt = [ -fb.Jabs;    % 1  Ag
             (p.c*fb.Jabs - fb.Jel - fb.Jex)/p.Vc;    % 2  Cp
             fb.Jex/p.Vp;    % 3  Cper
             fd.VL + fd.VTH - fd.VMAT + fd.VDAT - fd.MAO;    % 4  cda
             fd.VMAT - fd.Rel;    % 5  vda
             fs.VL5 - fs.VA5;    % 6  ld5
             fs.VA5 - fs.VM5 - fs.MAO5;    % 7  cda5
             fs.VM5 - fs.R5;    % 8  vda5
             fd.Rel + (p.f5/f)*fs.R5 - fd.VDAT - fd.VCOM - fd.Kclr ];  % 9  eda
end
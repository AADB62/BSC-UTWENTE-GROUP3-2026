function fl=serotonin_fluxes(Cp,ld5,cda5,vda5,p)
    fl.VL5=mm(p.VmL5, p.Kap, Cp);
    fl.VA5=mm(p.VmA5, p.KA5, ld5);
    fl.VM5=mm(p.VmM5, p.KM5, cda5) - p.kout5*vda5;
    fl.MAO5=p.kMAO5*cda5;
    fl.R5=p.RF*vda5;
end
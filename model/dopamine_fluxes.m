function fl=dopamine_fluxes(Cp,eda,cda,vda,f,p)
    fl.VL= 	mm(p.VmL, p.Kap, Cp);
    g     = g_auto(eda, p.eds);
    num   = p.VmTH * p.tyr * p.bh4;
    den   = p.tyr*p.bh4 + p.Ktyr*p.bh4 + p.Ktyr*p.Kbh4*(1 + cda/p.Kicda);
    fl.VTH = 0.56/(1 + p.tyr/p.Kityr) * g * num/den;
    fl.VMAT= mm(p.VmM, p.KM, cda) - p.kout*vda;
    fl.VDAT = mm(p.VmD, p.KD, eda);
    fl.MAO = 	p.kMAO*cda;
    fl.Rel= p.fire*vda;
    fl.VCOM=mm(p.VmC, p.KC, eda);
    fl.Kclr = 	(p.k0/f)*eda;
end
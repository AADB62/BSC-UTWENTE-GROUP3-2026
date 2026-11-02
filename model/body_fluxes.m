function fl=body_fluxes(Ag,Cp,Cper,p)
    fl.Jabs=p.ka*Ag;
    fl.Jel=p.CL*Cp;
    fl.Jex=p.Q*(Cp-Cper);
end
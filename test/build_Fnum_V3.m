function [Fnum, Jnum] = build_Fnum_V3()
% Reference model from MATLABscriptbcs_V3.m (PART 1-2, eds = 0.0019746).
% Returns numeric F(x,f) and its exact Jacobian J(x,f). Used only by tests.
    syms x1 x2 x3 x4 x5 x6 x7 x8 x9 real
    x = [x1;x2;x3;x4;x5;x6;x7;x8;x9];

    ka   = sym('ka','positive');    c    = sym('c','positive');
    CL   = sym('CL','positive');    Q    = sym('Q','positive');
    Vc   = sym('Vc','positive');    Vp   = sym('Vp','positive');
    VmL  = sym('VmL','positive');   VmL5 = sym('VmL5','positive');
    Kap  = sym('Kap','positive');
    eds  = sym('eds','positive');   tyr  = sym('tyr','positive');
    bh4  = sym('bh4','positive');   Ktyr = sym('Ktyr','positive');
    Kbh4 = sym('Kbh4','positive');  Kicda = sym('Kicda','positive');
    Kityr = sym('Kityr','positive'); VmTH = sym('VmTH','positive');
    VmM  = sym('VmM','positive');   KM   = sym('KM','positive');
    kout = sym('kout','positive');  fire = sym('fire','positive');
    VmD  = sym('VmD','positive');   KD   = sym('KD','positive');
    VmC  = sym('VmC','positive');   KC   = sym('KC','positive');
    kMAO = sym('kMAO','positive');  k0   = sym('k0','positive');
    f5   = sym('f5','positive');    RF   = sym('RF','positive');
    f    = sym('f','positive');
    VmA5 = sym('VmA5','positive');  KA5  = sym('KA5','positive');
    VmM5 = sym('VmM5','positive');  KM5  = sym('KM5','positive');
    kout5 = sym('kout5','positive'); kMAO5 = sym('kMAO5','positive');

    VL   = VmL *x2/(Kap+x2);
    VL5  = VmL5*x2/(Kap+x2);
    g    = 4.5/(8*(x9/eds)^4+1) + 0.5;
    VTH  = 0.56/(1+tyr/Kityr) * g * VmTH*tyr*bh4 / ...
           (tyr*bh4 + Ktyr*bh4 + Ktyr*Kbh4*(1+x4/Kicda));
    VMAT = VmM*x4/(KM+x4) - kout*x5;
    VDAT = VmD*x9/(KD+x9);
    VCOM = VmC*x9/(KC+x9);
    VA5  = VmA5*x6/(KA5+x6);
    VM5  = VmM5*x7/(KM5+x7) - kout5*x8;

    F = [ -ka*x1;
          (c*ka*x1 - CL*x2 - Q*(x2-x3))/Vc;
          Q*(x2-x3)/Vp;
          VL + VTH - VMAT + VDAT - kMAO*x4;
          VMAT - fire*x5;
          VL5 - VA5;
          VA5 - VM5 - kMAO5*x7;
          VM5 - RF*x8;
          fire*x5 + (f5/f)*RF*x8 - VDAT - VCOM - (k0/f)*x9 ];
    J = jacobian(F, x);

    S    = [ka, c, CL, Q, Vc, Vp, VmL, VmL5, Kap, ...
            eds, tyr, bh4, Ktyr, Kbh4, Kicda, Kityr, VmTH, ...
            VmM, KM, kout, fire, ...
            VmD, KD, VmC, KC, kMAO, k0, f5, RF, ...
            VmA5, KA5, VmM5, KM5, kout5, kMAO5];
    vals = [2.4, 5.071, 24.8, 6.8, 58.5, 72.9, 1840, 1750, 238.4, ...
            0.0019746, 126, 319, 46, 60, 110, 160, 125, ...
            236, 3, 0.4, 1, ...
            8000, 0.2, 30, 3, 10, 400, 0.70, 1, ...
            400, 131.8, 40, 10.5, 0.4, 10];
    assert(numel(S) == numel(vals))
    Fnum = matlabFunction(subs(F, S, vals), 'Vars', {x, f});
    Jnum = matlabFunction(subs(J, S, vals), 'Vars', {x, f});
end
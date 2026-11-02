function p = params()
% PARAMS_SUB1  Parameters of subsystem 1 (oral L-Dopa dose -> eda), V-3 calibration.
% Time in hours. Concentrations in uM, amounts in mg, volumes in L.
% f (surviving DA terminals) is NOT here: it is the scenario, passed separately.

    % --- Body (absorption, plasma, peripheral) ---
    p.ka    = 2.4;        % h^-1      gastric emptying / absorption
    p.c     = 5.071;      % umol/mg   mg -> umol conversion
    p.CL    = 24.8;       % L/h       clearance
    p.Q     = 6.8;        % L/h       inter-compartmental flow
    p.Vc    = 58.5;       % L         central volume
    p.Vp    = 72.9;       % L         peripheral volume

    % --- L-Dopa transport into the terminals ---
    p.VmL   = 1840;       % uM/h      uptake, DA terminal
    p.VmL5  = 1750;       % uM/h      uptake, 5-HT terminal
    p.Kap   = 238.4;      % uM        DERIVED (apparent K with tyr/trp competition)

    % --- TH step and D2 autoreceptor ---
    p.eds   = 0.0019746;  % uM        eda_ss, V-3 calibrated (Best et al.: 0.002024)
    p.tyr   = 126;        % uM        constant (A18)
    p.bh4   = 319;        % uM        constant
    p.Ktyr  = 46;         % uM
    p.Kbh4  = 60;         % uM
    p.Kicda = 110;        % uM        product inhibition by cda
    p.Kityr = 160;        % uM        substrate inhibition by tyr
    p.VmTH  = 125;        % uM/h      TH Vmax

    % --- DA terminal: VMAT2, leak, release ---
    p.VmM   = 236;        % uM/h      VMAT2 Vmax
    p.KM    = 3;          % uM
    p.kout  = 0.4;        % h^-1      vesicular leak
    p.fire  = 1;          % h^-1      release rate

    % --- Extracellular DA: DAT, COMT, MAO, clearance, 5-HT contribution ---
    p.VmD   = 8000;       % uM/h      DAT Vmax
    p.KD    = 0.2;        % uM
    p.VmC   = 30;         % uM/h      COMT Vmax
    p.KC    = 3;          % uM
    p.kMAO  = 10;         % h^-1      MAO (linear)
    p.k0    = 400;        % h^-1      clearance, k_rem = k0/f
    p.f5    = 0.70;       % -         surviving 5-HT terminals (constant)
    p.RF    = 1;          % h^-1      RFIRE, release rate 5-HT terminal

    % --- 5-HT terminal ---
    p.VmA5  = 400;        % uM/h      AADC Vmax
    p.KA5   = 131.8;      % uM        DERIVED
    p.VmM5  = 40;         % uM/h      VMAT2 Vmax
    p.KM5   = 10.5;       % uM        DERIVED (Reed erratum corrected, A24)
    p.kout5 = 0.4;        % h^-1      vesicular leak
    p.kMAO5 = 10;         % h^-1      MAO (linear)
end
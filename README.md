# system_simulation — Subsystem 1: oral L-DOPA dose → extracellular dopamine (eda)

MATLAB implementation of the nonlinear model of subsystem 1 (Group 3, Biological
Control Systems, University of Twente). Input: oral L-DOPA dose D [mg]. Output:
extracellular dopamine in the striatum, `eda`, as a function of the fraction `f`
of surviving dopaminergic terminals.

The model equations are frozen (calibration V-3, `eda_ss = 0.0019746 µM`).
This code does not modify them; it simulates and analyses them.

## State vector and units

`x = [Ag; Cp; Cper; cda; vda; ld5; cda5; vda5; eda]`

| State | Meaning | Units |
|-------|---------|-------|
| Ag | L-DOPA amount in the gut | mg |
| Cp, Cper | L-DOPA in plasma / peripheral compartment | µM |
| cda, vda | cytosolic / vesicular dopamine, dopaminergic terminal | µM |
| ld5, cda5, vda5 | L-DOPA, cytosolic and vesicular dopamine, serotonergic terminal | µM |
| eda | extracellular dopamine | µM |

Time is in hours. A dose is applied as the initial condition `Ag(0) = D` on top of
the equilibrium of the same scenario. `f` is constant during a simulation;
`f5 = 0.70` is fixed in `params.m`.

## Folder structure

| Folder / file | Content |
|---------------|---------|
| `parameters/params.m` | `params()` returns the struct `p` with the 35 parameters (names as in the symbol table). `f` is not part of `p`. |
| `model/` | `rhs.m` (right-hand side), `body_fluxes.m`, `serotonin_fluxes.m`, `dopamine_fluxes.m`, `mm.m` (Michaelis–Menten), `g_auto.m` (autoreceptor feedback). Each flux is computed once; each derivative is a signed sum of fluxes. |
| `analysis/equilibrium.m` | Equilibrium `x_e(f)` with `fsolve` (accepted only if flag > 0, max abs residual < 1e-8 and states ≥ −1e-9). |
| `analysis/simulate.m` | Nonlinear simulation with `ode15s` from `x_e(f) + D·e1`. |
| `analysis/metrics_eda.m` | Metrics of δeda: peak, time of peak, AUC (0–24 h), time above θ, half-peak duration, peak/eda_eq. |
| `analysis/lineal.m`, `eigen_table.m` | Jacobian (complex step) at `x_e(f)`, eigenvalues and time constants per `f`. |
| `analysis/run_sensitivity.m`, `sensitivity_oat.m` | Local one-at-a-time normalised sensitivity (±10 %). |
| `analysis/run_prcc.m` | Global sensitivity: Latin hypercube sampling + PRCC. |
| `analysis/saturation_check.m`, `saturation_test.m` | Transporter occupancies and linearity index of the dose–response. |
| `analysis/make_figures.m` | Figures of the report → `figures/`. |
| `test/` | Verification scripts (see below). |
| `main_simulation.m` | Grid of scenarios, simulations, metrics → `results.mat`. |
| `main.m` | Runs the whole chain in order. |

## How to run

1. Open MATLAB with the project root as current folder (each script adds its own paths).
2. Run `main.m` to regenerate all results and figures, or run the steps separately:
   `main_simulation` → `analysis/make_figures` → `analysis/run_sensitivity` → `analysis/run_prcc`.
3. Results are saved to `results.mat`, `sensitivity.mat` and `prcc.mat`; figures to `figures/`.

Requirements: MATLAB with the Optimization Toolbox (`fsolve`). `ode15s` is part of
base MATLAB. The PRCC computation does not use the Statistics Toolbox.

## Scenarios

- Doses: D = 0, 10, 100, 250, 400, 1800 mg. D = 0 is a control (the system must stay at
  equilibrium); 100 mg is the reference of the body block; 1800 mg is a non-clinical
  dose used only as a diagnostic.
- Fraction of surviving terminals: f = 1 to 0.05 in steps of 0.05, plus 0.02 and 0.01
  (f = 0 cannot be simulated because `k0/f` and `f5/f` diverge).
- θ (reference for "time above threshold") = midpoint between the baseline `eda` of
  f = 1 and f = 0.2, following Reed et al. (2012) Sec. 3.3. It is a comparative
  criterion, not a clinical level.

## Tests (`test/`)

| Test | What it checks |
|------|----------------|
| `test_rhs_vs_Fnum.m` (with `build_Fnum_V3.m`) | `rhs` against the symbolic V-3 model at 60 random points (worst relative error 8.4e-14). |
| `testequilibrium.m`, `testequilibriumvsV3.m` | Equilibrium against the V-3 solution (`fsolve`, `fzero`). |
| `testsimulation.m` | D = 0 stays at equilibrium; no negative states. |
| `testlineal.m`, `testeigentable.m` | Complex-step Jacobian against finite differences; eigenvalues real and negative. |
| `test_bfluxes.m`, `test_dafluxes.m` | Individual flux blocks. |

Test files start with `test` so they never shadow a model function.

## Known limitations

- Constants 0.56 (in `VTH`) and 4.5, 8, 0.5 (in `g_auto`) are hard-coded as in the V-3 script.
- Local sensitivity covers 34 parameters (not `c`, a unit conversion) and not the hard-coded constants.
- PRCC is computed only at D = 100 mg with log-uniform ranges ×/÷2 (f5 limited to 0.35–1);
  these ranges are a modelling assumption.
- `f` and `fs` are constant; only single doses are simulated.

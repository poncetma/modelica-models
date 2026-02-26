model ReganLinearStirling
 "KRUSTY/ASC-E2 style linear Stirling model following Regan & Lewandowski structure.
   - 2 gas nodes (expansion/compression) + lumped regenerator volume (dead)
   - 4 mechanical states (displacer/piston positions & velocities)
   - single thermal integrator for hot expansion node
   - simple alternator (L-R) and electrical damping equivalent (B_elec)
   - per-cycle update of linear pressure-factor coefficients (simple heuristic)
   Note: replace 'pressure-factor' heuristics with your SDM linearization for exact match.
  "

  // 1) PARAMETERS (KRUSTY / ASC-E2 style values)
  parameter Real R = 2077.0 "J/(kg K) (helium)";
  parameter Real T_guess = 700.0 "K (reference temp for m_gas calc)";

  // geometry / volumes (mid-stroke references)
  parameter Real V_e0 = 4.0e-5 "expansion mid-stroke vol (m^3)";
  parameter Real V_c0 = 6.0e-5 "compression mid-stroke vol (m^3)";
  parameter Real V_reg = 2.0e-5 "regenerator void vol (m^3)";
  parameter Real V_dead = 1.0e-5 "other dead vol (m^3)";

  // effective areas (estimates for ASC-class)
  parameter Real A_p = 1.539e-4 "power piston effective area (m^2)";
  parameter Real A_d = 4.91e-4  "displacer effective face area (m^2)";

  // charge pressure (ASC/KRUSTY reported)
  parameter Real p_mean = 3.52e6 "Mean charge pressure (Pa)";

  // compute gas mass from mean conditions and nominal total volume
  parameter Real Vtot0 = V_e0 + V_c0 + V_reg + V_dead "nominal total dead volume (m^3)";
  parameter Real m_gas = p_mean * Vtot0 / (R * T_guess) "working gas mass (kg)";

  // mechanical masses and linear mechanical terms
  parameter Real m_p = 0.28 "power piston mass (kg)";
  parameter Real m_d = 0.18 "displacer mass (kg)";
  parameter Real k_p = 800.0 "piston centering stiffness (N/m)";
  parameter Real k_d = 200.0 "displacer centering stiffness (N/m)";
  parameter Real c_p = 15.0 "piston viscous damping (N s/m)";
  parameter Real c_d = 5.0  "displacer viscous damping (N s/m)";

  // alternator / electrical
  parameter Real K_e = 10.0 "EMF constant (V per m/s)";
  parameter Real L_alt = 1e-3 "Alternator inductance (H)";
  parameter Real R_alt = 0.5 "Alternator coil resistance (Ohm)";
  parameter Real R_load = 10.0 "External resistive load (Ohm)";
  // electrical damping equivalent (B = Ke^2 / R_load) - keep as parameter for quick tuning:
  parameter Real B_elec = K_e*K_e / (R_alt + R_load) "electrical damping (N s/m)";

  // thermal / regenerator lumped
  parameter Real UA_e = 200.0 "hot-side UA (W/K) (tunable)";
  parameter Real UA_c = 80.0  "cold-side UA (W/K) (tunable)";
  parameter Real UA_r = 300.0 "regenerator coupling UA (W/K) (tunable)";
  parameter Real C_r = 20.0   "regenerator thermal capacitance (J/K) (physics-based estimate)";

  // some solver / numeric parameters
  parameter Real eps_time = 1e-9;

  // 2) VARIABLES / STATES
  // mechanical states
  Real x_p(start=0.0) "power piston displacement (m) (+ increases compression volume)";
  Real v_p(start=0.0) "power piston velocity (m/s)";
  Real x_d(start=1e-3) "displacer displacement (m) (+ increases expansion volume)";
  Real v_d(start=0.0) "displacer velocity (m/s)";

  // thermal states
  Real T_hg(fixed=false) "hot/expansion gas node temperature (K)";
  Real T_cg(fixed=false) "cold/compression gas node temperature (K)";
  Real T_r(start=(700.0)) "regenerator temp (K)";

  // electrical state
  Real i_alt(start=0.0) "alternator current (A)";

  // derived algebraic quantities
  Real V_e "instantaneous expansion volume (m^3)";
  Real V_c "instantaneous compression volume (m^3)";
  Real V_total "total instantaneous vol (m^3)";
  Real Vdot "dV_total/dt (m^3/s)";
  Real p "instantaneous pressure (Pa)";
  Real T_mix "volume-weighted mixed temp (K)";

  // pressure perturbation (linearized)
  Real p_pert "small pressure perturbation about mean (Pa)";

  // pressure-factor coefficients used in linearized mechanical eqns
  Real a21(start=0.0), a22(start=0.0), a23(start=0.0), a24(start=0.0);
  Real a41(start=0.0), a43(start=0.0);

  // instantaneous powers & heat flows
  Real P_pv_inst "instantaneous p*Vdot (W)";
  Real Qh_inst "hot-side heat flux (W)";
  Real Qc_inst "cold-side heat flux (W)";
  Real Qrh_inst "hot->reg heat (W)";
  Real Qrc_inst "cold->reg heat (W)";

  // integrals / averages
  Real E_elec_int(start=0.0);
  Real P_elec_avg "cycle-averaged electrical power (W)";

  // small references
  parameter Real Thot_ref = 800.0 "reference hot temp for linearization";

  // 3) INITIAL EQUATION - compute an equilibrium consistent initial point
initial equation
  // start velocities = 0
  v_p = 0.0;
  v_d = 0.0;
  i_alt = 0.0;

  // set thermal initial guesses consistent with KRUSTY environment
  T_hg = 1123.0; // hot boundary temp (K), typical ASC heater-head
  T_cg = 363.0;  // cold sink temp (K)

  // volumes at mid-stroke
  V_e = V_e0 + A_d * x_d;
  V_c = V_c0 + A_p * x_p - A_d * x_d;
  V_total = V_e + V_c + V_reg + V_dead;

  // initial pressure consistent with mean charge (we use p_mean as our initial)
  p = p_mean;

  // initialize simple pressure-factor guesses (tune later with proper linearization)
  a21 = 0.0; a22 = -1.0e3; a23 = -1.0e6; a24 = -1.0e3;
  a41 = 0.0; a43 =  1.0e6;

  // integrators
  E_elec_int = 0.0;

  // regenerator initial temp
  T_r = (T_hg + T_cg)/2.0;

  // ensure m_gas computed from p_mean and Vtot0 is consistent (parameter already computed)

  // end initial equation

  // 4) EQUATIONS
equation
  // --- volumes & kinematics ---
  V_e = V_e0 + A_d * x_d;
  V_c = V_c0 + A_p * x_p - A_d * x_d;
  V_total = V_e + V_c + V_reg + V_dead;
  Vdot = A_p * v_p + A_d * v_d;

  // --- mixture temperature and pressure (single-pressure assumption) ---
  T_mix = (V_e * T_hg + V_c * T_cg) / max(V_e + V_c, 1e-12);
  // ideal gas: p = m_gas * R * T_mix / V_total
  p = m_gas * R * T_mix / max(V_total, 1e-12);

  // --- linearized pressure perturbation used in mechanical linear equations ---
  // simple linear proxy: p_pert = Kp_th*(T_hg - Thot_ref) + Kp_xd * x_d + Kp_xp * x_p
  // (these Kp_* are implicit in the a21.. coefficients in SDM - here we keep structure simple)
  // We'll fold these into "a21" when using full SDM linearization. For now set p_pert = p - p_mean:
  p_pert = p - p_mean;

  // --- mechanical dynamics (Regan linear-style, with pressure-factor coefficients) ---
  der(x_d) = v_d;
  der(x_p) = v_p;

  // displacer acceleration: uses linear pressure-factor couplings (a21..a24)
  der(v_d) = ( a21 * p_pert + a22 * v_d + a23 * x_p + a24 * v_p - c_d * v_d - k_d * x_d ) / m_d;

  // piston acceleration: includes alternator reaction force as -B_elec * v_p (linear damping)
  der(v_p) = ( a41 * p_pert + a43 * x_d - c_p * v_p - k_p * x_p - B_elec * v_p ) / m_p;

  // --- electrical alternator (L-R circuit) ---
  // induced emf = K_e * v_p  =>  L di/dt + (R_alt + R_load) i = emf
  der(i_alt) = ( K_e * v_p - (R_alt + R_load) * i_alt ) / L_alt;

  // --- thermal balances: two-node + regenerator (lumped UA exchanges) ---
  // hot node heat in from heater
  Qh_inst = UA_e * (1123.0 - T_hg); // heater head at 1123 K for KRUSTY baseline
  // cold sink heat out
  Qc_inst = UA_c * (363.0 - T_cg);  // cold sink ~363 K for KRUSTY baseline
  // regenerator coupling
  Qrh_inst = UA_r * (T_hg - T_r);
  Qrc_inst = UA_r * (T_cg - T_r);

  // pV work apportioned by node volume fractions (approx)
  P_pv_inst = p * Vdot;

  // expansion node energy ODE (approximate lumped form)
  der(T_hg) = ( Qh_inst - Qrh_inst - p * Vdot * (V_e / max(V_e + V_c,1e-12)) ) / ( m_gas * (1.5*R) * (V_e / max(V_e + V_c,1e-12)) );

  // compression node energy ODE
  der(T_cg) = ( - Qrc_inst - p * Vdot * (V_c / max(V_e + V_c,1e-12)) + Qc_inst ) / ( m_gas * (1.5*R) * (V_c / max(V_e + V_c,1e-12)) );

  // regenerator thermal capacitance
  der(T_r) = ( Qrh_inst + Qrc_inst ) / C_r;

  // --- electrical power & integrator for averaged power ---
  // instantaneous electrical power dissipated in load (approx)
  // here we use electrical circuit: P_load = i^2 * R_load
  P_pv_inst = p * Vdot; // already defined above
  E_elec_int = der(E_elec_int); // ? der(E_elec_int) : 0; // placeholder to satisfy DAE checks (no-op)
  der(E_elec_int) = i_alt * i_alt * R_load;
  P_elec_avg = E_elec_int / max(time, eps_time);

  // --- per-cycle pressure-factor update (event) ---
  // simple zero-cross detection on displacer velocity to update linear pressure factors
  when edge(v_d) then
    // Very simple heuristic update: scale a23/a43 with amplitudes (placeholders)
    a23 = -1.0e6 * sign(v_p + 1e-12);
    a43 =  1.0e6 * sign(v_d + 1e-12);
    // a21,a22,a24,a41 left as nominal unless you implement full SDM linearization
  end when;

  annotation (
    experiment(StartTime=0, StopTime=1.0, Tolerance=1e-6, Interval=1e-4),
    Documentation(info="
      Regan-style KRUSTY-flavoured linear Stirling model.
      - Replace the per-cycle pressure-factor heuristic with your SDM linearization routine for fidelity.
      - Parameter values (volumes, areas, p_mean) set to KRUSTY/ASC-style starting points; tune UA and C_r to match data.
    ")
  );
end ReganLinearStirling;
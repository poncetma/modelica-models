model FreePistonStirling_Phasing
  "Minimal two-node FPSE Model (no initial-equation, no regenerator thermal mass).
   - Expansion (hot) and Compression (cold) gas nodes (T_hg, T_cg)
   - Regenerator present only as a void volume (V_reg) used for mass/pressure allocation
   - Forchheimer (linear+quadratic) pressure-drop across regenerator (DeltaP_reg)
   - Single algebraic pressure computed from ideal-gas partition (no extra mass eqns)
   - Piston & displacer dynamics, alternator (L-R + load)
   - Starts from variable start=... values; no steady-state solve attempted
  "

  // ---------------- PARAMETERS ----------------
  parameter Real R = 2077.0 "J/(kg K) (He)";
  parameter Real cv = 1.5*R "J/(kg K)";

  // geometry & volumes (mid-stroke)
  parameter Real V_e0  = 4.0e-5 "expansion mid-stroke vol (m^3)";
  parameter Real V_c0  = 6.0e-5 "compression mid-stroke vol (m^3)";
  parameter Real V_reg  = 2.0e-5 "regenerator void vol (m^3)";
  parameter Real V_dead = 1.0e-5 "other dead vol (m^3)";

  parameter Real A_p = 1.539e-4 "power piston area (m^2)";
  parameter Real A_d = 4.91e-4  "displacer area (m^2)";

  // charge conditions (KRUSTY / ASC)
  parameter Real p_mean = 3.52e6 "mean charge pressure (Pa)";
  parameter Real T_guess = 700.0 "reference temp for gas mass estimate (K)";
  parameter Real Vtot0 = V_e0 + V_c0 + V_reg + V_dead "nominal total vol";
  parameter Real m_gas = p_mean * Vtot0 / (R * T_guess) "fixed gas mass (kg)";

  // mechanical
  parameter Real m_p = 0.28 "piston mass (kg)";
  parameter Real m_d = 0.18 "displacer mass (kg)";
  parameter Real k_p = 100.0 "piston spring (N/m)";
  parameter Real k_d = 20000.0 "displacer spring (N/m)";
  parameter Real c_p = 0.0  "piston viscous damping (N*s/m)";
  parameter Real c_d = 1   "displacer viscous damping (N*s/m)";

  // alternator + circuit
  parameter Real K_e = 10.0 "EMF constant (V per m/s)";
  parameter Real L_alt = 1e-3 "H";
  parameter Real R_alt = 0.5 "Ohm (coil)";
  parameter Real R_load = 10.0 "Ohm (external)";

  // regenerator Forchheimer coefficients (tune from geometry)
  parameter Real K_flow1 = 10e3 "linear coeff (Pa·s/kg)";
  parameter Real K_flow2 = 1e5 "quadratic coeff (Pa·s^2/kg^2)";

  // heat transfer (lumped UAs)
  parameter Real UA_e = 200.0 "hot-side UA (W/K)";
  parameter Real UA_c = 80.0  "cold-side UA (W/K)";

  // boundary temps
  parameter Real T_hot = 1123.0 "heater temperature (K)";
  parameter Real T_cold = 363.0 "cold sink temp (K)";

  parameter Real eps = 1e-12 "small number to avoid division by zero";

  // ---------------- STATES ----------------
  // Mechanics
  Real x_p(start=0.0) "piston displacement (m), + increases compression vol";
  Real v_p(start=0.0) "piston velocity (m/s)";
  Real x_d(start=1e-3) "displacer displacement (m), + increases expansion vol";
  Real v_d(start=0.0) "displacer velocity (m/s)";

  // Electrical
  Real i_coil(start=0.0) "alternator current (A)";

  // Thermodynamics (two gas nodes)
  Real T_hg(start=T_hot)  "hot/expansion gas temp (K)";
  Real T_cg(start=T_cold) "cold/compression gas temp (K)";

  // ---------------- ALGEBRAIC / DERIVED ----------------
  Real V_e "instantaneous expansion volume (m^3)";
  Real V_c "instantaneous compression volume (m^3)";
  Real V_total "total instantaneous vol (m^3)";
  Real Vdot "dV_total/dt (m^3/s)";

  Real T_mix "approx mixture temp (K)";
  Real p "instantaneous pressure (Pa)";

  Real Vdot_reg "equivalent volumetric flow through regenerator (m^3/s)";
  Real rho_mix "approx mixture density (kg/m^3)";
  Real mdot_reg "mass flow through regenerator (kg/s)";
  Real DeltaP_reg "pressure drop across regenerator (Pa)";
  Real p_eff "effective pressure acting on mechanical elements (Pa)";

  Real Qh_inst "hot-side heat into expansion node (W)";
  Real Qc_inst "cold-side heat out of compression node (W)";

  Real P_pv_inst "instantaneous p*Vdot (W)";

  Real P_elec_inst "instantaneous electrical power (W)";
  Real E_elec_int(start=0.0) "integrated electrical energy (J)";
  Real P_elec_avg "averaged electrical power (W)";

equation
  // ---------- volumes & kinematics ----------
  V_e = V_e0 + A_d * x_d;
  V_c = V_c0 - A_d * x_d + A_p * x_p;
  V_total = V_e + V_c + V_reg + V_dead;

  Vdot = der(V_total); //A_p * v_p + A_d * v_d;

  // ---------- approximate mixture temperature & density ----------
  // Note: regenerator has no temperature state; mix T from the two gas nodes only
  T_mix = (V_e * T_hg + V_c * T_cg) / max(V_e + V_c, eps);
  rho_mix = p / (R * max(T_mix, 1.0));

  // ---------- pressure (single algebraic ideal-gas relation) ----------
  // Using total gas mass m_gas and a volume-weighted mixture temperature:
  //   p = m_gas * R * T_mix / V_total
  p = m_gas * R * T_mix / max(V_total, eps);

  // ---------- regenerator volumetric & mass flow (approx) ----------
  // project the global Vdot into the regenerator region (fractional apportionment)
  Vdot_reg = Vdot * ( V_reg / max(V_e + V_c, eps) );
  mdot_reg = rho_mix * Vdot_reg;

  // Forchheimer-like pressure drop across regenerator (mdot in kg/s)
  DeltaP_reg = K_flow1 * mdot_reg + K_flow2 * mdot_reg * abs(mdot_reg);

  // effective pressure used by mechanics (apportion drop)
  p_eff = p - DeltaP_reg * ( V_reg / max(V_total, eps) );

  // ---------- mechanical dynamics ----------
  der(x_p) = v_p;
  der(v_p) = ( A_p * (p_eff - p_mean) - c_p * v_p - k_p * x_p - K_e * i_coil ) / m_p;

  der(x_d) = v_d;
  der(v_d) = ( A_d * (p_eff - p_mean) - c_d * v_d - k_d * x_d ) / m_d;

  // ---------- alternator electrical ----------
  der(i_coil) = ( K_e * v_p - (R_alt + R_load) * i_coil ) / L_alt;

  // ---------- thermal: two-node lumped balances ----------
  // heat flows to nodes (simple UA model)
  Qh_inst = UA_e * ( T_hot - T_hg );    // heater -> hot gas
  Qc_inst = UA_c * ( T_cg - T_cold );   // cold gas -> sink (positive if gas hotter than sink)

  // allocate pV work between nodes by instantaneous volumes (pragmatic approximation)
  // note: m_gas * (V_e/(V_e+V_c)) is used as effective "mass in node" for heat capacity
  P_pv_inst = p * Vdot;

  der(T_hg) = ( Qh_inst - p * Vdot * (V_e / max(V_e + V_c, eps)) ) / ( m_gas * (V_e / max(V_e + V_c, eps)) * cv );

  der(T_cg) = ( - Qc_inst - p * Vdot * (V_c / max(V_e + V_c, eps)) ) / ( m_gas * (V_c / max(V_e + V_c, eps)) * cv );

  // ---------- electrical power bookkeeping ----------
  P_elec_inst = i_coil * i_coil * R_load;
  der(E_elec_int) = P_elec_inst;
  P_elec_avg = E_elec_int / max(time, 1e-6);

  annotation (experiment(StartTime=0, StopTime=1.0, Tolerance=1e-6, Interval=1e-4));

end FreePistonStirling_Phasing;
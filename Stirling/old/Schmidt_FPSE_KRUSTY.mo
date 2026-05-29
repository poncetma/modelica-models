model Schmidt_FPSE_KRUSTY
  "Schmidt-style two-node FPSE model adapted for KRUSTY/ASC-E2.
   - Isothermal hot (expansion) and cold (compression) gas nodes (Schmidt assumption)
   - Displacer motion prescribed sinusoidally (control swept volume & phase)
   - Power piston dynamic, alternator as viscous damping B_elec
   - No regenerator pressure drop; regenerator represented either as dead volume or as UA (optional)
   - Mass computed from mean pressure and mean volumes under isothermal assumption
   - Suitable for rapid, robust exploration of swept volume / phase space
   References: Schmidt theory & practical FPSE literature, KRUSTY/ASC documentation.
  "

  // --------------------------- PARAMETERS (physical / KRUSTY-flavored) ---------------------------
  parameter Real R = 2077.0 "J/(kg K)  (helium specific gas constant — use working gas)";
  parameter Real cv = 1.5 * R "J/(kg K) (approx for monoatomic He)";

  // KRUSTY reference conditions
  parameter Real p_mean = 3.52e6  "Mean charge pressure (Pa) — ASC/ KRUSTY documented value";
  parameter Real T_h = 1123.0     "Hot (expansion) isothermal temperature (K)";
  parameter Real T_c = 363.0      "Cold (compression) isothermal temperature (K)";

  // geometric / nominal volumes (starting guesses — tune as required)
  // V_e0, V_c0 are *mean* volumes (mid-stroke) of expansion and compression spaces.
  parameter Real V_e0 = 4.0e-5    "Expansion mean volume (m^3) — initial guess (KRUSTY-ish)";
  parameter Real V_c0 = 6.0e-5    "Compression mean volume (m^3) — initial guess";
  parameter Real V_reg = 2.0e-5    "Regenerator void volume (m^3) - assumed isothermal, included in total";
  parameter Real V_dead = 1.0e-5   "Other dead volume (m^3)";

  // piston & displacer geometry (areas)
  parameter Real A_p = 1.539e-4    "Power piston area (m^2)";
  parameter Real A_d = 4.91e-4     "Displacer effective cross-sectional area (m^2)";

  // power piston mechanical
  parameter Real m_p = 0.49        "Power piston mass (kg) - tune for 100 Hz)";
  parameter Real k_p = 2.0e5       "Piston centering spring (N/m) - rough starting value";
  parameter Real c_p = 2.0         "Mechanical viscous damping (N·s/m)";

  // alternator (simple viscous extraction)
  parameter Real B_elec = 0.0      "Electrical damping (N·s/m). start at 0 for measuring v_rms";
  // optionally map B_elec to electrical R via R_tot = K_e^2 / B_elec if needed.

  // displacer prescribed motion
  parameter Real Xd = 0.002        "Displacer amplitude (m) (2 mm starting guess)";
  parameter Real f_d = 100.0       "Prescribed displacer frequency (Hz)";
  parameter Real omega = 2*Modelica.Constants.pi*f_d;
  parameter Real phi = 1.05        "Phase offset of compression vs expansion (rad) - tune for p·V work";

  parameter Real eps = 1e-12;

  // --------------------------- STATES & VARIABLES ---------------------------
  // power piston mechanical states
  Real x_p(start=0.0) "power piston displacement (m)";
  Real v_p(start=0.0) "power piston velocity (m/s)";

  // gas temperatures are constant (isothermal Schmidt model)
  // instantaneous displacer position (prescribed)
  Real x_d "displacer position (m) - algebraic (prescribed sinusoid)";
  Real v_d "displacer velocity (m/s) - algebraic";

  // algebraic / derived quantities
  Real V_e "instantaneous expansion volume (m^3)";
  Real V_c "instantaneous compression volume (m^3)";
  Real V_total "total instantaneous volume";
  Real Vdot_e "dV_e/dt";
  Real Vdot_c "dV_c/dt";
  Real p "instantaneous pressure (Pa)";
  Real m_gas "total gas mass (kg)";

  Real P_pv_inst "instantaneous p·Vdot power (W)";
  Real P_elec_inst "instantaneous electrical extraction (W)";
  Real P_mech_loss_inst "mechanical viscous loss (W)";

  // bookkeeping for averages (simple integrals / user can postprocess)
  Real E_elec_int(start=0) "integrated electrical energy (J)";
  //Real E_act_int(start=0)  "integrated actuator energy (J) (work to drive displacer)";

equation
  // --------------------------- prescribed displacer (sinusoid) ---------------------------
  x_d = Xd * sin(omega * time + phi);
  v_d = der(x_d); // Modelica computes derivative of sinusoid algebraically

  // --------------------------- instantaneous volumes ---------------------------
  // Consistent sign convention: positive x_d increases expansion volume (moves away from hot)
  V_e = V_e0 + A_d * x_d;
  // compression volume decreases with positive x_d (displacer moves into cold space)
  V_c = V_c0 - A_d * x_d + A_p * x_p; // includes piston contribution to compression volume
  V_total = max(V_e + V_c + V_reg + V_dead, eps);

  Vdot_e = A_d * v_d;
  Vdot_c = -A_d * v_d + A_p * v_p;

  // --------------------------- mass consistent with mean pressure (isothermal) ---------------------------
  // compute total gas mass from mean pressure and mean geometry (Schmidt-style mean)
  m_gas = p_mean * ( V_e0/(R*T_h) + V_c0/(R*T_c) + V_reg/(R*((T_h+T_c)/2)) + V_dead/(R*T_c) );

  // instantaneous pressure from ideal gas partition (isothermal nodes)
  // p = m_gas / (Sum_i (V_i/(R*T_i)))  -> direct algebraic solve
  p = m_gas / ( V_e/(R*T_h) + V_c/(R*T_c) + V_reg/( R * ((T_h + T_c)/2.0) ) + V_dead/(R*T_c) );

  // --------------------------- piston mechanics (power piston free) ---------------------------
  der(x_p) = v_p;
  der(v_p) = ( A_p * (p - p_mean)    // pressure force relative to mean
               - c_p * v_p
               - k_p * x_p
               - B_elec * v_p        // electrical damping (simple)
             ) / m_p;

  // --------------------------- power bookkeeping ---------------------------
  P_pv_inst = p * (Vdot_e + Vdot_c); // p * dV_total/dt (positive when gas does work on surroundings)
  P_elec_inst = B_elec * v_p * v_p;  // instantaneous mechanical power removed by alternator
  P_mech_loss_inst = c_p * v_p * v_p;

  der(E_elec_int) = P_elec_inst;
  // compute actuator force required to produce prescribed displacer motion:
  // F_act = m_d * a_d + c_d * v_d + k_d * x_d + A_d * (p - p_mean)
  // The model does not need explicit m_d/c_d/k_d unless you want actuator energy; user can compute externally if desired.

annotation (experiment(StartTime=0, StopTime=0.5, Tolerance=1e-6, Interval=1e-4));
end Schmidt_FPSE_KRUSTY;
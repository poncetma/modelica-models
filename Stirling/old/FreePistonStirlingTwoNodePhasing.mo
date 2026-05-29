model FreePistonStirlingTwoNodePhasing
"Simplified two-node FPSE model:
   - No initial-equation steady-state solving
   - No regenerator thermal mass
   - Two gas thermal nodes (expansion, compression)
   - Regenerator only provides void volume + pressure drop
   - Linear alternator, L-R electrical circuit
   - Mechanics: piston + displacer
   - Basic UA heat flows to hot/cold boundaries"

  // ---------- PARAMETERS ----------
  parameter Real R = 2077.0 "J/(kg K) helium gas constant";
  parameter Real cv = 1.5 * R "J/(kg K) helium specific heat at constant volume";

  // Geometry & volumes (mid-stroke)
  parameter Real V_e0 = 4.0e-5 "Expansion mid-stroke volume (m3)";
  parameter Real V_c0 = 6.0e-5 "Compression mid-stroke volume (m3)";
  parameter Real V_reg = 2.0e-5 "Regenerator void volume (m3)";
  parameter Real V_dead = 1.0e-5 "Other dead volume (m3)";

  parameter Real A_p = 1.539e-4 "Power piston area (m2)";
  parameter Real A_d = 4.91e-4  "Displacer face area (m2)";

  // Gas charge
  parameter Real p_mean = 3.52e6 "Mean charge pressure (Pa)";
  parameter Real T_guess = 700 "Initial estimate for mean gas temperature (K)";
  parameter Real Vtot0 = V_e0 + V_c0 + V_reg + V_dead;
  parameter Real m_gas = p_mean * Vtot0 / (R * T_guess);

  // Mechanical
  parameter Real m_p = 0.28 "Piston mass (kg)";
  parameter Real m_d = 0.18 "Displacer mass (kg)";
  parameter Real k_p = 800 "Piston spring (N/m)";
  parameter Real k_d = 200 "Displacer spring (N/m)";
  parameter Real c_p = 15 "Piston damping (N*s/m)";
  parameter Real c_d = 5  "Displacer damping (N*s/m)";

  // Alternator
  parameter Real K_e = 10 "EMF constant (V/(m/s))";
  parameter Real L_alt = 1e-3 "Inductance (H)";
  parameter Real R_alt = 0.5 "Coil resistance (ohm)";
  parameter Real R_load = 10 "Load resistance (ohm)";

  // Flow resistance coefficients
  parameter Real K_flow1 = 1000 "Linear flow loss coefficient";
  parameter Real K_flow2 = 1e5 "Quadratic flow loss coefficient";

  // Heat transfer (UA terms)
  parameter Real UA_e = 200 "Hot-side UA";
  parameter Real UA_c = 80  "Cold-side UA";

  // Boundary temperatures
  parameter Real T_hot = 1123.0 "Heater head temperature (K)";
  parameter Real T_cold = 363.0 "Cooler temperature (K)";

  parameter Real eps = 1e-12 "Numerical small value";

  // ---------- STATES ----------
  // Mechanics
  Real x_p(start=0) "Piston position";
  Real v_p(start=0) "Piston velocity";
  Real x_d(start=1e-3) "Displacer position";
  Real v_d(start=0)   "Displacer velocity";

  // Electrical
  Real i_coil(start=0) "Alternator current";

  // Thermodynamics
  Real T_hg(start=T_hot)  "Expansion gas temperature";
  Real T_cg(start=T_cold) "Compression gas temperature";

  // ---------- ALGEBRAIC VARIABLES ----------
  Real V_e, V_c, V_total, Vdot;
  Real m_e, m_c, m_r;
  Real T_mix, rho_mix;
  Real mdot_reg, Vdot_reg;
  Real DeltaP_reg, p, p_eff;

  Real Qh_inst, Qc_inst;
  Real P_elec_inst;
  Real E_elec_int(start=0);
  Real P_elec_avg;

equation
  // --- Volumes ---
  V_e = V_e0 + A_d * x_d;
  V_c = V_c0 + A_p * x_p - A_d * x_d;
  V_total = V_e + V_c + V_reg + V_dead;

  Vdot = A_p * v_p + A_d * v_d;

  // --- Mass partition ---
  m_e = p * V_e  / (R * max(T_hg,1));
  m_r = p * V_reg / (R * max(T_mix,1));
  m_c = p * V_c  / (R * max(T_cg,1));

  // Total mass must equal charge mass
  m_e + m_r + m_c = m_gas;

  // --- Mixture properties used for mdot estimation ---
  T_mix = (V_e*T_hg + V_c*T_cg) / max(V_e+V_c, eps);
  rho_mix = p / (R * max(T_mix,1));

  // --- Regenerator mass flow (projected from Vdot) ---
  Vdot_reg = Vdot * (V_reg / max(V_e + V_c, eps));
  mdot_reg = rho_mix * Vdot_reg;

  // --- Pressure drop across regenerator ---
  DeltaP_reg = K_flow1*mdot_reg + K_flow2*mdot_reg*abs(mdot_reg);

  p_eff = p - DeltaP_reg * (V_reg / max(V_total, eps));

  // --- Mechanics ---
  der(x_p) = v_p;
  der(v_p) = ( A_p*(p_eff - p_mean) - c_p*v_p - k_p*x_p - K_e*i_coil ) / m_p;

  der(x_d) = v_d;
  der(v_d) = ( A_d*(p_eff - p_mean) - c_d*v_d - k_d*x_d ) / m_d;

  // --- Alternator ---
  der(i_coil) = ( K_e*v_p - (R_alt + R_load)*i_coil ) / L_alt;

  // --- Pressure from mass conservation ---
  p = m_gas / ( V_e/(R*T_hg) + V_reg/(R*T_mix) + V_c/(R*T_cg) );

  // --- Heat transfer to boundary plates ---
  Qh_inst = UA_e * (T_hot - T_hg);
  Qc_inst = UA_c * (T_cg - T_cold);

  // --- Thermal ODEs ---
  der(T_hg) = ( Qh_inst - p*Vdot*(V_e/max(V_e+V_c,eps)) ) / max(m_e*cv, eps);
  der(T_cg) = ( -Qc_inst - p*Vdot*(V_c/max(V_e+V_c,eps)) ) / max(m_c*cv, eps);

  // --- Electrical power ---
  P_elec_inst = i_coil*i_coil*R_load;

  der(E_elec_int) = P_elec_inst;
  P_elec_avg = E_elec_int / max(time,1e-6);

  annotation(
    experiment(StartTime=0,StopTime=1,Tolerance=1e-6,Interval=1e-4)
  );

end FreePistonStirlingTwoNodePhasing;
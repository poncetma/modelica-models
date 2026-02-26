model FreePistonStirlingTwoNode_Power
  // ============================================================
  // PARAMETERS
  // ============================================================

  // Geometry
  parameter Real A_p = 0.001 "Piston area (m2)";
  parameter Real A_d = 0.001 "Displacer hot face area (m2)";

  parameter Real V_c0 = 3e-5 "Cold dead volume (m3)";
  parameter Real V_e0 = 2e-5 "Hot dead volume (m3)";

  // Gas
  parameter Real R = 2077.0 "Gas constant for He (J/kg/K)";
  parameter Real m_gas = 0.005 "Total gas mass (kg)";

  // Temperatures fixed at boundaries
  parameter Real T_hot = 1050 "Hot wall temperature (K)";
  parameter Real T_cold = 525 "Cold wall temperature (K)";

  // Mechanical
  parameter Real m_p = 0.2 "Piston mass (kg)";
  parameter Real m_d = 0.1 "Displacer mass (kg)";
  parameter Real k_p = 3500 "Piston centering spring (N/m)";
  parameter Real k_d = 1800 "Displacer centering spring (N/m)";
  parameter Real c_p = 40 "Piston damping (N·s/m)";
  parameter Real c_d = 25 "Displacer damping (N·s/m)";

  // Electrical (linear alternator)
  parameter Real R_coil = 6.0 "Coil resistance (ohm)";
  parameter Real K_e = 12.0 "Back-EMF constant (V·s/m)";

  // ============================================================
  // STATE VARIABLES
  // ============================================================

  Real x_p(start=0), v_p(start=0) "Piston position & velocity";
  Real x_d(start=0), v_d(start=0) "Displacer position & velocity";

  Real i_coil(start=0) "Alternator current";

  Real p(start=2e6) "Mean pressure";

  // Thermal states (gas nodes)
  Real T_hg(fixed=false) "Hot gas node temperature";
  Real T_cg(fixed=false) "Cold gas node temperature";

  // Derived volumes
  Real V_e "Expansion volume";
  Real V_c "Compression volume";

  // Mass in each zone
  Real m_e, m_c;

  // Heat capacities (simplified)
  parameter Real Cp = 5193 "He Cp (J/kg/K)";
  parameter Real tau = 0.25 "Characteristic thermal time (s)";

equation
  // ============================================================
  // VOLUMES
  // ============================================================

  V_e = V_e0 + A_d*x_d;
  V_c = V_c0 + A_p*x_p - A_d*x_d;

  // ============================================================
  // THERMODYNAMIC RELATIONS
  // ============================================================

  m_e = p*V_e/(R*T_hg);
  m_c = p*V_c/(R*T_cg);
  m_e + m_c = m_gas;

  // Simple thermal dynamics toward wall temperatures
  der(T_hg) = (T_hot - T_hg)/tau;
  der(T_cg) = (T_cold - T_cg)/tau;

  // ============================================================
  // MECHANICS
  // ============================================================

  der(x_p) = v_p;
  der(x_d) = v_d;

  der(v_p) =
    ( A_p*p - k_p*x_p - c_p*v_p - K_e*i_coil )/m_p;

  der(v_d) =
    ( A_d*p - k_d*x_d - c_d*v_d )/m_d;

  // ============================================================
  // ALTERNATOR ELECTRICAL
  // ============================================================

  der(i_coil) = ( -R_coil*i_coil - K_e*v_p )/0.01;

  // ============================================================
  // INITIAL EQUILIBRIUM
  // ============================================================

initial equation
  // At equilibrium, velocity = 0
  v_p = 0;
  v_d = 0;
  i_coil = 0;

  // Gas temperatures equal wall temps at steady state
  T_hg = T_hot;
  T_cg = T_cold;

  // Static force balance (defines x_p0, x_d0, p0)
  A_p*p - k_p*x_p = 0;
  A_d*p - k_d*x_d = 0;

  // Mass constraint handled automatically
end FreePistonStirlingTwoNode_Power;
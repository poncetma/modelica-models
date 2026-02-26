model FreePistonStirling

  // --- Gas constants ---
  parameter Real R = 8.314; //J/mol/K
  parameter Real M_He = 0.004; //molar mass
  parameter Real R_He = R/M_He; //Shorthand to use PV = (m)RT
  parameter Real cv = 1.5*R_He; //@15 C, 1 bar

  // --- Boundary temperatures ---
  parameter Real T_hot = 1123;
  parameter Real T_cold = 363;

  // --- Volumes --- Are these remotely accurate?
  parameter Real V_e0 = 4e-5;
  parameter Real V_c0 = 6e-5;
  parameter Real V_reg = 2e-5;
  parameter Real V_dead = 1e-5;

  // Areas
  parameter Real A_p = 1.5e-4;
  parameter Real A_d = 5.0e-4;

  // Mechanical -- Set all damping coeffs to 0 at first.
  parameter Real m_p = 0.28;
  parameter Real m_d = 0.18;
  parameter Real c_p = 0.1;//15
  parameter Real c_d = 0.1;//5
  parameter Real k_p = 800;
  parameter Real k_d = 800;
  parameter Real B_elec = 0.0;

  // UA
  parameter Real UAh = 200;
  parameter Real UAc = 80;
  parameter Real UAr = 300;

  // Regenerator capacitance
  //parameter Real C_r = 1000;

  // Gas charge
  parameter Real p_ref = 2e6; //20 bar?
  parameter Real T_guess = 700; //0.5*(T_hot + T_cold);
  parameter Real Vtot0 = V_e0 + V_c0 + V_reg + V_dead;
  parameter Real m_gas = p_ref*Vtot0/(R_He*T_guess);  

  // States
  Real x_p(start=0), v_p(start=0); //Give a slight perturbation here
  Real x_d(start=0), v_d(start=0);

  // Temps
  Real T_g(start=T_guess);
  //Real T_r(start=T_guess);

  // Volumes
  Real V_e, V_c, V_total, Vdot;

  // Pressure + mixture temperature
  Real p;  

  // Heat flows
  Real Qh, Qc;
  //Real Qrh;
  

equation
  // --- Kinematics ---
  der(x_p) = v_p;
  der(x_d) = v_d;

  // --- Volumes ---
  V_e = max(V_e0 + A_d*x_d, 1e-15);
  V_c = max(V_c0 - A_p*x_p - A_d*x_d, 1e-15);
  V_total = V_e + V_c + V_reg + V_dead;
  Vdot   = -A_p*v_p;

  // --- Pressure and temperature mixing ---
  p = m_gas*R_He*T_g / V_total;

  // --- Mechanical ---
  der(v_p) = ( A_p*(p - p_ref) - c_p*v_p - k_p*x_p - B_elec*v_p)/m_p;
  //der(v_d) = ( A_d*(p - p_ref) - c_d*v_d - k_d*x_d )/m_d;
  der(v_d) = ( A_d*(p - p_ref) - c_d*v_d - k_d*x_d )/m_d; //force from gas should be opposite?

  // --- Thermal ---
  Qh = UAh*(T_hot - T_g);
  Qc = UAc*(T_cold - T_g);
  //Qrh = UAr*(T_r - T_g);
  

  //der(T_g) = (Qh + Qc - Qrh - p*Vdot)/(m_gas*cv);
  der(T_g) = (Qh + Qc - p*Vdot)/(m_gas*cv);
  //der(T_r) = (Qrh)/C_r;
    
  annotation(
    experiment(StartTime=0, StopTime=0.5, Tolerance=1e-7, Interval=1e-4)
  );


end FreePistonStirling;
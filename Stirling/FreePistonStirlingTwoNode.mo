model FreePistonStirlingTwoNode

  /* Model of a free-piston Stirling engine with two connected gas volumes and an ideal regenerator */

  // --- Gas constants ---
  parameter Real R = 8.314; //J/mol/K
  parameter Real M_He = 0.004; //molar mass
  parameter Real R_He = R/M_He; //Shorthand to use PV = mRT
  parameter Real cv = 1.5*R_He; //@15 C, 1 bar

  // --- Boundary temperatures ---
  parameter Real T_hot = 850 + 273.0; 
  //850 C is the design temp, real KRUSTY condition was 650. Should tune the system based on the designed temp + power output.
  parameter Real T_cold = 363; //90 C spec in official ASC testing (NASA)


  // -- Area estimates 
  parameter Real D_d = 0.065/2.0; //Displacement cylinder outer diameter (visual estimate from diagram)
  parameter Real A_d = 3.14159*(D_d)^2./4.; //Leaving room for heat exchangers in the displacement cylinder 
  parameter Real A_p = A_d*0.7; //Shared outer diameter with the displacer shaft going through the power piston. 

  // --- Volume estimates
  //parameter Real V_e0 = 4e-5;
  //parameter Real V_c0 = 6e-5;  
  parameter Real V_swept_d = A_d*(0.1864/2.0)*0.5; //Volume swept by displacer piston, crude estimate
  parameter Real V_swept_p = 0.25*V_swept_d; //Crude estimate
  //Baseline volumes corresponding to mid-stroke piston positions (x=0)
  parameter Real V_e_dead = V_swept_d*0.2; 
  parameter Real V_e0 = V_swept_d*0.5 + V_e_dead; 
  //The dead volume on the power piston side seems very low (call it 10% of the displacer side dead volume)
  parameter Real V_c_dead = V_swept_p*0.1;
  parameter Real V_c0 = V_swept_d*0.5 + V_swept_p*0.5 + V_c_dead; 
  parameter Real V_reg = V_e_dead;
  parameter Real V_dead_other = V_reg*0.1;  
  
  // Mechanical -- Set all damping coeffs to 0 at first.
  //Note the total engine mass is 0.975 kg. 
  parameter Real m_p = 0.25;
  parameter Real m_d = 0.15;
  parameter Real c_p = 0.01;//15--these higher damping values really knock down the power output
  parameter Real c_d = 0.001; //
  parameter Real k_p = 1000.0; 
  parameter Real k_d = 200.0;
  parameter Real B_elec = 20.0;

  // UA
  parameter Real UAh = 300;
  parameter Real UAc = 300;
  
  parameter Real UAr = 500; //300;

  // Regenerator capacitance
  parameter Real C_r = 10; //10 J/K is in the right physical ballpark, but it doesn't seem that a regenerator is compatible with a two-node model.
  // Gas charge
  parameter Real p_ref = 3.52e6; //corrected as per Wood et al.  
  parameter Real T_guess = 700; //Tune this value so it roughly equals the average of T_mix. 
  parameter Real Vtot0 = V_e0 + V_c0 + V_reg + V_dead_other;
  parameter Real m_gas = p_ref*Vtot0/(R_He*T_guess);  

  // States
  //x_p,x_d = 0 correspond to mid-stroke
  Real x_p(start=0.00), v_p(start=0.001); //Give a slight perturbation here
  Real x_d(start=0.00), v_d(start=-0.001);

  // Temps
  Real T_hg(start=T_hot);
  Real T_cg(start=T_cold);
  Real T_r(start=750);

  // Volumes
  Real V_e, V_c, V_total, V_c_dot, V_e_dot;
  Real frac_e, frac_c;
  Real Vdot_pp; //power piston displacement rate

  // Pressure + mixture temperature
  Real p;  
  Real T_mix;

  // Heat flows
  Real Qh, Qc;
  //Real Qr;
  Real Qrh, Qrc;
  //Real pV_e, pV_c;
  
  Real P_elec_inst;
  Real E_elec_integral;
  Real P_elec_avg;
  parameter Real eps_time = 1e-9;


equation
  // --- Kinematics ---
  der(x_p) = v_p;
  der(x_d) = v_d;

  // --- Volumes ---
  V_e = max(V_e0 + A_d*x_d, 1e-15);
  V_c = max(V_c0 - A_d*x_d + A_p*x_p , 1e-15); //use a consistent convention where positive displacement of the "local" piston increases the local volume
  V_total = V_e + V_c + V_reg + V_dead_other;
  //Vdot = der(V_total); //keep it general  
  V_e_dot =  der(V_e); //A_d*v_d;
  V_c_dot =  der(V_c); //-A_d*v_d + A_p*v_p;  

  // --- Pressure and temperature mixing ---
  T_mix = (V_e * T_hg + V_c * T_cg) / V_total;
  p = m_gas * R_He * T_mix / V_total; //pressure computed from weighted average temperature
  
  // --- Mechanical ---
  der(v_p) = ( A_p*(p - p_ref) - c_p*v_p - k_p*x_p - B_elec*v_p)/m_p;  //higher pressure should force the power piston in +ve direction (pushing alternator)
  der(v_d) = (-A_d*(p - p_ref) - c_d*v_d - k_d*x_d )/m_d; //higher pressure should push the displacer backward to the dead space

  // --- Thermal ---
  Qh = UAh*(T_hot - T_hg);
  Qc = UAc*(T_cold - T_cg);  
  
  //Qr = UAr*(T_r - 0.5*(T_hg + T_cg)); 
  //Qr = UAr*(T_r - T_hg);   
  Qrh = UAr*(T_r - T_hg);
  Qrc = UAr*(T_r - T_cg);
  
  //der(T_hg) = (Qh - Qr - p*( V_e_dot )) / (m_gas * cv * (V_e / V_total));
  //der(T_cg) = (Qc + Qr - p*( V_c_dot )) / (m_gas * cv * (V_c / V_total));
  frac_e = V_e / V_total;  
  frac_c = V_c / V_total;   
  //der(T_hg) = (Qh + Qrh - p*( V_e_dot )) / (m_gas * cv * frac_e); //account for work done on boundary
  //der(T_cg) = (Qc + Qrc - p*( V_c_dot )) / (m_gas * cv * frac_c);
  Vdot_pp = A_p*v_p;
  der(T_hg) = (Qh + Qrh ) / (m_gas * cv * frac_e); //only account for real piston work on the compression side
  der(T_cg) = (Qc + Qrc - p*( Vdot_pp )) / (m_gas * cv * frac_c);
  
  //der(T_r) = (Qr - Qr)/C_r; //simplified regenerator, just shuttles heat
  der(T_r) = (-Qrh - Qrc)/C_r; //genuine regenerator 
  
  // Electrical
  P_elec_inst = B_elec*v_p*v_p; //instantaneous electrical power = F*v = (B*v)*v
  der(E_elec_integral) = P_elec_inst;
  P_elec_avg = E_elec_integral / max(time, 1e-6);
    
  annotation(
    experiment(StartTime=0, StopTime=10.0, Tolerance=1e-7, Interval=1e-4, Method=cvode)
  );


end FreePistonStirlingTwoNode;
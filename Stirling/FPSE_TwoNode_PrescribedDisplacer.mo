model TwoNode_PrescribedDisplacer
"Two-node free-piston Stirling (minimal/regenerator-simple).
   - Expansion (hot) and compression (cold) gas nodes
   - Regenerator represented as a simple heat-exchange link (UA_r) only
   - NO pressure drop through regenerator
   - Correct sign convention: pressure pushes piston and displacer in opposite directions
   - Work terms in energy balances expressed as p * dV_e/dt and p * dV_c/dt
   - No initial-equation block; use start= values
  "

  // ------------- PARAMETERS -------------
  parameter Real R = 2077.0 "J/(kg K), helium";
  parameter Real cv = 1.5 * R "J/(kg K)";

  // geometric / volume parameters (mid-stroke reference)
  parameter Real D_d = 0.065/2.0; //Displacement cylinder outer diameter (visual estimate from diagram)
  parameter Real A_d = 3.14159*(D_d)^2./4.; //Leaving room for heat exchangers in the displacement cylinder 
  parameter Real A_p = A_d*0.7; //Shared outer diameter with the displacer shaft going through the power piston.   
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

  // charge and derived mass
  parameter Real p_mean = 3.52e6 "mean charge pressure (Pa)";
  parameter Real T_guess = 700; //0.5*(T_hot + T_cold) "reference temp for m_gas calc (K)";
  parameter Real Vtot0 = V_e0 + V_c0 + V_reg + V_dead_other "nominal total vol (m^3)";
  parameter Real m_gas = p_mean * Vtot0 / (R * T_guess) "working gas mass (kg)";

  // mechanical parameters
  parameter Real m_p = 0.25 "power piston mass (kg)";
  parameter Real m_d = 0.025 "displacer mass (kg)";
  parameter Real k_p = 200000 "piston centering spring (N/m)";
  //parameter Real k_d = 100 "displacer centering spring (N/m)";
  parameter Real c_p = 0.00  "piston viscous damping (N·s/m)";
  //parameter Real c_d = 0.0001   "displacer viscous damping (N·s/m)";

  // alternator & electric load
  parameter Real B_elec = 20;
  /*
  parameter Real K_e = 0.01    "EMF constant (V per m/s)";
  parameter Real L_alt = 1e-3  "alternator inductance (H)";
  parameter Real R_alt = 0.0   "coil resistance (Ohm)";
  parameter Real R_load = 1000.0 "external load resistance (Ohm)";
  */

  // heat transfer (including 'regenerator' as UA link)
  parameter Real UA_e = 300.0 "hot-side UA (W/K)";
  parameter Real UA_c = 300.0  "cold-side UA (W/K)";
  parameter Real UA_r = 0.0 "regenerator-equivalent UA (hot<->cold) (W/K)";
  parameter Real C_r = 10 "regenerator capacitance (J/K)";
  
  // boundary temps
  parameter Real T_hot = 1123; //1123.0 "heater temperature (K)";
  parameter Real T_cold = 363.0 "cold sink temperature (K)";

  // numerical
  parameter Real eps = 1e-12 "numerical small number";
  
  
  // Prescribed displacer parameters
  parameter Real Xd = (0.5*V_swept_d)/A_d; //0.025    "displacer amplitude (m) e.g. 2 cm";
  parameter Real f = 105.0     "prescribed frequency (Hz)";
  parameter Real phi = 0.0; //180*3.1415/180.0     "prescribed phase (rad)";
  

  // ------------- STATES -------------
  // mechanical states
  Real x_p(start=0.0) "piston displacement (m); choose sign so + increases compression volume";
  Real v_p(start=0.0) "piston velocity (m/s)";

  Real x_d "displacer displacement (m); choose sign so + increases expansion volume";
  Real v_d "displacer velocity (m/s)";

  // electrical state 
  //Real i_coil(start=0.0) "alternator coil current (A)";

  // thermal states (two gas nodes)
  Real T_hg(start=T_hot) "hot/expansion gas temperature (K)";
  Real T_cg(start=T_cold) "cold/compression gas temperature (K)";
  Real T_r(start = 0.5*(T_hot + T_cold)) "regenerator temperature";

  // bookkeeping for energy/power
  Real E_elec_int(start=0.0) "integrated electrical energy delivered to load (J)";

  // ------------- ALGEBRAIC / DERIVED -------------
  Real V_e "instantaneous expansion volume (m^3)";
  Real V_c "instantaneous compression volume (m^3)";
  Real V_total "total instantaneous volume (m^3)";
  Real Vdot_e "dV_e/dt (m^3/s)";
  Real Vdot_c "dV_c/dt (m^3/s)";
  Real Vdot_pp;
  //Real Vdot "total dV/dt (m^3/s)";

  Real T_mix "volume-weighted mixture temperature (K)";
  Real p "instantaneous pressure (Pa)";

  Real Qh "hot-side heat input (W)";
  Real Qc "cold-side heat input (W)";
  //Real Qr_inst "regenerator-equivalent heat from hot to cold (W)";
  Real Qrh;
  Real Qrc;

  Real P_pv_e "pV work rate associated with expansion node (W)";
  Real P_pv_c "pV work rate associated with compression node (W)";

  Real P_elec_inst "instantaneous electrical power to load (W)";
  Real P_elec_avg "averaged electrical power (W)";
  
  Real frac_e;
  Real frac_c;
  
  Real F_actuator; //Factor in work done by a virtual actuator/controller
  Real P_actuator;

equation
  // -------- volumes & kinematics --------
  V_e = V_e0 + A_d * x_d;                    // expansion node volume
  V_c = V_c0 + A_p * x_p - A_d * x_d;        // compression node volume
  V_total = V_e + V_c + V_reg + V_dead_other;
  // individual volume derivatives (explicit)
  Vdot_e = A_d * v_d; //der(V_e); //
  Vdot_c = A_p * v_p - A_d * v_d; // der(V_c); //
  //Vdot = Vdot_e + Vdot_c;

  // -------- pressure (single algebraic ideal-gas relation) --------
  // form: p = m_gas * R * T_mix / V_total, with T_mix volume-weighted of the two gas nodes
  T_mix = (V_e * T_hg + V_c * T_cg) / max(V_e + V_c, eps);
  p = m_gas * R * T_mix / max(V_total, eps);

  // -------- mechanical dynamics (correct sign convention) --------
  
  // --- Displacer:
  //   coordinate x_d increases upward, away from the hot space. 
  // Prescribing the motion explicitly as a function of time
  x_d = Xd * sin(2*Modelica.Constants.pi*f*time + phi);
  v_d = der(x_d);
  
  // --- Power piston:
  //   coordinate x_p increases so that compression volume rises when x_p increases.
  //   gas pressure pushes the piston in the +x_p direction (force = +A_p * p)
  //   pressure-related restoring around reference: A_p*(p - p_mean) drives dynamics.
  der(x_p) = v_p;
  der(v_p) = ( A_p * (p - p_mean)   // pressure pushing piston in +x_p direction
               - c_p * v_p
               - k_p * x_p
               - B_elec*v_p
               //- K_e * i_coil      // electromagnetic reaction force (opposes motion when current flows)
             ) / m_p;

  
  // -------- alternator electrical circuit (linear L-R + load) --------
  // emf = K_e * v_p  (generator)
  //der(i_coil) = ( K_e * v_p - (R_alt + R_load) * i_coil ) / L_alt;  

  // -------- heat transfers (UA) and regenerator as simple UA link --------
  Qh = UA_e * ( T_hot - T_hg );   // heater → hot gas (positive when T_hot > T_hg)
  Qc = UA_c * ( T_cold - T_cg );  // cold gas → sink (positive when T_cg > T_cold)  
  Qrh = UA_r * (T_r - T_hg);
  Qrc = UA_r * (T_r - T_cg);

  // -------- p·V work of each node explicitly (requested) --------
  // Work done *on* the gas in the expansion node = p * dV_e/dt (note sign convention: if V_e increases, gas does work on surroundings)
  // We follow the convention: in the node energy balance we subtract p * dV_node/dt when it is work done by the node.
  P_pv_e = p * Vdot_e;
  P_pv_c = p * Vdot_c;

  // -------- energy balances for each gas node (explicit work terms) --------
  //
  // Expansion node (hot):
  //   m_e * cv * dT_hg/dt = Qh_inst - P_pv_e - Qr_inst
  //     - Qh_inst: heat input from heater
  //     - P_pv_e: p * dV_e/dt (work associated with expansion node)
  //     - Qr_inst: heat transferred from hot node into regenerator (positive if T_hg > T_cg)
  //
  // Compression node (cold):
  //   m_c * cv * dT_cg/dt = Qr_inst - P_pv_c - Qc_inst
  //     - Qr_inst: heat received from hot node via regenerator-like UA
  //     - P_pv_c: p * dV_c/dt (work associated with compression node)
  //     - Qc_inst: heat rejected to cold sink (positive when gas > sink)
  //
  // We allocate the total gas mass between nodes proportionally to instantaneous volumes:
  //   m_e = m_gas * (V_e/(V_e+V_c));   m_c = m_gas * (V_c/(V_e+V_c))
  // This makes the capacity terms proportional to node volume fraction.
  //
  // (All divisions protected with max(...,eps) to avoid numerical issues)
  
  frac_e = V_e / max(V_e + V_c, eps);  
  frac_c = V_c / max(V_e + V_c, eps);

  // energy ODEs:  
  //der(T_hg) = ( Qh + Qrh - P_pv_e ) / ( (m_gas * frac_e) * cv );
  //der(T_cg) = ( Qc + Qrc - P_pv_c ) / ( (m_gas * frac_c) * cv );
  Vdot_pp = A_p*v_p;
  der(T_hg) = (Qh + Qrh ) / (m_gas * cv * frac_e); //only account for real piston work on the compression side
  der(T_cg) = (Qc + Qrc - p*( Vdot_pp )) / (m_gas * cv * frac_c);
  der(T_r)  = (-Qrh -Qrc)/C_r;

  // -------- electrical power bookkeeping (load power) --------
  P_elec_inst = B_elec*v_p*v_p; //instantaneous electrical power
  //P_elec_inst = i_coil * i_coil * R_load;  
  der(E_elec_int) = P_elec_inst;
  // average power computed externally (user can post-process or use E_elec_int/time)
  P_elec_avg = E_elec_int / max(time, 1e-6);
  
  F_actuator = m_d * der(v_d) - ( -A_d*(p-p_mean) ); 
  P_actuator = abs(F_actuator*v_d); 

  annotation (experiment(StartTime=0, StopTime=1, Tolerance=1e-6, Interval=1e-4));
end TwoNode_PrescribedDisplacer;
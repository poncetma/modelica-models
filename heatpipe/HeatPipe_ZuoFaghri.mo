model HeatPipe_ZuoFaghri
"Six-node thermal network tailored to KRUSTY-style sodium heat pipes.
 Uses Haynes-230 wall material and liquid sodium working fluid (explicit masses & properties).
 Sources: NASA KRUSTY report (Haynes 230, sodium heat pipes), INL microreactor heat-pipe docs,
            Haynes datasheets, MatLib (PNNL) sodium property correlations. See chat citations."

  import Modelica.Constants.pi;

  //----------------------------
  // Geometry (KRUSTY-like)
  //----------------------------
  parameter Real Le = 0.30 "Evaporator length (m) -- KRUSTY-like";
  parameter Real La = 0.40 "Adiabatic length (m)";
  parameter Real Lc = 0.30 "Condenser length (m)";
  parameter Real L = Le + La + Lc;

  // outer diameter ~ 5/8 in = 15.875 mm (value reported in INL microreactor example)
  parameter Real do = 0.015875 "Outer diameter (m) (5/8 in)";
  parameter Real twall = 2e-3 "Wall thickness (m) (typical reactor heat pipe)";
  parameter Real ro = do/2 "Outer radius (m)";
  parameter Real wick_thickness = 1e-3 "Wick thickness nominal (m)";
  parameter Real ri = ro - twall - wick_thickness "Inner radius (m); assume wick annulus ~1 mm thick inside wall";  

  parameter Real A_cond_outer = 2*pi*ro*Lc "condenser outer convective area (m2)";
  parameter Real A_evap_outer = 2*pi*ro*Le "evaporator outer convective area (m2)";
  parameter Real A_wick_outer_e = 2*pi*(ri+wick_thickness)*Le;  
  parameter Real A_wick_outer_c = 2*pi*(ri+wick_thickness)*Lc;    
  parameter Real A_wick_outer_a = 2*pi*(ri+wick_thickness)*La;  
  parameter Real A_wick_inner_e = 2*pi*(ri)*Le;
  parameter Real A_wick_inner_c = 2*pi*(ri)*Lc;
  parameter Real A_wick_inner_a = 2*pi*(ri)*La;
  parameter Real A_wall_inner_e = A_wick_outer_e;
  parameter Real A_wall_inner_c = A_wick_outer_c;
  parameter Real A_wall_inner_a = A_wick_outer_a;

  //----------------------------
  // Materials: Haynes-230 (wall) -- properties representative at ~800 C
  // (Haynes datasheets: k ≈ 24 W/mK at 800C, cp ≈ 590–600 J/kgK, density ≈ 9050 kg/m3)
  //----------------------------
  parameter Real rho_wall = 9050 "kg/m3 (Haynes-230 approx)";
  parameter Real cp_wall  = 595  "J/kg.K (Haynes-230 at ~800 C, datasheet)";
  parameter Real k_wall   = 24   "W/m.K (Haynes-230 at ~800 C)";

  //----------------------------
  // Wick: screen-covered groove (effective properties)
  //----------------------------
  parameter Real eps = 0.60 "wick porosity (assumed)";
  // effective wick solid+liquid properties computed below (need cp/k/rho of void+solid)
  // choose effective solid (screen) assumed steel-like: use stainless-like numbers for screen base:
  parameter Real rho_screen = 8000 "kg/m3 (approx stainless screen)";
  parameter Real cp_screen  = 500  "J/kg.K (screen)";
  parameter Real k_screen   = 20   "W/m.K (screen)";

  // effective thermal conduction of wick (porous + liquid)
  //parameter Real k_wick = 10 "W/m.K (effective porous wick thermal conductivity, conservative)";
  parameter Real k_wick = min(k_screen, k_Na_const) "W/m.K (effective porous wick thermal conductivity, conservative)";

  //----------------------------
  // Liquid sodium property helpers (use MatLib/PNNL correlations)
  // We include simple temperature-dependent functions for rho and cp based on the MatLib forms.
  // Reference: MatLib / PNNL sodium property expressions (Fink & Leibowitz) (see chat citations).
  //----------------------------
  function rho_Na_T
    input Real T "Temperature (K)";
    output Real rho "density of liquid sodium (kg/m3)";
  algorithm
    // simple approximate value in reactor T-range (use MatLib recommended approximations in practice)
    // Use a representative constant for clarity: ~850 kg/m3 at ~1000K (MatLib gives T-dependent formula).
    rho := 850.0;
  end rho_Na_T;

  function cp_Na_T
    input Real T "Temperature (K)";
    output Real cp "isobaric specific heat (J/kg.K)";
  algorithm
    // Use MatLib/PNNL fit: Cp ~ 1658.2 - 0.8479*T + ... is valid up to 2000 K (Fink & Leibowitz)
    // For stability and brevity we use the primary polynomial piece (valid 372K..2000K)
    // cp = 1658.2 - 0.8479*T + 4.4541e-4*T^2 - 2.9926e6 / T^2  (units J/kg-K)
    cp := 1658.2 - 0.8479*T + 4.4541e-4*T*T - 2.9926e6/(T*T);
  end cp_Na_T;

  parameter Real k_Na_const = 70 "W/m.K (approx liquid sodium thermal conductivity at high T)";

  //----------------------------
  // Volumes and masses (explicit)
  //----------------------------
  // annular wick (approximated area at mean radius)
  parameter Real r_wick_mid = ri + wick_thickness/2;
  parameter Real A_wick_radial = pi*( (ri + wick_thickness)^2 - ri^2 ); // annulus XS area 
  parameter Real V_wick_e = A_wick_radial * Le;
  parameter Real V_wick_a = A_wick_radial * La;
  parameter Real V_wick_c = A_wick_radial * Lc;

  // wall volume (cylindrical shell)
  parameter Real A_wall_radial = pi*((ro)^2 - (ro-twall)^2); //XS area
  parameter Real V_wall_e = A_wall_radial * Le;
  parameter Real V_wall_a = A_wall_radial * La;
  parameter Real V_wall_c = A_wall_radial * Lc;

  // masses (solid wall + wick porous mass = screen solid fraction + liquid sodium fraction within wick)
  // approximate wick solid volume fraction = (1 - eps)
  parameter Real m_wall_e = rho_wall * V_wall_e;
  parameter Real m_wall_a = rho_wall * V_wall_a;
  parameter Real m_wall_c = rho_wall * V_wall_c;

  // wick composite mass (solid screen + sodium)
  parameter Real V_wick_solid_e = (1 - eps) * V_wick_e;
  parameter Real V_wick_liquid_e = eps * V_wick_e;
  parameter Real V_wick_solid_a = (1 - eps) * V_wick_a;
  parameter Real V_wick_liquid_a = eps * V_wick_a;
  parameter Real V_wick_solid_c = (1 - eps) * V_wick_c;
  parameter Real V_wick_liquid_c = eps * V_wick_c;

  parameter Real m_wick_solid_e = rho_screen * V_wick_solid_e;
  parameter Real m_wick_liq_e   = rho_Na_T(900) * V_wick_liquid_e; // approx T ~ 900K initial
  parameter Real m_wick_solid_a = rho_screen * V_wick_solid_a;
  parameter Real m_wick_liq_a   = rho_Na_T(900) * V_wick_liquid_a;
  parameter Real m_wick_solid_c = rho_screen * V_wick_solid_c;
  parameter Real m_wick_liq_c   = rho_Na_T(800) * V_wick_liquid_c;

  parameter Real m_wick_e = m_wick_solid_e + m_wick_liq_e;
  parameter Real m_wick_a = m_wick_solid_a + m_wick_liq_a;
  parameter Real m_wick_c = m_wick_solid_c + m_wick_liq_c;

  //----------------------------
  // Thermal capacitances (explicit Cp * mass)
  //----------------------------
  // wall heat capacity
  parameter Real C1 = m_wall_e * cp_wall;  
  parameter Real C4 = m_wall_c * cp_wall;
  parameter Real C6 = m_wall_a * cp_wall;

  // wick heat capacity (composite: solid + liquid sodium)
  parameter Real cp_wick_eff_e = (m_wick_solid_e*cp_screen + m_wick_liq_e*cp_Na_T(900)) / m_wick_e;
  parameter Real cp_wick_eff_a = (m_wick_solid_a*cp_screen + m_wick_liq_a*cp_Na_T(900)) / m_wick_a;
  parameter Real cp_wick_eff_c = (m_wick_solid_c*cp_screen + m_wick_liq_c*cp_Na_T(800)) / m_wick_c;

  parameter Real C2 = m_wick_e * cp_wick_eff_e;  
  parameter Real C3 = m_wick_c * cp_wick_eff_c;
  parameter Real C5 = m_wick_a * cp_wick_eff_a;

  //----------------------------
  // Heat transfer boundary (condenser convection)
  //----------------------------
  parameter Real Qe = 4000.0/8 "W applied heat to evaporator (W) -- KRUSTY-scale reactor ~5 kW thermal divided by 8 heat pipes";
  parameter Real h_cond = 50 "W/m2.K condenser convection HTC (estimate)";
  parameter Real T_inf = 273.15 "cold sink (Stirling hot side )";


  // Define the exact parameters used by Faghri  
  
  parameter Real k_1 = k_wall;
  parameter Real k_2 = k_wick;
  parameter Real k_3 = k_wick;
  parameter Real k_4 = k_wall;
  parameter Real k_5 = k_wick;
  parameter Real k_6 = k_wall;   
  parameter Real A_1 = A_wall_inner_e;
  parameter Real A_2 = A_wick_inner_e; 
  parameter Real A_3 = A_wick_inner_c;
  parameter Real A_4 = A_wall_inner_c;
  parameter Real A_5 = A_wick_radial;
  parameter Real A_6 = A_wall_radial;
  parameter Real L_1 = twall;
  parameter Real L_2 = wick_thickness;
  parameter Real L_3 = wick_thickness;
  parameter Real L_4 = twall;
  parameter Real L_5 = La;
  parameter Real L_6 = La;  
  
  parameter Real rho_wick_e = m_wick_e/V_wick_e;
  parameter Real rho_wick_c = m_wick_c/V_wick_c;
  parameter Real rho_wick_a = m_wick_a/V_wick_a;
  
  parameter Real alpha_1 = k_1/rho_wall/cp_wall;
  parameter Real alpha_2 = k_2/rho_wick_e/cp_wick_eff_e; 
  parameter Real alpha_3 = k_3/rho_wick_c/cp_wick_eff_c;
  parameter Real alpha_4 = k_4/rho_wall/cp_wall;
  parameter Real alpha_5 = k_5/rho_wick_a/cp_wick_eff_a;
  parameter Real alpha_6 = k_6/rho_wall/cp_wall;
  
  
  parameter Real C_1 = k_1*A_1/L_1; //conductances
  parameter Real C_2 = k_2*A_2/L_2;
  parameter Real C_3 = k_3*A_3/L_3;
  parameter Real C_4 = k_4*A_4/L_4;
  parameter Real C_5 = k_5*A_5/L_5;
  parameter Real C_6 = k_6*A_6/L_6;
  parameter Real Xi_12 = (C_1)/(C_1 + C_2);
  parameter Real Xi_21 = (C_2)/(C_2 + C_1);
  parameter Real Xi_23 = (C_2)/(C_2 + C_3);
  parameter Real Xi_32 = (C_3)/(C_3 + C_2);
  parameter Real Xi_34 = (C_3)/(C_3 + C_4);
  parameter Real Xi_43 = (C_4)/(C_4 + C_3);
  parameter Real eta_1 = (C_1)/(C_1 + C_5 + C_6);
  parameter Real eta_4 = (C_4)/(C_1 + C_5 + C_6);
  parameter Real eta_5 = (C_5)/(C_1 + C_5 + C_6);
  parameter Real eta_6 = (C_6)/(C_1 + C_5 + C_6);
  parameter Real eta_prime_4 = (C_4)/(C_4 + C_5 + C_6 + h_cond*A_cond_outer/2.);
  parameter Real eta_prime_5 = (C_5)/(C_4 + C_5 + C_6 + h_cond*A_cond_outer/2.);
  parameter Real eta_prime_6 = (C_6)/(C_4 + C_5 + C_6 + h_cond*A_cond_outer/2.);
   
  //----------------------------
  // State variables (temperatures)
  //----------------------------
  Real T1(start=1073.15) "Evaporator wall (K) ~800 C";
  Real T2(start=1073.15) "Evaporator wick (K)";  
  Real T3(start=873.15)  "Condenser wick (K) ~600 C";
  Real T4(start=873.15)  "Condenser wall (K)";
  Real T5(start=973.15) "Adiabatic wick (K)";
  Real T6(start=973.15) "Adiabatic wall (K)";

  //----------------------------
  // Local conductance helper: radial conduction through wall (cylindrical shell)  
  //----------------------------
  function k_radial_shell
    input Real k_mat "material conductivity";
    input Real r_in;
    input Real r_out;
    input Real Lax "axial length (m)";
    output Real G "conductance W/K";
  algorithm
    // using conduction of cylindrical shell between r_in and r_out along axial length:
    // R_th = (ln(r_out/r_in))/(2*pi*k_mat*Lax)
    // G = 1/R_th
    G := 1.0 / ( (log(r_out/r_in))/(2*pi*k_mat*Lax) );
  end k_radial_shell;

  //----------------------------
  // Radial conductances used in equations (explicit)
  //----------------------------
  Real G1_2 = k_radial_shell(k_wall, ro-twall, ro, Le);
  Real G5_6 = k_radial_shell(k_wall, ro-twall, ro, Lc);

  // axial conductance through wick or wall between sections (approx k*A/L)
  Real G2_3 = k_wick * A_wick_radial / Le;
  Real G3_4 = k_wall * (2*pi*(ri+wick_thickness+twall/2)*La) / La; // approximate axial through wall area / La -> reduces to k*2pi*r
  Real G4_5 = k_wick * A_wick_radial / Lc;

equation
/*
  // Eqn (4) Evaporator wall energy balance
  C1*der(T1) = G1_2*(T2 - T1) + Qe;

  // Eqn (5) Evaporator wick
  C2*der(T2) = G1_2*(T1 - T2) + G2_3*(T3 - T2);

  // Eqn (6) Adiabatic wall
  C3*der(T3) = G2_3*(T2 - T3) + G3_4*(T4 - T3);

  // Eqn (7) Adiabatic wick
  C4*der(T4) = G3_4*(T3 - T4) + G4_5*(T5 - T4);

  // Eqn (8) Condenser wick
  C5*der(T5) = G4_5*(T4 - T5) + G5_6*(T6 - T5);

  // Eqn (9) Condenser wall w/ convective sink to T_inf
  C6*der(T6) = G5_6*(T5 - T6) - h_cond * A_cond_outer * (T6 - T_inf);
*/


  der(T1) = 2*alpha_1/L_1^2*( (Xi_12 + eta_1 - 2)*T1 + Xi_21*T2 + eta_5*T5 + eta_6*T6 + (Qe/2.)/(C_1 + C_5 + C_6) );
  
  der(T2) = 2*alpha_2/L_2^2*( Xi_12*T1 + (Xi_21 + Xi_23 - 2)*T2 + Xi_32*T3 );
  
  der(T3) = 2*alpha_3/L_3^2*( Xi_23*T2 + (Xi_32 + Xi_34 - 2)*T3 + Xi_43*T4 );
  
  der(T4) = 2*alpha_4/L_4^2*(Xi_34*T3 + (Xi_43 + eta_4 - 2)*T4 + eta_prime_5*T5 + eta_prime_6*T6 + (h_cond*A_cond_outer/2.*T_inf)/(C_4 + C_5 + C_6 + h_cond*A_cond_outer/2.) );
  
  der(T5) = 2*alpha_5/L_5^2*(eta_1*T1 + eta_prime_4*T4 + (eta_5 + eta_prime_5 - 2)*T5 + (eta_6 + eta_prime_6)*T6 + (Qe/2.)/(C_1 + C_5 + C_6) + (h_cond*A_cond_outer/2.*T_inf)/(C_4 + C_5 + C_6 + h_cond*A_cond_outer/2.) );
  
  der(T6) = 2*alpha_6/L_6^2*(eta_1*T1 + eta_prime_4*T4 + (eta_5 + eta_prime_5)*T5 + (eta_6 + eta_prime_6 - 2)*T6 + (Qe/2.)/(C_1 + C_5 + C_6) + (h_cond*A_cond_outer/2.*T_inf)/(C_4 + C_5 + C_6 + h_cond*A_cond_outer/2.) );
  
  
  
  /*
  der(T1) = 2*alpha_1/L_1^2*( (Xi_12 + eta_1 - 2)*T1 + Xi_21*T2) + 2*alpha_1/L_1^2*(Qe/2.)/(C_1);
  
  der(T2) = 2*alpha_2/L_2^2*(Xi_12*T1 + (Xi_21 + Xi_23 - 2)*T2 + Xi_32*T3);
  
  der(T3) = 2*alpha_3/L_3^2*(Xi_23*T2 + (Xi_32 + Xi_34 - 2)*T3 + Xi_43*T4);
  
  der(T4) = 2*alpha_4/L_4^2*(Xi_34*T3 + (Xi_43 + eta_4 - 2)*T4) + 2*alpha_4/L_4^2*(h_cond*A_cond_outer/2.*T_inf)/(C_4 + h_cond*A_cond_outer/2.);
  */
  
//annotation (uses(Modelica));



end HeatPipe_ZuoFaghri;
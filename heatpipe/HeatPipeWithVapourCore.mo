model HeatPipeWithVapourCore
/*
This model, inspired by Zhang et al. ("Comparative study of two quick-analysis models for frozen startup of high-temperature heat pipes"),
uses a dynamic thermal resistance to represent the vapour core of the KRUSTY sodium heat pipe. 
This resistance mimics the effect of the heat pipe "activating" during start-up, which is a function of the vapour pressure and temperature. This function is highly non-linear due to the sharp rarefied->continuous flow transition occuring at the transition temeprature. This temperature
is computed from the critical Knudsen number based on kinetic theory. 

The effect of this is that, during cold/frozen start-up, there is initially almost no heat transfer to the condenser.
The rest of the heat pipe should behave the same way as the Zuo&Faghri-inspired model.

Note that the melting phase of start-up is neglected as the total enthalpy of fusion is negligible, as also shown here. 

Inputs: Number of heat pipes being modelled together, heat flux into evaporator, condenser/Stirling boundary condition
Outputs: Evaporator wall temperature, 
*/
  import Modelica.Constants.pi;

  //----------------------------
  // Geometry (KRUSTY)
  //----------------------------
  parameter Real Le = 0.25 "evaporator length [m]"; 
  parameter Real Lc = 0.10 "condenser length [m]"; 
  parameter Real La = 0.65 "adiabatic length [m]";

  parameter Real do_wall = 0.0127 "wall OD [m]"; 
  parameter Real ro_wall = do_wall/2;
  parameter Real ri_wall = ro_wall - twall;
  parameter Real ro_wick = ri_wall;
  parameter Real ri_wick = ro_wick - twick;
  parameter Real di_wick = 2*ri_wick;
  
  parameter Real twall = 0.00089 "wall thickness [m]";  
  parameter Real twick = 1e-3 "wick thickness, estimated [m]";
  
  parameter Real A_v_wicked = pi*ri_wick^2. "cross-sectional area of the vapour core in the wicked region (evaporator)";
  parameter Real A_v_unwicked = pi*ri_wall^2. "cross-sectional area of the vapour core in the unwicked regions (adiabatic, condenser)";
  
  parameter Real r_vapour_eff = sqrt(V_vcore/(Le + La + Lc)/pi); 
  //----------------------------
  // Material properties 
  //----------------------------
  
  // Wall
  parameter Real k_wall = 24 "wall conductivity [W/m/k]";
  parameter Real rho_wall = 9050 "wall density [kg/m3]";
  parameter Real cp_wall = 595 "wall specific heat [J/kg/K]";

  // Wick (screen + liquid)
  parameter Real rho_screen = 8000 "kg/m3 (approx stainless screen)";
  parameter Real cp_screen = 500 "J/kg.K (screen)";
  parameter Real k_screen = 20 "W/m.K (screen)";
  parameter Real rho_Na_l = rho_l_correlation(T_melt); //fixed at the melting point
  parameter Real cp_Na_l = cp_l_correlation(T_melt);
  parameter Real k_Na_l = 70 "conductivity of sodium liquid [W/m/K]";  
  parameter Real porosity = 0.7;  
  
  parameter Real m_wick_tot = (rho_screen*V_wick*(1-porosity) + rho_Na_l*V_wick*porosity);
  parameter Real rho_wick = m_wick_tot/V_wick;    
  parameter Real cp_wick = ( (1-porosity)*rho_screen*cp_screen + porosity*rho_Na_l*cp_Na_l)/rho_wick; //Chi model (mass-weighted eff. cp)
  parameter Real k_wick =
    k_Na_l*((k_Na_l + k_screen) - (1 - porosity)*(k_Na_l - k_screen)) /
    ((k_Na_l + k_screen) + (1 - porosity)*(k_Na_l - k_screen)) "wick effective conductivity, Chi model [W/m/K]";

  // Vapour
  Real cp_v "sodium vapour heat capacity at constant pressure";

  parameter Real d_v  = 4.54E-10 "effective molecular diameter of sodium vapour, using Van der Waals radius [m]";
  parameter Real kappa = 1.380649E-23 "Boltzmann constant [J/K]";
  
 
  parameter Real M_g = 0.02299 "molar mass of monatomic sodium [kg/mol]";
  parameter Real R_g = 8.31446/M_g "specific gas constant [J/kg/K]";
   
  parameter Real h_lv = 4237E3 "latent heat of vapourisation [J/kg]";
  
  // Liquid (testing)
  parameter Real h_sl = 113E3 "latent heat of fusion [J/kg]";
  
  parameter Real T_melt = 370.97 "sodium melting point [K]";
  parameter Real cp_l_thresh = cp_l_correlation(T_melt) "liquid sodium heat capacity at melting point";
  parameter Real rho_l_thresh =  rho_l_correlation(T_melt) "liquid sodium density at melting point";
  parameter Real H_thresh_melt = rho_l_thresh*cp_l_thresh*T_melt*V_liquid "internal enthalpy (energy) of sodium at the melting point [J]"; 
  
  parameter Real m_melt_total = rho_l_thresh*V_liquid;
  parameter Real H_melt_total = h_sl*m_melt_total;
  parameter Real eq_temp_rise = H_melt_total/cp_l_thresh;
  
  
  //----------------------------
  // Volumes & capacitances
  //----------------------------
  parameter Real V_wall_e = pi*(ro_wall^2 - ri_wall^2)*Le;
  parameter Real V_wall_c = pi*(ro_wall^2 - ri_wall^2)*Lc;
  parameter Real V_wall_a = pi*(ro_wall^2 - ri_wall^2)*La;
  parameter Real V_wick = pi*(ro_wick^2 - ri_wick^2)*(Le); //KRUSTY heat pipes were only wicked in the evaporator (thermosiphon)
  parameter Real V_liquid = V_wick*porosity; //crude estimate, assume the wick is exactly saturated with sodium
  parameter Real V_vcore = A_v_wicked*Le + A_v_unwicked*(La + Lc) "volume of the vapour core [m^3]"; 

  parameter Real C_evap_wall = rho_wall * V_wall_e * cp_wall "evaporator wall capacitance [J/K]";
  //parameter Real C_evap_wall = 0.1*(rho_wall * V_wall_e * cp_wall) "evaporator wall capacitance [J/K]";
  parameter Real C_cond_wall = rho_wall * V_wall_c * cp_wall "condenser wall capacitance [J/K]";
  parameter Real C_wick = rho_wick * cp_wick * V_wick "wick effective capacitance [J/K]";
  parameter Real C_adiab_wall = rho_wall * V_wall_a * cp_wall "adiabatic wall capacitance [J/K]";
  //parameter Real C_adiab_wall = 0.0 "adiabatic wall capacitance [J/K]";
  
  Real C_vapour "vapour capacitance, a function of instantaneous density/concentration [J/K]"; 
    
  parameter Real C_Stirling = 3*C_cond_wall "Estimated/tuned Stirling capacitance [J/K]"; 

  //----------------------------
  // Resistances  
  //----------------------------
  parameter Real R_wall_e = 1/(2*pi*k_wall*Le)*log(ro_wall/ri_wall); 
  parameter Real R_wick_e = 1/(2*pi*k_wick*Le)*log(ro_wick/ri_wick); 

  parameter Real R_wall_c = 1/(2*pi*k_wall*Lc)*log(ro_wall/ri_wall); 
  //There is no wick in the condenser or adiabitic region in KRUSTY
  parameter Real R_wick_c = 0.0;   //1/(2*pi*k_wick*Lc)*log(ro_wick/ri_wick);  
  //There is no wick in the condenser or adiabitic region in KRUSTY
  parameter Real R_wick_a = 0.0;   //1/(2*pi*k_wick*Lc)*log(ro_wick/ri_wick);
  parameter Real R_evap_radial = R_wall_e + R_wick_e;
  parameter Real R_cond_radial = R_wall_c + R_wick_c;
  
  parameter Real R_wall_axial = (La + 0.5*Le + 0.5*Lc )/(k_wall*pi*(ro_wall^2. - ri_wall^2)) "effective axial thermal resistance of the wall"; 
  //parameter Real R_wall_axial = 100*(La + 0.5*Le + 0.5*Lc )/(k_wall*pi*(ro_wall^2. - ri_wall^2)) "effective axial thermal resistance of the wall"; 
  
  parameter Real R_wv = 1E-6 "wick-vapour interface resistance";
  
  Real R_vapour_ax "Axial resistance of the vapour core, dependent on vapour pressure and flow transition"; 
  
  parameter Real Q_cond_nominal = 2350/8;
  parameter Real R_stirling_hp_interface = 165/Q_cond_nominal;
  parameter Real T_stirling_nominal = 630 + 273.15 "Stirling hot-side temperature, Poston et al., Fig. 10 [K]"; 
  parameter Real T_cond_nominal = T_stirling_nominal + Q_cond_nominal*R_stirling_hp_interface; 
  parameter Real HTC = Q_cond_nominal/(T_cond_nominal - T_stirling_nominal);
  
  parameter Real T_stirling_cold_nominal  = 65 + 273.15; 
  parameter Real HTC_cold = Q_cond_nominal/(T_stirling_nominal - T_stirling_cold_nominal);
  
  //----------------------------
  // Helper correlation functions
  //----------------------------
  
  /*vapour pressure correlation*/
  function pv_correlation
    input Real T "Temperature [K]";    
    output Real pv "vapor pressure [Pa]";  
    protected constant Real h_lv = 4237.0E3;
    protected constant Real R_g = 8.31446/0.02299; 
  algorithm
    //Theoretical (Clausius-Clapeyron), relating to the lowest experimentally known vapour pressure (at 0.49 atm) 
    pv := 49.64925e3*exp(h_lv/R_g*(1/1072. - 1/T));
    
    //Extrapolating from non-vacuum vapour pressures    
    //pv := 10^(4.51961 - 5202.12/T) * 101325.; //source: https://ntrs.nasa.gov/api/citations/19650014783/downloads/19650014783.pdf
    //pv := 10^( 2.46077 - 1873.728/(T - 416.372) ) * 100e3; //source: NIST https://webbook.nist.gov/cgi/inchi?ID=C7440235&Mask=4
  end pv_correlation;
  
  /*heat capacity correlation. source: https://webbook.nist.gov/cgi/cbook.cgi?ID=C7440235&Mask=1E9F#:~:text=Table_title:%20Gas%20Phase%20Heat%20Capacity%20(Shomate%20Equation),Data%20last%20reviewed%20in%20June%2C%201962%20%7C */
  function cp_v_correlation
    input Real T "Temperature [K]";
    output Real cp "heat capacity at constant pressure [J/kg/K]";
  algorithm 
    cp := ( 20.80573 + 0.277206*(T/1000.) - 0.392086*(T/1000.)^2. + 0.119634*(T/1000.)^3. - 0.008879/(T/1000.)^2 ) / 0.023;
  end cp_v_correlation;  
  
  function mu_v_correlation 
    input Real T "Temperature [K]";
    output Real mu "Dynamic viscosity of sodium vapour [Pa*s]";
  algorithm 
    mu := 8.4491E-9*(T-273.15) + 4.8898E-6; //correlation derived from table by Dunning https://www.osti.gov/servlets/purl/4120472
    //mu := 2.872E-5 //Source: https://scispace.com/pdf/thermophysical-properties-of-sodium-2v306hkn25.pdf";
  end mu_v_correlation;
  
  function cp_l_correlation
    input Real T "Temperature [K]";
    output Real cp "heat capacity at constant pressure [J/kg/K]";
  algorithm 
    cp := ( 40.25707 - 28.23849*(T/1000) + 20.69402*(T/1000)^2. - 3.641872*(T/1000)^3. - 0.079874/(T/1000)^2. ) / 0.023;
  end cp_l_correlation;      
  
  function rho_l_correlation
    input Real T "Temperature [K]";
    output Real rho "density [kg/m3]";
  algorithm 
    rho := (0.927 - 0.238E-3*(T-273.15-100))*1000.; 
  end rho_l_correlation;  
  
  
  //----------------------------
  // States
  //----------------------------
  Real T_monolith (start=300, fixed=false "fuel monolith avg temp [K]"); 
  Real T_evap (start=300, fixed = true) "evaporator wall temp [K]";
  Real T_wick (start=300, fixed = true) "wick temp (evap) [K]";
  Real T_vap (start=300, fixed = true) "vapour temperature [K]";  
  Real T_cond (start=300, fixed = true) "condenser wall temp [K]";  
  
  Real T_Stirling (start=300, fixed = false) "Stirling PCS lumped temperature [K]"; //may need to switch 'fixed' true/false  
  
  Real p_v  "vapour pressure at instantaneous temperature [Pa]";
  Real rho_v "instantaneous vapour density [kg/m3]";  
  Real mu_v "instantaneous dynamic viscosity [Pa*s]";
  Real DeltaP "pressure drop in vapour core"; 
  
  //----------------------------
  // Toy model of the core -- assume single lump, i.e. no internal temperature gradient
  // This model is only meant to be used for quick tests without external FMU coupling.
  //----------------------------  
  parameter Real rho_monolith = 17110;
  parameter Real cp_monolith = 189;
  parameter Real V_monolith = 0.00186454 "monolith volume, from mesh [m^3]"; 
  parameter Real C_monolith = rho_monolith*cp_monolith*V_monolith;    
  parameter Real w_contact = 0.015;
  parameter Real A_contact = w_contact*Le;
  parameter Real L_mono_evap = 0.01; //take ~half the monolith thickness as an approx.
  parameter Real k_monolith = 37.5; 
  parameter Real R_mono_evap = L_mono_evap/A_contact/k_monolith;
  
  //----------------------------
  // Flow transition 
  //----------------------------
  parameter Real Knuds_crit = 0.01 "critical Knudsen number [-] for transition to continuum flow";  
  Real T_tr (start=1000., fixed=false) "Transition temperature (rarefied->continuum flow) [K]";  
  Real p_v_tr "vapour pressure at the computed transition temperature";
  Real rho_v_tr "vapour density at the computed transition temperature";    
  Real mu_v_tr "dynamic viscosity at the computed transition temperature";
  
  Real R_vapour_ax_tr "vapour resistance at the computed transition temperature";
  Real R_vapour_ax_posttransition; 
  
  
  //Smoothing function 
  Real S; 
  parameter Real beta = 1; //tunable sigmoid parameter: higher value gives sharper transition which means the "stall" in evaporator temperature happens closer to the transition temperature
  parameter Real R_high = 1e6; //arbitrary high resistance value 
  
  //----------------------------
  // Heat flows
  //----------------------------
  Real Q_ew "evaporator->wick";
  Real Q_wv "wick->vapour";
  Real Q_vc "vapour->condenser";  
  Real Q_ec "evaporator->condenser";
  Real Q_cs "condenser->Stirling CPS";
  Real R_vc "combined thermal resistance from vapour to condenser"; 
  
  //----------------------------
  // Inputs
  //----------------------------
  input Integer N_HPs_input;  
  input Real Q_evap_input;
  input Real Q_cond_input;
  input Real Q_stirling_input;
  Real Q_evap_bc "heat addition rate to evaporator [W]";
  Real Q_cond_bc "heat removal rate from condenser [W]";
  Real Q_Stirling_bc "heat removal rate from Stirling PCS [W]";
  Real N_HPs;
  
  
  //----------------------------
  // Model options -- these are not intended to be selectable by the FMU handler
  //----------------------------
  //Option to account for heat pipe activation or not 
  parameter Boolean MODEL_HP_STARTUP = true;
  //Option to model the core as a single lump or simply take Q_evap as a boundary condition
  parameter Boolean MODEL_CORE_INTERNALLY = false; 
  //Give the option to model the thermal mass of the Stirling PCS
  parameter Boolean MODEL_STIRLING_MASS = true; 
    
  //----------------------------
  // Latches (flags with persistence)
  //----------------------------
  Real STIRLING_ACTIVATED (start=0, fixed=true);
    
  //----------------------------
  // Outputs / signals
  //----------------------------
  Boolean HP_activated "Flag raised when the heat pipe has effective heat transfer";
  
  //track the Reynolds number of the vapour throughout the simulation
  Real v_v "vapour velocity [m/s]" ;
  Real mdot_v "vapour mass flow rate [kg/s]";
  Real Re; 

equation
  cp_v = cp_v_correlation(T_vap);
  
  //The factor 1.051/sqrt(2) comes from the corrected mean free path accounting for the motion of other particles during the flight of a an average particle
  T_tr = ( sqrt(2)*pi*d_v^2.*pv_correlation(T_tr)*(di_wick) )/( 1.051*kappa )*Knuds_crit; //non-linear equation solved automatically with initial guess. 
  
  //Note: The experimentally-derived heat pipe activation temperature in KRUSTY was ~773 K as opposed to ~700 K computed here.

  
  p_v = pv_correlation(T_vap);
  p_v_tr = pv_correlation(T_tr);
  rho_v = p_v/R_g/T_vap;
  rho_v_tr = p_v_tr/R_g/T_tr;
  mu_v = mu_v_correlation(T_vap);
  mu_v_tr = mu_v_correlation(T_tr); 
  DeltaP = 8*mu_v_correlation(T_vap)*(La+Le/2+Lc/2)/(h_lv*rho_v*r_vapour_eff^4.)*Q_wv;// just to inspect the pressure drop 
  
  
  /*axial convective thermal resistance derived from Clausius-Clapeyron (Guo et al.) */
  R_vapour_ax_tr = 8*(R_g)*T_tr^2.*mu_v_tr*(La+Le/2.+Lc/2.)/(h_lv^2.*p_v_tr*pi*rho_v_tr*r_vapour_eff^4.);
  R_vapour_ax_posttransition = 8*(R_g)*T_vap^2.*mu_v*(La+Le/2.+Lc/2.)/(h_lv^2.*p_v*pi*rho_v*r_vapour_eff^4.);  
  //alternative formulation based on integral Claus-Clapeyron--numerically unstable
  //R_vapour_ax_posttransition = (1/(R_g/h_lv*log(p_v/(p_v - DeltaP)) + 1/T_vap) - T_vap)/Q_wv;     
  
  
  S = 1 / (1 + exp(-beta*(T_vap - T_tr))); //Smoothing sigmoid function   
  if MODEL_HP_STARTUP then 
    R_vapour_ax = (1 - S)*R_high + S*(R_vapour_ax_posttransition); 
    //R_vapour_ax = (1 - S)*R_high + S*(R_vapour_ax_tr); //the resistance at the transition temp is low enough that it may as well be used throughout.
  else 
    R_vapour_ax = 1e-6; //stick to some arbitrarily low resistance (short-circuit) at all times, reverting to Zuo-Faghri type model 
  end if;
  
  mdot_v = (Q_wv/h_lv);
  v_v = mdot_v/rho_v/(pi*r_vapour_eff^2);
  Re = rho_v*v_v*(2*r_vapour_eff)/mu_v; 
  
  /* Basic test for heat pipe activation: is much more heat conducted through the vapour than through the wall? */  
  if R_vapour_ax < R_wall_axial/10. then 
    HP_activated = true;
  else 
    HP_activated = false;
  end if;
    
  
  //----------------------------
  // Heat transfer relations
  //----------------------------
  
  if N_HPs_input > 0 then
    N_HPs = N_HPs_input;
  else 
    N_HPs = 8;
  end if; 
  
  
  if MODEL_CORE_INTERNALLY then 
    Q_evap_bc = (T_monolith - T_evap)/R_mono_evap;     
  else
    /*
    if abs(Q_evap_input) > 1E-9 then //need the absolute value in case there is an initially higher temp in the heat pipes than in the core
    //if Q_evap_input <> 0 then 
      Q_evap_bc = Q_evap_input/N_HPs;
    else 
      Q_evap_bc = Q_cond_nominal;
    end if;
    */
    Q_evap_bc = Q_evap_input/N_HPs;  
  end if;
  
  if MODEL_STIRLING_MASS then 
    Q_cond_bc = 0;
    if Q_stirling_input > 1E-9 then
      Q_Stirling_bc = Q_stirling_input/N_HPs;
    else
      Q_Stirling_bc = HTC_cold*(T_Stirling - T_stirling_cold_nominal); 
      //Q_Stirling_bc = Q_cond_nominal;
      //Q_Stirling_bc = (Q_cond_nominal*8-640)/8; //trigger a load rejection after a long time (t=8h transient: 640 W of load rejection in total)          
    end if;
  else 
    if abs(Q_cond_input) > 1E-9 then 
      Q_cond_bc = Q_cond_input/N_HPs;
    else       
      Q_cond_bc = max(0, HTC*(T_cond - T_stirling_nominal));  //Completely avoid adding heat to the condenser based on this HTC
    end if;
    Q_Stirling_bc = 0;
  end if;  
  
  /*Heat transfer rates*/
  
  // Evaporator->wick
  Q_ew  = (T_evap - T_wick) / R_evap_radial;

  // Wick->vapor 
  Q_wv = (T_wick - T_vap) / R_wv;  

  // Vapor->condenser 
  R_vc = max(R_cond_radial, R_cond_radial + R_vapour_ax);   
  Q_vc = (T_vap - T_cond) / R_vc;  
  
  // Evaporator->condenser
  Q_ec = (T_evap - T_cond) / R_wall_axial;
  
  // Condenser->Stirling 
  Q_cs = (T_cond - T_Stirling) / R_stirling_hp_interface; //thermal resistance according to prescribed nominal DeltaT

  //----------------------------
  // Energy balances
  //----------------------------
  
  // Core monolith - based on the full geometry. 
  if MODEL_CORE_INTERNALLY then 
    C_monolith * der(T_monolith) = 2350 - 8*Q_evap_bc;
    //C_monolith * der(T_monolith) = 2350 - 2*Q_evap_bc; //only two Stirling convertors are active at first
  else 
    T_monolith = 0; //placeholder value
  end if; 
  
  // Evaporator 
  //(C_evap_wall + C_adiab_wall) * der(T_evap) = Q_evap_bc - Q_ew - Q_ec; //with thermal mass of adiabatic section lumped in  
  C_evap_wall * der(T_evap) = Q_evap_bc - Q_ew - Q_ec;

  // Wick
  C_wick * der(T_wick) = Q_ew - Q_wv; //Neglecting axial conductance in the wick

  // Vapour
  C_vapour = rho_v*cp_v*V_vcore; //Update vapour capacitance based on density (rarefied/continuous state)
  C_vapour * der(T_vap) = Q_wv - Q_vc;
  
  // Condenser and Stirling PCS
  if MODEL_STIRLING_MASS then    
    //C_cond_wall * der(T_cond) = Q_vc + Q_ec - Q_cs;    
    (C_cond_wall + C_adiab_wall) * der(T_cond) = Q_vc + Q_ec - Q_cs;    
    
    //latch trigger
    if (T_Stirling > T_stirling_nominal) then       
      der(STIRLING_ACTIVATED) = 1; //positive derivative that gets integrated to activate the latch
    else       
      der(STIRLING_ACTIVATED) = 0;
    end if;
    
    if (STIRLING_ACTIVATED > 0) then 
      C_Stirling * der(T_Stirling) = Q_cs - Q_Stirling_bc;       
    else 
      C_Stirling * der(T_Stirling) = Q_cs; //The Stirling engine can only heat up until it's turned on
    end if;
  else 
    //C_cond_wall * der(T_cond) = Q_vc + Q_ec - Q_cond_bc;
    (C_cond_wall + C_adiab_wall) * der(T_cond) = Q_vc + Q_ec - Q_cond_bc;
    T_Stirling = 0; //placeholder value
    
    der(STIRLING_ACTIVATED) = 0; //Don't need the latch in this simplified case. 
  end if; 
  
  

end HeatPipeWithVapourCore;
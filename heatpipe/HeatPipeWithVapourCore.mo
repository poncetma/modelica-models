model HeatPipeWithVapourCore
/*
This model, inspired by Zhang et al. ("Comparative study of two quick-analysis models for frozen startup of high-temperature heat pipes"),
uses a dynamic thermal resistance to represent the vapour core of the KRUSTY sodium heat pipe. 
This resistance mimics the effect of the heat pipe "activating" during start-up, which is a function of the vapour pressure and temperature. This function is highly non-linear due to the sharp rarefied->continuous flow transition occuring at the transition temeprature. This temperature
is computed from the critical Knudsen number based on kinetic theory. The effect of this is that, during cold/frozen start-up, there is initially almost no heat transfer to the condenser. The rest of the heat pipe should behave the same way as the Zuo&Faghri-inspired model. Note that the melting phase of start-up is neglected as the total enthalpy of fusion is negligible, as also shown here. 

Throughput limits from KRUSTY's heat pipes (viscous limit, flooding limit) have also been added, derived from a high-fidelity/advective model (Poston et al., "Results of the KRUSTY Nuclear System Test). 

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
  //parameter Real k_wall = 24 "wall conductivity [W/m/k]";
  //parameter Real rho_wall = 9050 "wall density [kg/m3]";
  parameter Real rho_wall = 8970 "wall density [kg/m3] - official HAYNES 230 data";
  //parameter Real cp_wall = 595 "wall specific heat [J/kg/K]";

  // Wick (screen + liquid)
  parameter Real rho_screen = 8000 "kg/m3 (approx stainless screen)";
  parameter Real cp_screen = 601.5 "J/kg.K (screen)"; //SS316 @800C (Kim (ANL), tech rep)
  parameter Real k_screen = 26.1 "W/m.K (screen)"; //SS316 @800C (Kim (ANL), tech rep)
  parameter Real rho_Na_l = rho_l_correlation(T_melt); //fixed at the melting point
  parameter Real cp_Na_l = cp_l_correlation(T_melt);
  parameter Real k_Na_l = 58 "conductivity of sodium liquid [W/m/K]"; //@800C IAEA data collection  
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

  Real C_evap_wall;// = rho_wall * V_wall_e * cp_wall "evaporator wall capacitance [J/K]";
  Real C_cond_wall;// = rho_wall * V_wall_c * cp_wall "condenser wall capacitance [J/K]";
  parameter Real C_wick = rho_wick * cp_wick * V_wick "wick effective capacitance [J/K]";
  Real C_adiab_wall;// = rho_wall * V_wall_a * cp_wall "adiabatic wall capacitance [J/K]";
  
  Real C_vapour "vapour capacitance, a function of instantaneous density/concentration [J/K]"; 
    
  //The Stirling capacitance is tuned based on Fig. 10 of Poston et al. It should have the same initial temperature drop once the heat rmeoval is activated.
  parameter Real C_stirling = 17.780965 "Stirling capacitance [J/K]"; //tuned multiple based on start-up datag ;// = 15*C_cond_wall ; 

  //----------------------------
  // Resistances  
  //----------------------------
  parameter Real R_wall_e = 1/(2*pi*k_wall(800)*Le)*log(ro_wall/ri_wall); 
  parameter Real R_wick_e = 1/(2*pi*k_wick*Le)*log(ro_wick/ri_wick); 

  parameter Real R_wall_c = 1/(2*pi*k_wall(800)*Lc)*log(ro_wall/ri_wall); 
  parameter Real R_wall_a = 1/(2*pi*k_wall(800)*La)*log(ro_wall/ri_wall); 
  //There is no wick in the condenser or adiabitic region in KRUSTY
  parameter Real R_wick_c = 0.0;   //1/(2*pi*k_wick*Lc)*log(ro_wick/ri_wick);  
  //There is no wick in the condenser or adiabitic region in KRUSTY
  parameter Real R_wick_a = 0.0;   //1/(2*pi*k_wick*Lc)*log(ro_wick/ri_wick);
  parameter Real R_evap_radial = R_wall_e + R_wick_e;
  parameter Real R_cond_radial = R_wall_c + R_wick_c;
  parameter Real R_adiab_radial = R_wall_a + R_wick_a; 
  
  //parameter Real R_wall_axial = (La + 0.5*Le + 0.5*Lc )/(k_wall(800)*pi*(ro_wall^2. - ri_wall^2)) "effective axial thermal resistance of the wall"; 
  parameter Real R_wall_axial_ea = (0.5*Le + 0.5*La)/(k_wall(800)*pi*(ro_wall^2. - ri_wall^2.));
  parameter Real R_wall_axial_ac = (0.5*La + 0.5*Lc)/(k_wall(800)*pi*(ro_wall^2. - ri_wall^2.));
  
  parameter Real R_wv = 1E-4 "wick-vapour interface resistance"; //Setting this too low leads to instabilities (but it should be << 1e-3) //1E-6
  
  Real R_vapour_ax "Axial resistance of the vapour core, dependent on vapour pressure and flow transition"; 
  
  
  parameter Real T_stirling_activation = 650 + 273.15 "Temperature at which Stirlings are turned on in start-up run of KRUSTY"; 
  parameter Real T_stirling_nominal = 650 + 273.15 "Stirling hot-side temperature, Poston et al., Fig. 10 [K]"; //630 + 273.15
  parameter Real T_stirling_cold_nominal  = 65 + 273.15; 
  
  /*The below parameters don't need to be computed based on the number of heat pipes as we can assume each heat pipe has its own nominal power draw equal to 2250/8 */
  parameter Real Q_draw_nominal = 2250/8 "nominal power draw [W]";
  parameter Real R_stirling_hp_interface = ( (T_stirling_nominal + 145) - T_stirling_nominal )/Q_draw_nominal;     
  parameter Real T_cond_nominal = T_stirling_nominal + Q_draw_nominal*R_stirling_hp_interface;   
  parameter Real HTC = Q_draw_nominal/(T_cond_nominal - T_stirling_nominal);  
  parameter Real HTC_cold = Q_draw_nominal/(T_stirling_nominal - T_stirling_cold_nominal);
  
  //----------------------------
  // Helper correlation functions
  //----------------------------
  
  /*Haynes 230 heat capacity correlation*/
  function cp_wall
    input Real T "Temperature [K]";
    output Real cp "heat capacity [J/kg/K]";
    protected Real T_c = T - 273.15; 
  algorithm 
    cp := 0.2403*T_c + 381; 
  end cp_wall;
  
  /*Haynes 230 conductivity correlation*/
  function k_wall
    input Real T "Temperature [K]";
    output Real k "heat capacity [J/kg/K]";
    protected Real T_c = T - 273.15; 
  algorithm
    k := 1.998E-2*T_c + 8.413; 
  end k_wall;
  
  /*vapour pressure correlation*/
  function pv_correlation
    input Real T "Temperature [K]";    
    output Real pv "vapour pressure [Pa]";  
    protected constant Real h_lv = 4237.0E3;
    protected constant Real R_g = 8.31446/0.02299; 
  algorithm
    //Theoretical (Clausius-Clapeyron), relating to the lowest experimentally known vapour pressure (at 0.49 atm) 
    pv := 49.64925e3*exp(h_lv/R_g*(1/1072. - 1/max(T,1))); //made numerically safe
    
    //Empirical correlations
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
  
  /*************/
  /*Power throughput limit correlations (flooding limit, viscous limit)*/
  /*************/
  //Data from Poston et al. (high-fidelity model), Fig. 7
  function viscous_limit
    input Real T "Sodium temperature [K]"; 
    output Real Q_max "Power limit [W]";
    protected Real T_c;
  algorithm
    T_c := T - 273.15; 
    if T_c > 400. then //Restrict to range of validity
      //Q_max := 1.068E-07*exp(0.03991*T_c); //own correlation from digitised data (R^2=0.993)
      Q_max := 8.441E-53*T_c^(1.994E+1); //own correlation from digitised data (R^2 = 0.9990)
    elseif T_c > 0 then 
      //Q_max := 1.068E-07*exp(0.03991*400.); //Evaluate at lowest bound
      //Q_max := 8.441E-53*400^(1.994E+1); //Evaluate at lowest bound
      Q_max := 8.441E-53*T_c^(1.994E+1); //Extroplate
    else
      Q_max := 0.;
    end if; 
    
  end viscous_limit;
  
  function flooding_limit
    input Real T "Sodium temperature [K]"; 
    output Real Q_max "Power limit [W]";
    protected Real T_c;
  algorithm
    T_c := T - 273.15; 
    if T_c > 400.0 then 
      Q_max := 3.404E-3*T_c^2. -2.202*T_c + 3.800E+2; //own correlation from digitised data (R^2 = 0.9999)
    else
      //Q_max := 3.404E-3*400.^2. -2.202*400. + 3.800E+2; //Evaluated at lowest bound
      Q_max := 3.404E-3*T_c^2. -2.202*T_c + 3.800E+2; //Extrapolate
    end if; 
  end flooding_limit;
  
  //----------------------------
  // States
  //----------------------------
  Real T_monolith (start=1081, fixed=false) "fuel monolith avg temp [K]"; 
  Real T_evap  (start=1073.15, fixed=false)"evaporator wall temp [K]";
  Real T_wick  (start=1073.15, fixed=false)"wick temp (evap) [K]";
  Real T_vap (start=1073.15, fixed=false)  "vapour temperature [K]";  
  Real T_cond (start=1073.15, fixed=false) "condenser wall temp [K]";  
  
  Real T_adiab (start = 1073.15, fixed = false) "adiabatic wall temp [K]"; 
  
  Real T_stirling (start=T_stirling_nominal, fixed = false) "Stirling PCS lumped temperature [K]"; //may need to switch 'fixed' true/false  
  
  Real p_v (start=1) "vapour pressure at instantaneous temperature [Pa]";
  Real rho_v (start=1) "instantaneous vapour density [kg/m3]";  
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
  Real R_vapour_ax_posttransition (start=1); 
  
  
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
  //Real Q_ec "evaporator->condenser";
  
  Real Q_va "vapour->adiabatic section";
  Real Q_ea "evaporator->adiabatic section";
  Real Q_ac "adiabatic section->condenser";
  
  output Real Q_cs "condenser->Stirling CPS";
  Real R_vc "combined thermal resistance from vapour to condenser"; 
  Real R_va "combined thermal resistance from vapour to adiabatic section"; 
  //----------------------------
  // Inputs
  //----------------------------
  input Integer N_HPs_input;  
  input Real Q_evap_input;
  input Real Q_cond_input;
  input Real Q_stirling_input; //if modelled internally
  input Real T_stirling_input; //if modelled externally
  
  Real Q_evap_bc "heat addition rate to evaporator [W]";
  Real Q_cond_bc "heat removal rate from condenser [W]";
  Real Q_stirling_bc "heat removal rate from Stirling PCS [W]"; 
  Real N_HPs;
  
  
  //----------------------------
  // Model options -- these are not intended to be selectable by the FMU handler, but baked-in at compilation
  //----------------------------
  //Option to make set a fixed cold initial temperature in the whole heat pipe.
  parameter Boolean COLD_START = false;
  //Option to dynamically model the vapour core (giving a temp-dependent thermal resistance) mainly for startup behaviour.
  parameter Boolean MODEL_HP_STARTUP = true;
  //Option to model the core as a single lump or simply take Q_evap as a boundary condition (coupling) 
  parameter Boolean MODEL_CORE_INTERNALLY = false; 
  //Option to model the thermal mass of the Stirling PCS (whether internally or externally) or just take it as a BC on the condenser 
  parameter Boolean MODEL_STIRLING = true;   
  //Follow-up option to do so internally or via an external FMU 
  parameter Boolean MODEL_STIRLING_INTERNALLY = false; //leaving this false with MODEL_STIRLING true can cause the initialisation to fail within the OM environment as T_Stirling defaults to 0
  //Option to make the Stirling engine only activate above a set temperature  
  parameter Boolean MODEL_STIRLING_ACTIVATION = true;
      
  //----------------------------
  // Latches (flags with persistence)
  //----------------------------
  //Real STIRLING_ACTIVATED (start=0, fixed=true);
  Boolean STIRLING_ACTIVATED(start=false,fixed=true); //changed to a flag using 'when' statement.
    
  //----------------------------
  // Outputs / signals
  //----------------------------
  Boolean HP_activated "Flag raised when the heat pipe has effective heat transfer";
  
  //track the Reynolds number of the vapour throughout the simulation
  Real v_v "vapour velocity [m/s]" ;
  Real mdot_v "vapour mass flow rate [kg/s]";
  Real Re; 
  
  Real Q_lim_viscous; 
  Real Q_lim_flood; 
  
  input Real T_cond_init_input;   
  
  output Real T_evap_interface "actual evaporator interface temperature [K]";
  
  parameter Real R_evap_halfinternal = log((ro_wall-twall/2)/ri_wall)/(2*pi*Le*k_wall(800)); //roughly but not exactly half of R_wall_e  

  Real R_contact_eff; 

initial equation 

if COLD_START then 
  T_evap = 15 + 273.15;
  T_wick = 15 + 273.15;
  T_vap = 15 + 273.15;
  T_cond = 15 + 273.15;
  T_adiab = 15 + 273.15;
  
  T_stirling = 15 + 273.15;

else  //change the starting conditions to speed up convergence to steady-state in pseudotransient, full power cases
  
  der(T_evap) = 0;
  der(T_wick) = 0;
  der(T_vap) = 0;
  der(T_adiab) = 0;  
  
  if T_cond_init_input > 1E-9 then 
    T_cond = T_cond_init_input;
  else        
    //T_cond = T_stirling + R_stirling_hp_interface*(Q_evap_input/N_HPs); //fix the delta T as per definition of R_stirling_hp_interface
    der(T_cond) = 0;  
  end if;   
  
  if MODEL_STIRLING_INTERNALLY then     
    der(T_stirling) = 0;  
  end if; 
  
  
end if;

equation  
  //compute temp-dependent properties
  C_evap_wall = rho_wall * V_wall_e * cp_wall(T_evap);
  C_cond_wall = rho_wall * V_wall_c * cp_wall(T_cond);  
  //C_adiab_wall = rho_wall * V_wall_a * cp_wall(T_evap); //lumped in with evaporator  
  //C_adiab_wall = rho_wall * V_wall_a * cp_wall(T_cond); //lumped in with condenser  
  C_adiab_wall = rho_wall * V_wall_a * cp_wall(T_adiab); //Not lumped in 
  
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
  if R_vapour_ax < (R_wall_axial_ea + R_wall_axial_ac)/10. then 
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
    Q_evap_bc = Q_evap_input/N_HPs;  //must be allowed to cross zero and go negative 
  end if;
  
  if MODEL_STIRLING then 
    Q_cond_bc = 0;
    
    if MODEL_STIRLING_INTERNALLY and Q_stirling_input > 1E-9 then
      Q_stirling_bc = Q_stirling_input/N_HPs;
    else
      if MODEL_STIRLING_INTERNALLY and STIRLING_ACTIVATED then 
        Q_stirling_bc = HTC_cold*(T_stirling - T_stirling_cold_nominal);
        
      else 
        Q_stirling_bc = 0;
      end if;          
    end if;
    
  else 
    if abs(Q_cond_input) > 1E-9 then 
      Q_cond_bc = Q_cond_input/N_HPs;
    else       
      Q_cond_bc = max(0, HTC*(T_cond - T_stirling_nominal));  //Completely avoid adding heat to the condenser based on this HTC
      //Q_cond_bc = max(0, HTC*(T_cond^4 - T_stirling_nominal^4));  
    end if;
    Q_stirling_bc = 0;
  end if;  
  
  /*Heat transfer rates*/
  
  // Evaporator->wick
  Q_ew  = (T_evap - T_wick) / R_evap_radial;

  // Wick->vapour 
  Q_wv = (T_wick - T_vap) / R_wv;  

  // Vapour->adiabatic section 
  Q_lim_viscous = viscous_limit(T_vap);
  Q_lim_flood = flooding_limit(T_vap);
  
  R_va = max(R_adiab_radial, R_adiab_radial + R_vapour_ax/(La+Le/2.+Lc/2.)*(Le/2. + La/2.)); 
  Q_va = min(Q_lim_flood,min(Q_lim_viscous, (T_vap - T_adiab)/R_va)); 

  // Vapour->condenser 
  R_vc = max(R_cond_radial, R_cond_radial + R_vapour_ax);       
  Q_vc = min(Q_lim_flood, min(Q_lim_viscous, (T_vap - T_cond) / R_vc));  //Impose physical throughput limits
  
  /* //old simplified model 
  Evaporator->condenser
  Q_ec = (T_evap - T_cond) / R_wall_axial;
  */
  
  // Evaporator-> adiabatic section
  Q_ea = (T_evap - T_adiab) / (R_wall_axial_ea); 
  
  // Adiabatic section -> condenser
  Q_ac = (T_adiab - T_cond) / (R_wall_axial_ac); 
  
  // Condenser->Stirling. 
  //with thermal resistance according to prescribed nominal DeltaT
  Q_cs = max(1E-8,  (T_cond - T_stirling) / R_stirling_hp_interface );  
  //This max() protection is strictly needed to avoid a crash and I'm not sure why. Maybe because Q_internal in the Stirling model is strictly positive. Have not yet checked if the same error happens with the simpler Stirling model  
  
  
  
  
  //----------------------------
  // Energy balances
  //----------------------------
  
  // Core monolith - based on the full geometry. 
  if MODEL_CORE_INTERNALLY then 
    C_monolith * der(T_monolith) = 2350 - 8*Q_evap_bc;    
  else 
    T_monolith = 0; //placeholder value
  end if; 
  
  // Evaporator 
  //(C_evap_wall + C_adiab_wall) * der(T_evap) = Q_evap_bc - Q_ew - Q_ec; //with thermal mass of adiabatic section lumped in  
  //C_evap_wall * der(T_evap) = Q_evap_bc - Q_ew - Q_ec; 
  C_evap_wall * der(T_evap) = Q_evap_bc - Q_ew - Q_ea; 
  

  // Wick
  C_wick * der(T_wick) = Q_ew - Q_wv; //No axial conduction as there is no wick in the adiabatic/condenser sections in KRUSTY

  // Vapour
  C_vapour = rho_v*cp_v*V_vcore; //Update vapour capacitance based on density (rarefied/continuous state)
  C_vapour * der(T_vap) = Q_wv - Q_va - Q_vc; 
  
  // Adiabatic section 
  C_adiab_wall * der(T_adiab) = Q_ea + Q_va - Q_ac; 
  
  if MODEL_STIRLING_ACTIVATION then 
    when MODEL_STIRLING and T_stirling > T_stirling_activation then 
      STIRLING_ACTIVATED = true;
      Modelica.Utilities.Streams.print("Stirling engine activated!");
    end when;
  else 
    STIRLING_ACTIVATED = true;
  end if; 
  
  // Condenser and Stirling PCS
  if MODEL_STIRLING then        
    //C_cond_wall * der(T_cond) = Q_vc + Q_ec - Q_cs;    
    C_cond_wall * der(T_cond) = Q_vc + Q_ac - Q_cs;    
    //(C_cond_wall + C_adiab_wall) * der(T_cond) = Q_vc + Q_ec - Q_cs; //with thermal mass of adiabatic section lumped in     
    
      
    if MODEL_STIRLING_INTERNALLY then 
      C_stirling * der(T_stirling) = Q_cs - Q_stirling_bc; //Q will jump from zero when the engine is turned on . 
    else
      if T_stirling_input > 1E-9 then
        T_stirling = T_stirling_input; 
      else 
        T_stirling = T_stirling_nominal; //avoids a crash at initialisation if no Stirling temp received.         
      end if; 
    end if;  
  else 
    //C_cond_wall * der(T_cond) = Q_vc + Q_ec - Q_cond_bc;
    //(C_cond_wall + C_adiab_wall) * der(T_cond) = Q_vc + Q_ec - Q_cond_bc;
    C_cond_wall * der(T_cond) = Q_vc + Q_ac - Q_cond_bc;
    
    T_stirling = 0; //placeholder value    
    
  end if; 
  
  //Projection to the wall surface based on internal resistance
  //T_evap_interface = T_evap + R_evap_halfinternal; 
  //add effective contact resistance (will just push T_hp_wall higher on the other side). 
  T_evap_interface = T_evap + R_evap_halfinternal + R_contact_eff; 
  
  //R_contact_eff = 10*R_wall_e *exp(-(T_evap - 400) / 50); //R_wall_e
  R_contact_eff = 0;
end HeatPipeWithVapourCore;
/*
Implementation of Zuo & Faghri's seminal heat pipe model using standard MSL components, modified for a sodium heat pipe. 

Inputs: Integrated heat flux through evaporator wall, integrated heat flux out of condenser
Outputs: Evaporator wall temperature, condenser wall temperature 
*/
model HeatPipe_ZuoFaghri_MSLcomponents_alt
  import Modelica.Constants.pi;
  //----------------------------
  // Geometry (KRUSTY-like)
  //----------------------------
  parameter Integer N_HPs_nominal = 8;
  input Integer N_HPs_input;
  Real N_HPs;
  
  parameter Real Le = 0.25 "Evaporator length (m) -- KRUSTY";
  parameter Real Lc = 0.10 "Condenser length (m) -- KRUSTY, approx";
  parameter Real La = 1 - Le - Lc "Adiabatic length (m)";
  parameter Real L = Le + La + Lc;
  //Total length stated as approx. 100 cm
  parameter Real do = 0.0127 "Outer diameter (m) -- KRUSTY";
  parameter Real twall = 0.089*0.01 "Wall thickness (m), KRUSTY (just under 1 mm)";
  parameter Real ro = do/2 "Outer radius (m)";
  parameter Real wick_thickness = 1e-3 "Wick thickness m), guesstimate";
  parameter Real ri = ro - twall - wick_thickness "Inner radius (m)";
  parameter Real A_cond_outer = 2*pi*ro*Lc "condenser outer convective area (m2)";
  parameter Real A_evap_outer = 2*pi*ro*Le "evaporator outer convective area (m2)";
  parameter Real A_wick_outer_e = 2*pi*(ri + wick_thickness)*Le;
  parameter Real A_wick_outer_c = 2*pi*(ri + wick_thickness)*Lc;
  parameter Real A_wick_outer_a = 2*pi*(ri + wick_thickness)*La;
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
  parameter Real cp_wall = 595 "J/kg.K (Haynes-230 at ~800 C, datasheet)";
  parameter Real k_wall = 24 "W/m.K (Haynes-230 at ~800 C)";
  //----------------------------
  // Wick: nickel screen based (effective properties)
  //----------------------------
  parameter Real porosity = 0.70 "wick porosity (assumed)";
  // effective wick solid+liquid properties computed below (need cp/k/rho of void+solid)
  // choose effective solid (screen) assumed steel-like: use stainless-like numbers for screen base:
  parameter Real rho_screen = 8000 "kg/m3 (approx stainless screen)";
  parameter Real cp_screen = 500 "J/kg.K (screen)";
  parameter Real k_screen = 20 "W/m.K (screen)";
  // effective thermal conduction of wick (porous + liquid)
  parameter Real k_Na_const = 70 "W/m.K (approx liquid sodium thermal conductivity at high T)";
  //parameter Real k_wick = min(k_screen, k_Na_const) "W/m.K (effective porous wick thermal conductivity, conservative)";
  //parameter Real k_wick = k_Na_const*porosity + k_screen*(1-porosity); //simple average conductivity
  parameter Real k_wick = k_Na_const*((k_Na_const + k_screen) - (1 - porosity)*(k_Na_const - k_screen))/((k_Na_const + k_screen) + (1 - porosity)*(k_Na_const - k_screen));
  // Chi model (heat pipe theory)
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

  //----------------------------
  // Volumes and masses (explicit)
  //----------------------------
  // annular wick (approximated area at mean radius)
  parameter Real r_wick_mid = ri + wick_thickness/2;
  parameter Real A_wick_radial = pi*((ri + wick_thickness)^2 - ri^2);
  // annulus XS area
  parameter Real V_wick_e = A_wick_radial*Le;
  parameter Real V_wick_a = A_wick_radial*La;
  parameter Real V_wick_c = A_wick_radial*Lc;
  // wall volume (cylindrical shell)
  parameter Real A_wall_radial = pi*((ro)^2 - (ro - twall)^2);
  //XS area
  parameter Real V_wall_e = A_wall_radial*Le;
  parameter Real V_wall_a = A_wall_radial*La;
  parameter Real V_wall_c = A_wall_radial*Lc;
  // masses (solid wall + wick porous mass = screen solid fraction + liquid sodium fraction within wick)
  // approximate wick solid volume fraction = (1 - porosity)
  parameter Real m_wall_e = rho_wall*V_wall_e;
  parameter Real m_wall_a = rho_wall*V_wall_a;
  parameter Real m_wall_c = rho_wall*V_wall_c;
  // wick composite mass (solid screen + sodium)
  parameter Real V_wick_solid_e = (1 - porosity)*V_wick_e;
  parameter Real V_wick_liquid_e = porosity*V_wick_e;
  parameter Real V_wick_solid_a = (1 - porosity)*V_wick_a;
  parameter Real V_wick_liquid_a = porosity*V_wick_a;
  parameter Real V_wick_solid_c = (1 - porosity)*V_wick_c;
  parameter Real V_wick_liquid_c = porosity*V_wick_c;
  parameter Real m_wick_solid_e = rho_screen*V_wick_solid_e;
  parameter Real m_wick_liq_e = rho_Na_T(900)*V_wick_liquid_e;
  // approx T ~ 900K initial
  parameter Real m_wick_solid_a = rho_screen*V_wick_solid_a;
  parameter Real m_wick_liq_a = rho_Na_T(900)*V_wick_liquid_a;
  parameter Real m_wick_solid_c = rho_screen*V_wick_solid_c;
  parameter Real m_wick_liq_c = rho_Na_T(800)*V_wick_liquid_c;
  parameter Real m_wick_e = m_wick_solid_e + m_wick_liq_e;
  parameter Real m_wick_a = m_wick_solid_a + m_wick_liq_a;
  parameter Real m_wick_c = m_wick_solid_c + m_wick_liq_c;
  //----------------------------
  // Thermal capacitances (explicit Cp * mass)
  //----------------------------
  // wall heat capacitance
  parameter Real C1 = m_wall_e*cp_wall;
  parameter Real C4 = m_wall_c*cp_wall;
  parameter Real C6 = m_wall_a*cp_wall;
  // wick heat capacatance (composite: solid + liquid sodium)
  parameter Real cp_wick_eff_e = (m_wick_solid_e*cp_screen + m_wick_liq_e*cp_Na_T(900))/m_wick_e;
  parameter Real cp_wick_eff_a = (m_wick_solid_a*cp_screen + m_wick_liq_a*cp_Na_T(900))/m_wick_a;
  parameter Real cp_wick_eff_c = (m_wick_solid_c*cp_screen + m_wick_liq_c*cp_Na_T(800))/m_wick_c;
  parameter Real C2 = m_wick_e*cp_wick_eff_e;
  parameter Real C3 = m_wick_c*cp_wick_eff_c;
  parameter Real C5 = m_wick_a*cp_wick_eff_a;
  //----------------------------
  // Boundaries
  //----------------------------
  //Heat flux input, set as an input for the FMU. Note that if not set as a parameter, the provided value will be ignored.
  //parameter input Real Q_input = 5000.0/8 "W applied heat to evaporator (W) ";
  //Modelica.Blocks.Interfaces.RealInput Q_evap_interface;\
  input Real Q_cond_nominal_input;
  Real Q_cond_nominal;
  Real Q_evap_nominal;
  input Real Q_cond_input; 
  input Real Q_evap_input;
  Real Q_cond_current;
  Real Q_evap_current;
  //Modelica.Blocks.Interfaces.RealInput T_hot_input; 
  parameter Real T_inf = 50 + 273.15; //approximate from Poston et al Fig. 10
  Real R_stirling_hp_interface; //~170-200 C temperature drop from HP condenser to Stirling hot end
  parameter Real T_Stirling_nominal = 630 + 273.15; //Stirling hot end temperature (Poston et al Fig. 10)
  input Real T_Stirling_inst_input;
  Real T_cond_nominal; 
  //parameter Real HTC = -Q_cond_nominal/(T_cond_nominal - T_inf);
  Real HTC;
  //parameter Real HTC = 0.50; 
  
  /* Option A: model the Stirling engine with a fixed heat transfer coefficient */
  //parameter Real h_cond = 500 "W/m2.K condenser convection HTC (estimate)";
  //parameter Real T_inf = 20 + 273.15 "cold sink (Stirling hot side )";
  /* Option B: treat the hot end of the Stirling engine as a fixed temp boundary */
  //parameter Real h_cond = 1e8;
  //short-circuit the convective element
  //Define this BC as an "input" rather than a normal parameter, allowing it to be controlled externally
  //parameter input Real T_inf = 600.0 + 273.15;
  //Can use a normal input type (time-varying) if I set a prescribedTemperature boundary condition.
  //Also need to output the hot side temperature and the cold side heat flux.
  //output Real T_evap;
  output Real T_cond;
  output Real T_evap;
  //Real Q_hot; //no longer needs to be an output type
  //output Real Q_hot_allHPs; 
  //----------------------------
  // Ports
  //----------------------------
  //Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_a evapPort "Evaporator heat port";
  //Modelica.Thermal.HeatTransfer.Interfaces.HeatPort_b condPort "Condenser heat port";
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
  parameter Real C_1 = k_1*A_1/L_1;
  //conductances
  parameter Real C_2 = k_2*A_2/L_2;
  parameter Real C_3 = k_3*A_3/L_3;
  parameter Real C_4 = k_4*A_4/L_4;
  parameter Real C_5 = k_5*A_5/L_5;
  parameter Real C_6 = k_6*A_6/L_6;
  parameter Real R_1 = 1/C_1;
  parameter Real R_2 = 1/C_2;
  parameter Real R_3 = 1/C_3;
  parameter Real R_4 = 1/C_4;
  parameter Real R_5 = 1/C_5;
  parameter Real R_6 = 1/C_6;
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R1_1(R = R_1/2) annotation(
    Placement(transformation(origin = {-104, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R2_2(R = R_2/2) annotation(
    Placement(transformation(origin = {-22, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R3_1(R = R_3/2) annotation(
    Placement(transformation(origin = {4, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R4_1(R = R_4/2) annotation(
    Placement(transformation(origin = {56, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R5_1(R = R_5/2) annotation(
    Placement(transformation(origin = {-22, 34}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R6_1(R = R_6/2) annotation(
    Placement(transformation(origin = {-22, 70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor M1(C = C1, T(start = 1030.0, fixed = false, displayUnit = "degC"), der_T(start = 0, fixed = true)) annotation(
    Placement(transformation(origin = {-90, 8}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor M2(C = C2, T(start = 1030.0, fixed = false), der_T(start = 0, fixed = true)) annotation(
    Placement(transformation(origin = {-34, 8}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor M3(C = C3, T(start = 1030.0, fixed = false), der_T(start = 0, fixed = true)) annotation(
    Placement(transformation(origin = {18, 8}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor M4(C = C4, T(start = 1030.0, fixed = false), der_T(start = 0, fixed = true)) annotation(
    Placement(transformation(origin = {70, 6}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor M5(C = C5, T(start = 1030.0, fixed = false), der_T(start = 0, fixed = true)) annotation(
    Placement(transformation(origin = {-6, 50}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.HeatCapacitor M6(C = C6, T(start = 1030.0, fixed = false), der_T(start = 0, fixed = true)) annotation(
    Placement(transformation(origin = {-6, 86}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R1_2(R = R_1/2) annotation(
    Placement(transformation(origin = {-76, -13}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R2_1(R = R_2/2) annotation(
    Placement(transformation(origin = {-47, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R5_2(R = R_5/2) annotation(
    Placement(transformation(origin = {10, 34}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R6_2(R = R_6/2) annotation(
    Placement(transformation(origin = {11, 70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R3_2(R = R_3/2) annotation(
    Placement(transformation(origin = {31, -15}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R4_2(R = R_4/2) annotation(
    Placement(transformation(origin = {85, -13}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Q_cond annotation(
    Placement(transformation(origin = {119, -13}, extent = {{10, -10}, {-10, 10}})));
  Modelica.Thermal.HeatTransfer.Sources.PrescribedHeatFlow Q_evap annotation(
    Placement(transformation(origin = {-138, -14}, extent = {{-9, -9}, {9, 9}})));
equation
  if N_HPs_input > 0 then
    N_HPs = N_HPs_input;
  else 
    N_HPs = N_HPs_nominal;
  end if;

  if Q_cond_nominal_input > 1e-6 then
    Q_cond_nominal = -1.*Q_cond_nominal_input/N_HPs;
  else
    Q_cond_nominal = -2350.0/N_HPs;
  end if;  
  Q_evap_nominal = -1*Q_cond_nominal;  
  R_stirling_hp_interface = 165/(-Q_cond_nominal); //~145 C measured temperature drop from HP condenser to Stirling hot end
  //Can actually tune the temperature drop to get the desired T_avg
  T_cond_nominal = T_Stirling_nominal + (-Q_cond_nominal)*R_stirling_hp_interface; 
  HTC = -Q_cond_nominal/(T_cond_nominal - T_Stirling_nominal);

  //Q_cond_current = -Q_cond_input/N_HPs; //-1.*HTC*(T_cond - T_inf);
  //Q_cond_current = -1.*HTC*(T_cond - T_inf);  
  //Q_evap_current = 2350.0/N_HPs;
  
  if Q_evap_input > 1e-6 then
    Q_evap_current = Q_evap_input/N_HPs; 
  else
    Q_evap_current = Q_evap_nominal;
  end if;  

  //if Q_cond_input > 1e-6 then
  //  Q_cond_current = -Q_cond_input/N_HPs;
  /*
  if T_Stirling_inst_input > 0 then //instead of prescribing Q_cond directly, vary T_stirling         
    if (Q_cond_input > 0) then //But Q_cond_input takes precedent
      Q_cond_current = -Q_cond_input/N_HPs; 
    else
      Q_cond_current = -1.*HTC*(T_cond - T_Stirling_inst_input);    
    end if;
  else //in case there is no driven Q_cond flow, mimic the action of a PCS.     
    //Q_cond_current = -1.*HTC*(T_cond - T_inf); 
    if (Q_cond_input > 0) then 
      Q_cond_current = -Q_cond_input/N_HPs; 
    else
      Q_cond_current = -1.*HTC*(T_cond - T_Stirling_nominal); 
    end if;
  end if;
  */
  
  if (Q_cond_input > 0) then //Q_cond_input takes precednet
    Q_cond_current = -Q_cond_input/N_HPs; 
  else
    if T_Stirling_inst_input > 0 then //instead of prescribing Q_cond directly, vary T_stirling         
      Q_cond_current = -1.*HTC*(T_cond - T_Stirling_inst_input);    
    else
      Q_cond_current = -1.*HTC*(T_cond - T_Stirling_nominal); //mainly for steady-state initialisation
    end if;
  end if;
  
  T_evap = M1.T;
  T_cond = M4.T;
  //Q_hot = R1_1.Q_flow;
  //Q_hot_allHPs = Q_hot*N_HPs;
  
//manually adjust for the number of heat pipes
//copy the evaporator wall temperature to this output variable.
//Q_out = Convection.Q_flow;
//R4_2.Q_flow; //Get the Qflow out through the cold side.
  connect(R1_1.port_b, M1.port) annotation(
    Line(points = {{-94, -14}, {-94, -16}, {-90, -16}, {-90, -2}}, color = {191, 0, 0}));
  connect(R1_2.port_a, M1.port) annotation(
    Line(points = {{-86, -13}, {-86, -15}, {-90, -15}, {-90, -2}}, color = {191, 0, 0}));
  connect(R1_2.port_b, R2_1.port_a) annotation(
    Line(points = {{-66, -13}, {-56, -13}, {-56, -15}}, color = {191, 0, 0}));
  connect(R2_1.port_b, M2.port) annotation(
    Line(points = {{-37, -14}, {-37, -9}, {-34, -9}, {-34, -2}}, color = {191, 0, 0}));
  connect(R2_2.port_a, M2.port) annotation(
    Line(points = {{-32, -14}, {-32, -12}, {-34, -12}, {-34, -2}}, color = {191, 0, 0}));
  connect(R2_2.port_b, R3_1.port_a) annotation(
    Line(points = {{-12, -14}, {-6, -14}}, color = {191, 0, 0}));
  connect(R3_1.port_b, M3.port) annotation(
    Line(points = {{14, -14}, {14, -12}, {18, -12}, {18, -2}}, color = {191, 0, 0}));
  connect(R3_2.port_a, M3.port) annotation(
    Line(points = {{21, -15}, {21, -2}, {18, -2}}, color = {191, 0, 0}));
  connect(R3_2.port_b, R4_1.port_a) annotation(
    Line(points = {{41, -15}, {45, -15}}, color = {191, 0, 0}));
  connect(R4_1.port_b, M4.port) annotation(
    Line(points = {{66, -14}, {70, -14}, {70, -4}}, color = {191, 0, 0}));
  connect(R4_2.port_a, M4.port) annotation(
    Line(points = {{75, -13}, {69, -13}, {69, -5}}, color = {191, 0, 0}));
  connect(R6_1.port_b, M6.port) annotation(
    Line(points = {{-12, 70}, {-6, 70}, {-6, 76}}, color = {191, 0, 0}));
  connect(R6_2.port_a, M6.port) annotation(
    Line(points = {{2, 70}, {-6, 70}, {-6, 76}}, color = {191, 0, 0}));
  connect(R5_1.port_b, M5.port) annotation(
    Line(points = {{-12, 34}, {-6, 34}, {-6, 40}}, color = {191, 0, 0}));
  connect(R5_2.port_a, M5.port) annotation(
    Line(points = {{0, 34}, {-6, 34}, {-6, 40}}, color = {191, 0, 0}));
  connect(R6_1.port_a, M1.port) annotation(
    Line(points = {{-32, 70}, {-32, 68}, {-90, 68}, {-90, -2}}, color = {191, 0, 0}));
  connect(R5_2.port_b, M4.port) annotation(
    Line(points = {{20, 34}, {70, 34}, {70, -4}}, color = {191, 0, 0}));
  connect(R6_2.port_b, M4.port) annotation(
    Line(points = {{22, 70}, {70, 70}, {70, -4}}, color = {191, 0, 0}));
  connect(R5_1.port_a, M1.port) annotation(
    Line(points = {{-32, 34}, {-32, 32}, {-90, 32}, {-90, -2}}, color = {191, 0, 0}));
  connect(R4_2.port_b, Q_cond.port) annotation(
    Line(points = {{96, -12}, {110, -12}}, color = {191, 0, 0}));
  connect(Q_evap.port, R1_1.port_a) annotation(
    Line(points = {{-129, -14}, {-114, -14}}, color = {191, 0, 0}));
  
//connect(Q_cond_current, Q_cond.Q_flow);
  Q_cond.Q_flow = Q_cond_current;
  //connect(T_hot_input, MonolithWallTemp.T);
/*Switch to the below equations in case we want to simulate within the OpenModelica env*/
//Q_cond.Q_flow = -(5000./8);
//MonolithWallTemp.T = 600 + 273.15;
//old stuff
//connect(evapPort, Q_evap.port);
  Q_evap.Q_flow = Q_evap_current;
//connect(Q_evap_interface, Q_evap.Q_flow);
//Q_evap.Q_flow = Q_input;
//Q_evap_interface = Q_input;
  annotation(
    uses(Modelica(version = "4.0.0")));
end HeatPipe_ZuoFaghri_MSLcomponents_alt;
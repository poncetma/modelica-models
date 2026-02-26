model HeatPipe_ZuoFaghri_Massless
  import Modelica.Constants.pi;
  
  
    //----------------------------
  // Geometry (KRUSTY-like)
  //----------------------------
  parameter Real Le = 0.25 "Evaporator length (m) -- KRUSTY";  
  parameter Real Lc = 0.25 "Condenser length (m) -- KRUSTY, approx";
  parameter Real La = 1 - Le - Lc "Adiabatic length (m)";
  parameter Real L = Le + La + Lc; 
//Total length stated as approx. 100 cm
  // outer diameter ~ 5/8 in = 15.875 mm (value reported in INL microreactor example)
  parameter Real do = 0.0127; //0.015875 "Outer diameter (m) (5/8 in)";
  parameter Real twall = 0.089*0.01 "Wall thickness (m), KRUSTY (just under 1 mm)";
  parameter Real ro = do/2 "Outer radius (m)";
  parameter Real wick_thickness = 1e-3 "Wick thickness nominal (m)";
  parameter Real ri = ro - twall - wick_thickness "Inner radius (m); assume wick annulus ~1 mm thick inside wall, guesstimate";  

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
  parameter Real eps = 0.70 "wick porosity (assumed)";
  // effective wick solid+liquid properties computed below (need cp/k/rho of void+solid)
  // choose effective solid (screen) assumed steel-like: use stainless-like numbers for screen base:
  parameter Real rho_screen = 8000 "kg/m3 (approx stainless screen)";
  parameter Real cp_screen  = 500  "J/kg.K (screen)";
  parameter Real k_screen   = 20   "W/m.K (screen)";

  // effective thermal conduction of wick (porous + liquid)
  parameter Real k_Na_const = 70 "W/m.K (approx liquid sodium thermal conductivity at high T)";
  
  //parameter Real k_wick = min(k_screen, k_Na_const) "W/m.K (effective porous wick thermal conductivity, conservative)";
  parameter Real k_wick = k_Na_const*eps + k_screen*(1-eps); //simple average conductivity

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
  // wall heat capacitance
  parameter Real C1 = m_wall_e * cp_wall;  
  parameter Real C4 = m_wall_c * cp_wall;
  parameter Real C6 = m_wall_a * cp_wall;

  // wick heat capacatance (composite: solid + liquid sodium)
  parameter Real cp_wick_eff_e = (m_wick_solid_e*cp_screen + m_wick_liq_e*cp_Na_T(900)) / m_wick_e;
  parameter Real cp_wick_eff_a = (m_wick_solid_a*cp_screen + m_wick_liq_a*cp_Na_T(900)) / m_wick_a;
  parameter Real cp_wick_eff_c = (m_wick_solid_c*cp_screen + m_wick_liq_c*cp_Na_T(800)) / m_wick_c;

  parameter Real C2 = m_wick_e * cp_wick_eff_e;  
  parameter Real C3 = m_wick_c * cp_wick_eff_c;
  parameter Real C5 = m_wick_a * cp_wick_eff_a;

  //----------------------------
  // Heat transfer boundary (condenser convection)
  //----------------------------
  parameter Real Q_input = 5000.0/8 "W applied heat to evaporator (W) -- KRUSTY-scale reactor ~5 kW thermal divided by 8 heat pipes";
  parameter Real h_cond = 50 "W/m2.K condenser convection HTC (estimate)";
  parameter Real T_inf = 20 + 273.15 "cold sink (Stirling hot side )";

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

  parameter Real C_1 = k_1*A_1/L_1; //conductances
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


  Modelica.Thermal.HeatTransfer.Sources.FixedHeatFlow Q_evap(Q_flow = Q_input)  annotation(
    Placement(transformation(origin = {-132, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R1_1(R = 1/C_1/2)  annotation(
    Placement(transformation(origin = {-104, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R2_2(R = 1/C_2/2)  annotation(
    Placement(transformation(origin = {-22, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R3_1(R = 1/C_3/2)  annotation(
    Placement(transformation(origin = {4, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R4_1(R = 1/C_4/2)  annotation(
    Placement(transformation(origin = {56, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R5_1(R = 1/C_5/2)  annotation(
    Placement(transformation(origin = {-22, 34}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R6_1(R = 1/C_6/2)  annotation(
    Placement(transformation(origin = {-22, 70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Sources.FixedTemperature boundaryTemp(T = T_inf)  annotation(
    Placement(transformation(origin = {144, -12}, extent = {{8, -8}, {-8, 8}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor Convection(R = 1/(h_cond*A_cond_outer))  annotation(
    Placement(transformation(origin = {114, -12}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R1_2(R = 1/C_1/2) annotation(
    Placement(transformation(origin = {-76, -13}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R2_1(R = 1/C_2/2) annotation(
    Placement(transformation(origin = {-47, -14}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R5_2(R = 1/C_5/2) annotation(
    Placement(transformation(origin = {10, 34}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R6_2(R = 1/C_6/2) annotation(
    Placement(transformation(origin = {11, 70}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R3_2(R = 1/C_3/2) annotation(
    Placement(transformation(origin = {31, -15}, extent = {{-10, -10}, {10, 10}})));
  Modelica.Thermal.HeatTransfer.Components.ThermalResistor R4_2(R = 1/C_4/2) annotation(
    Placement(transformation(origin = {85, -13}, extent = {{-10, -10}, {10, 10}})));
equation
  connect(Convection.port_b, boundaryTemp.port) annotation(
    Line(points = {{124, -12}, {136, -12}}, color = {191, 0, 0}));
  connect(Q_evap.port, R1_1.port_a) annotation(
    Line(points = {{-122, -14}, {-114, -14}}, color = {191, 0, 0}));
  connect(R1_2.port_b, R2_1.port_a) annotation(
    Line(points = {{-66, -13}, {-56, -13}, {-56, -15}}, color = {191, 0, 0}));
  connect(R2_2.port_b, R3_1.port_a) annotation(
    Line(points = {{-12, -14}, {-6, -14}}, color = {191, 0, 0}));
  connect(R3_2.port_b, R4_1.port_a) annotation(
    Line(points = {{41, -15}, {45, -15}}, color = {191, 0, 0}));
  connect(R4_2.port_b, Convection.port_a) annotation(
    Line(points = {{95, -13}, {103, -13}}, color = {191, 0, 0}));
  connect(R1_1.port_b, R1_2.port_a) annotation(
    Line(points = {{-94, -14}, {-88, -14}, {-88, -12}, {-86, -12}}, color = {191, 0, 0}));
  connect(R5_1.port_a, R1_1.port_b) annotation(
    Line(points = {{-32, 34}, {-94, 34}, {-94, -14}}, color = {191, 0, 0}));
  connect(R6_1.port_a, R1_1.port_b) annotation(
    Line(points = {{-32, 70}, {-94, 70}, {-94, -14}}, color = {191, 0, 0}));
  connect(R5_2.port_b, R4_2.port_a) annotation(
    Line(points = {{20, 34}, {76, 34}, {76, -12}}, color = {191, 0, 0}));
  connect(R4_2.port_a, R4_1.port_b) annotation(
    Line(points = {{76, -12}, {66, -12}, {66, -14}}, color = {191, 0, 0}));
  connect(R6_2.port_b, R4_2.port_a) annotation(
    Line(points = {{22, 70}, {76, 70}, {76, -12}}, color = {191, 0, 0}));
  connect(R5_1.port_b, R5_2.port_a) annotation(
    Line(points = {{-12, 34}, {0, 34}}, color = {191, 0, 0}));
  connect(R6_1.port_b, R6_2.port_a) annotation(
    Line(points = {{-12, 70}, {2, 70}}, color = {191, 0, 0}));
  connect(R3_1.port_b, R3_2.port_a) annotation(
    Line(points = {{14, -14}, {22, -14}}, color = {191, 0, 0}));
  annotation(
    uses(Modelica(version = "4.0.0")));
end HeatPipe_ZuoFaghri_Massless;
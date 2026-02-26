model DummyTHKRUSTY
/*
 Node-based implementation of radial heat conduction in KRUSTY with a uniform source. 
 Inputs: integral power, HP wall heat flux
 Outputs: T_wall, T_avg
*/
import Modelica.Constants.pi;

  // --- Physical Parameters (KRUSTY U-10Mo) ---
  final constant Integer N = 20 "Number of radial shells";
  parameter Real r_inner = 0.02 "Inner fuel radius [m]";
  parameter Real r_outer = 0.055 "Outer fuel radius [m]";
  parameter Real L = 0.25 "Core height [m]";
  parameter Real fuel_volume = pi*(r_outer^2. - r_inner^2.)*L;
  
  parameter Real rho = 17110 "Density [kg/m3]";
  parameter Real cp = 189; //140 "Specific heat capacity [J/kg.K]";
  parameter Real k = 37.5; //20 "Thermal conductivity [W/m.K]";
  
// --- Heat Pipe / Boundary Parameters ---
  parameter Real T_HP_steady = 1073.5; //873.15 "Fixed HP temperature for steady-state [K]";
  parameter Integer n_hp = 8 "Number of heat pipes";
  parameter Real P_total_nom = 2350.0 "Nominal total thermal power [W]";
  Real P_integral; 
  input Real P_integral_input "Instantaneous integral power [W]";
  Real q_flux_hp;
  input Real Q_hp_input "Q flow through heat pipe boundary [W]";
  parameter Real w_contact = 0.015 "Effective contact width per heat pipe [m]";
  
  // --- Calculated Operational Values ---
  //parameter Real P_per_hp = P_integral / n_hp "Power per heat pipe [W]";
  parameter Real A_hp_eff = n_hp * w_contact * L "Total effective HP contact area [m2]";
  //parameter Real q_flux_hp = P_integral / A_hp_eff "Flux at the HP interface [W/m2]";
  Real q_gen "Required Q''' [W/m3]";
  
  // --- Variables ---
  //Real T[N](start=fill(800, N)) "Temperature at shell centers [K]";
  Real T[N] "Temperature at shell centers [K]";
  Real Q_flow[N+1] "Heat flow rate across shell faces [W]";
  Real V_shell[N] "Volume of each shell [m3]";
  
  //output Real T_outer_wall "Temperature at the heat pipe interface [K]";
  output Real T_mean "Volume-weighted mean temperature [K]";

  // --- Geometry Constants ---
  parameter Real dr = (r_outer - r_inner) / N;
  parameter Real r_face[N+1] = {r_inner + (i-1)*dr for i in 1:N+1};

initial equation
  // Force steady-state at t=0
  //T_outer_wall = T_HP_steady;    
  //q_flux_hp = P_total_nom/A_hp_eff;
  Q_flow[1] = 0; // Inner Boundary: Adiabatic (Control Rod Hole) 
  for i in 1:N loop //skip the first element since we set it adiabatic no matter what
    der(T[i]) = 0;
  end for;

equation
  P_integral = P_total_nom; 
  //P_integral = P_integral_input;   
  
  q_flux_hp = P_total_nom/A_hp_eff; 
  //q_flux_hp = Q_hp_input/A_hp_eff;  
  
  q_gen = P_integral / fuel_volume; //uniform heating

  // 1. Heat Flows at Faces
  Q_flow[1] = 0; // Inner Boundary: Adiabatic (Control Rod Hole) 
  
  // 1. Boundary Condition Logic
  // At t=0, we fix the Temperature. The solver calculates the required Flux.
  // At t>0, we fix the Flux. The solver calculates the resulting Temperature.  
  Q_flow[N+1] = q_flux_hp*A_hp_eff; 
  
  /* //This structure also breaks the FMU generation. Have to put the first equation in the actual initial equation block.
  if initial() then
    T_outer_wall = T_HP_steady;    
  else 
    T_outer_wall = T[N] - Q_flow[N+1]/A_hp_eff * (dr/2) / k;
  end if;
  */
  //T_outer_wall = T[N] - Q_flow[N+1]/A_hp_eff * (dr/2) / k;
  
  // 2. Energy Balance for each Shell
  //inner Q_flows
  for i in 2:N loop
    Q_flow[i] = -k * (2 * pi * r_face[i] * L) * (T[i] - T[i-1]) / dr;
  end for;
  for i in 1:N loop
    V_shell[i] = pi * (r_face[i+1]^2 - r_face[i]^2) * L;
    rho * V_shell[i] * cp * der(T[i]) = Q_flow[i] - Q_flow[i+1] + q_gen * V_shell[i];
  end for;  
 //This algorithm breaks the FMU generation, for some reason. 
 /*
algorithm
  // 4. Mean Temperature Calculation
  
  T_mean := 0;
  for i in 1:N loop
    T_mean := T_mean + T[i] * V_shell[i];
  end for;
  T_mean := T_mean / (pi * (r_outer^2 - r_inner^2) * L);
  */
  T_mean = sum(T[i] * V_shell[i] for i in 1:N) / fuel_volume;
  
  annotation(uses(Modelica(version="4.0.0")));
end DummyTHKRUSTY;
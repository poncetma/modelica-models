model SimpleTHKRUSTY_Lossy
/*
 Node-based implementation of radial heat conduction in KRUSTY with a uniform source. 
 This version includes fixed thermal losses out of the core (heat not going through the heat pipes).
 
 Inputs: integral power, T_wall
 Outputs: HP wall heat flux, T_avg
*/
import Modelica.Constants.pi;

  // --- Physical Parameters (KRUSTY U-10Mo) ---
  final constant Integer N = 50 "Number of radial shells"; //tested up to 500
  parameter Real r_inner = 0.02 "Inner fuel radius [m]";
  //parameter Real r_outer = 0.055 "Outer fuel radius [m]";
  parameter Real L = 0.25 "Core height [m]";
  //parameter Real fuel_volume = pi*(r_outer^2. - r_inner^2.)*L;
  
  //The fuel volume should  match exactly: 0.00186454 m3
  parameter Real fuel_volume = 0.00186454; //computed in OpenFOAM from KRUSTY mesh. 
  //with fixed inner radius, compute the appropriate equivalent outer radius
  parameter Real r_outer = sqrt(r_inner^2 + fuel_volume/L/pi);
  
  parameter Real rho = 17110 "Density [kg/m3]";
  parameter Real cp = 189; //140 "Specific heat capacity [J/kg.K]";
  parameter Real k = 37.5; //20 "Thermal conductivity [W/m.K]";
  
// --- Heat Pipe / Boundary Parameters ---
  parameter Real T_HP_nominal = 1083.3074675962812; //1073.5; //873.15 "Fixed HP temperature for steady-state [K]";
  parameter Integer n_hp = 8 "Number of heat pipes";
  parameter Real P_total_nom = 2350.0  "Nominal total thermal power [W]";
  input Real Q_loss_nominal_input;
  input Real Q_loss_input;
  Real Q_loss "Heat lost through MLI [W]";
  Real P_integral; 
  input Real P_integral_input "Instantaneous integral power [W]";
  //Real q_flux_hp;
  //input Real Q_hp_input "Q flow through heat pipe boundary [W]";
  input Real T_outer_wall_input "Temperature of heat pipe wall [K]";
  
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
  
  Real T_outer_wall "Temperature at the heat pipe interface [K]";
  output Real Q_evap_out "Heat flow out of heat pipe wall [W]";
  output Real T_mean "Volume-weighted mean temperature [K]";

  // --- Geometry Constants ---
  parameter Real dr = (r_outer - r_inner) / N;
  parameter Real r_face[N+1] = {r_inner + (i-1)*dr for i in 1:N+1};

initial equation
  // Force steady-state at t=0
  //T_outer_wall = T_HP_steady;   
  //Q_flow[1] = 0;   
  //T_outer_wall = T_outer_wall_input;
  
  for i in 1:N loop //Can no longer skip the first element
    der(T[i]) = 0;
  end for;

equation
  if Q_loss_input > 1e-6 then 
    Q_loss = Q_loss_input;
  else
    Q_loss = Q_loss_nominal_input;
  end if;

  //P_integral = P_total_nom; 
  //P_integral = P_integral_input;   
  if P_integral_input > 1e-6 then
    P_integral = P_integral_input;
  else
    P_integral = P_total_nom + Q_loss;
  end if;
  
  //q_flux_hp = P_total_nom/A_hp_eff; 
  //q_flux_hp = Q_hp_input/A_hp_eff;  

  q_gen = P_integral / fuel_volume; //uniform heating

  // 1. Heat Flows at Faces
  Q_flow[1] = 0; // Inner Boundary: Adiabatic (Control Rod Hole) 
  
  // 1. Boundary Condition Logic
  // At t=0, we fix the Temperature. The solver calculates the required Flux.
  // At t>0, we fix the Flux. The solver calculates the resulting Temperature.  
  //Q_flow[N+1] = q_flux_hp*A_hp_eff; 
  
  if T_outer_wall_input > 1e-6 then
    T_outer_wall = T_outer_wall_input;
  else 
    T_outer_wall = T_HP_nominal;
  end if;
  //Q_flow[N+1] = P_total_nom; 
  //T_outer_wall = T_HP_nominal;
  //T_outer_wall = T[N] - Q_flow[N+1]/A_hp_eff * (dr/2) / k;
  
  // 2. Energy Balance for each Shell
  //inner Q_flows --positive when heat flows towards +ve r (negative T gradient)
  for i in 2:N loop
    Q_flow[i] = -k * (2 * pi * r_face[i] * L) * (T[i] - T[i-1]) / dr; 
  end for;
  //Q_flow[N+1] = -k * (2 * pi * r_face[N+1] * L) * (T_outer_wall - T[N]) / (dr/2); //half-dr due to cell-centredness
  //T_outer_wall = T[N] - (Q_flow[N+1] / (2 * pi * r_face[N+1] * L)) * (dr/2) / k;
  //Q_flow[N+1] = -k * (A_hp_eff) * (T_outer_wall - T[N]) / (dr/2); //half-dr due to cell-centredness
  //Now there are two parallel heat flow paths at the outer wall. 
  Q_flow[N+1] = -k * (A_hp_eff) * (T_outer_wall - T[N]) / (dr/2) + Q_loss; 
  
  for i in 1:N loop
    V_shell[i] = pi * (r_face[i+1]^2 - r_face[i]^2) * L;
    rho * V_shell[i] * cp * der(T[i]) = Q_flow[i] - Q_flow[i+1] + q_gen * V_shell[i];
  end for;  
 //This algorithm breaks the FMU generation, for some reason. 

  T_mean = sum(T[i] * V_shell[i] for i in 1:N) / fuel_volume;
  Q_evap_out = Q_flow[N+1] - Q_loss;
  
  annotation(uses(Modelica(version="4.0.0")));
end SimpleTHKRUSTY_Lossy;
model SimpleTHKRUSTY_Lossy_nonuniformq
/*
 Control-volume-based implementation of 1D heat conduction in KRUSTY with a uniform source. 
 Currently based on Cartesian coordinates.
 This version includes fixed thermal losses out of the core (heat not going through the heat pipes).
 
 Inputs: integral power (heat deposition rate from fission + decay heat), T_wall
 Outputs: HP wall heat flux, T_avg
*/
import Modelica.Constants.pi;

  // --- Physical Parameters (KRUSTY U-10Mo) ---
  final constant Integer N = 65 "Number of radial shells"; //tested up to 500
  parameter Real r_inner = 0.02 "Inner fuel radius [m]";  
  parameter Real L = 0.25 "Core height [m]";   
  parameter Real fuel_volume = 0.00186454; //computed in OpenFOAM from KRUSTY mesh.   
  parameter Real r_outer = sqrt(r_inner^2 + fuel_volume/L/pi); //with fixed inner radius, compute the appropriate equivalent outer radius  
  
  //correlations for mat properties
  function cp_correlation
    input Real T_local "Temperature (K)";
    output Real cp "isobaric specific heat (J/kg.K)";
    protected Real T_celcius_local;
  algorithm
    T_celcius_local := T_local - 273.15;    
    //cp := (0.00007333333*(1091.45-273.5) + 0.134666667)*1000; //own correlation from Burkes data 
    cp := (0.00007333333*(T_celcius_local) + 0.134666667)*1000; //own correlation from Burkes data 
    //cp := 189;
  end cp_correlation;
  
  function k_correlation
    input Real T_local "Temperature (K)";
    output Real k "isobaric specific heat (J/kg.K)";
    protected Real T_celcius_local;
  algorithm
    T_celcius_local := T_local - 273.15;
    //k := 10.2 + (3.51E-2)*T_celcius_local;     //for some reason this is causing instabilities
    k := 10.2 + (3.51E-2)*(1091.45-273.15);     
  end k_correlation;
  
  function rho_correlation
    input Real T_local "Temperature (K)";
    output Real rho "isobaric specific heat (J/kg.K)";
    protected Real T_celcius_local;
  algorithm
    T_celcius_local := T_local - 273.15;        
    rho := (17.15 - (8.63E-4)*(T_celcius_local + 20.))*1000.0; //in kg/m3
    //rho := (17.15 - 8.63E-4*((1091.45-273.25) + 20.))*1000.0; //in kg/m3
    //rho := 17110; //Room temp value is actually 17320, 98.5% of theoretical 17580    
  end rho_correlation;
  
  // --- Heat Pipe / Boundary Parameters ---
  parameter Real T_HP_nominal = 1073.3074676166252 "Fixed HP temperature for steady-state [K]"; 
  parameter Integer n_hp = 8 "Number of heat pipes";
  parameter Real P_total_nom = 2350.0  "Nominal total thermal power [W]";
  input Real Q_loss_nominal_input;
  input Real Q_loss_input;
  Real Q_loss "Heat lost through MLI [W]";
  Real P_integral; 
  input Real P_integral_input "Instantaneous integral power [W]";  
  input Real T_outer_wall_input "Temperature of heat pipe wall [K]";
  
  parameter Real w_contact = 0.015 "Effective contact width per heat pipe [m]";  
  parameter Real A_hp_eff = n_hp * w_contact * L "Total effective HP contact area [m2]";
  Real q_gen[N] "Required Q''' [W/m3]";   
  Real q_gen_prof[N]; //Define q_gen based on a prescribed radial profile
  
  Real T[N] "Temperature at shell centers [K]";
  Real Q_flow[N+1] "Heat flow rate across shell faces [W]";
  Real V_shell[N] "Volume of each shell [m3]";
  
  Real T_outer_wall "Temperature at the heat pipe interface [K]";
  output Real Q_evap_out "Heat flow out of heat pipe wall [W]";
  output Real T_mean "Volume-weighted mean temperature [K]";

  //--Nodalisation
  parameter Real dr = (r_outer - r_inner) / N;
  final parameter Real r_face[N+1] = {r_inner + (i-1)*dr for i in 1:N+1};
  
  
  Real power_profile_integral;
  Real verified_power_integral;
  
initial equation
  for i in 1:N loop 
    der(T[i]) = 0;
  end for;
equation

  if Q_loss_input > 1e-6 then 
    Q_loss = Q_loss_input;
  else
    Q_loss = Q_loss_nominal_input;
  end if;  
  
  if P_integral_input > 1e-6 then
    P_integral = P_integral_input;
  else
    P_integral = P_total_nom + Q_loss;
  end if;

  for i in 1:N loop
    //follow radial profile and renormalise to nominal integral power
    if i < 9./10.*N then
      q_gen_prof[i] = 1;
    else
      q_gen_prof[i] = 1; //1 reverts back to normal profile   
    end if;
    V_shell[i] = pi * (r_face[i+1]^2 - r_face[i]^2) * L;
  end for;   
  power_profile_integral = sum(q_gen_prof[i] * V_shell[i] for i in 1:N);
  
  
  Q_flow[1] = 0; // Inner Boundary: Adiabatic (Control Rod Hole) 
  
  // The temperature at the evap wall is fixed, received from HP solver. 
  if T_outer_wall_input > 1e-6 then
    T_outer_wall = T_outer_wall_input;
  else 
    T_outer_wall = T_HP_nominal;
  end if;
  
  // Energy Balance for each Shell
  // inner Q_flows is positive when heat flows towards positive r (negative T gradient)
  for i in 2:N loop
    Q_flow[i] = -k_correlation((T[i-1] + T[i])/2.) * (2 * pi * r_face[i] * L) * (T[i] - T[i-1]) / dr; 
  end for;  
  //Assume Q_loss is constant (true for small variation of outer wall temp)
  //Q_flow[N+1] = -k_correlation((T[N]+T_outer_wall)/2.) * (A_hp_eff) * (T_outer_wall - T[N]) / (dr[N]/2) + Q_loss; 
  
  // Outer boundary: MUST use the actual interfacial area to have a consistent 2nd order scheme!
  // This affects the solution but have shown that it's not a radical difference
  Q_flow[N+1] = -k_correlation((T[N]+T_outer_wall)/2.) * (2*pi*r_face[N+1]*L) * (T_outer_wall - T[N]) / (dr/2);//+Q_loss; 

  for i in 1:N loop    
    q_gen[i] = q_gen_prof[i]/power_profile_integral*P_integral; //update q_gen
    rho_correlation(T[i]) * V_shell[i] * cp_correlation(T[i]) * der(T[i]) = Q_flow[i] - Q_flow[i+1] + q_gen[i] * V_shell[i];    
  end for;  
  verified_power_integral = sum(q_gen[i] * V_shell[i] for i in 1:N);
 
  
  T_mean = sum(T[i] * V_shell[i] for i in 1:N) / fuel_volume;
  Q_evap_out = Q_flow[N+1] - Q_loss;
  
  annotation(uses(Modelica(version="4.0.0")));
end SimpleTHKRUSTY_Lossy_nonuniformq;
model conductionKRUSTYRadialDisc
/*
 Finite-difference implementation of 1D radial heat conduction in KRUSTY with a possibly non-uniform source. 
 The cylindrical coordinate discretisation accounts for non-uniform and time-varying conductivity. 
 
 This version implicitly accounts thermal losses out of the core (heat not going through the heat pipes).
 
 Inputs: total integral power (heat deposition rate from fission + decay heat), evaporator wall temperature, heat loss rate 
 Outputs: HP wall heat flux, Core average temperature
*/
import Modelica.Constants.pi;
  //input Real Q_loss_nominal_input;
  input Real Q_loss_input;
  input Real P_integral_input "Instantaneous integral power [W]";  
  input Real T_outer_wall_input "Temperature of heat pipe wall [K]";  
  
  // --- Physical Parameters (KRUSTY U-10Mo) ---
  final constant Integer N = 80 "Number of radial shells"; //tested up to 500
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
    cp := (0.00007333333*(T_celcius_local) + 0.134666667)*1000; //own correlation from Burkes data 
    //cp := 189;
  end cp_correlation;
  
  function k_correlation
    input Real T_local "Temperature (K)";
    output Real k "isobaric specific heat (J/kg.K)";
    protected Real T_celcius_local;
  algorithm
    T_celcius_local := T_local - 273.15;
    k := 10.2 + (3.51E-2)*T_celcius_local; //With the old discretisation, this caused instabilities. Now with rigorously derived version, no issues. 
    //k := 10.2 + (3.51E-2)*(1091.45-273.15);         
  end k_correlation;
  
  function rho_correlation
    input Real T_local "Temperature (K)";
    output Real rho "isobaric specific heat (J/kg.K)";
    protected Real T_celcius_local;
  algorithm
    T_celcius_local := T_local - 273.15;        
    rho := (17.15 - (8.63E-4)*(T_celcius_local + 20.))*1000.0; //in kg/m3    
    //rho := 17110; //Room temp value is actually 17320, 98.5% of theoretical 17580    
  end rho_correlation;
  
  // --- Heat Pipe / Boundary Parameters ---
  parameter Real T_HP_nominal = 1073.3074676166252 "Fixed HP temperature for steady-state [K]"; 
  parameter Integer n_hp = 8 "Number of heat pipes";  
  
  Real Q_loss "Heat lost through MLI [W]";
  Real P_integral; 
  
  parameter Real w_contact = 0.015 "Effective contact width per heat pipe [m]";  
  parameter Real A_hp_eff = n_hp * w_contact * L "Total effective HP contact area [m2]";
  Real q_gen[N] "Required Q''' [W/m3]";   
  Real q_gen_prof[N]; //Define q_gen based on a prescribed radial profile
  
  Real T[N] "Temperature at shell centers [K]";
  Real q_flow[N+1] "Heat flux across shell faces [W/m^2]";
  Real V_shell[N] "Volume of each shell [m3]";
  
  Real T_outer_wall "Temperature at the heat pipe interface [K]";

  //--Nodalisation
  parameter Real dr = (r_outer - r_inner) / N;
  final parameter Real r_face[N+1] = {r_inner + (i-1)*dr for i in 1:N+1};  
  
  Real power_profile_integral;
  Real verified_power_integral;

  output Real Q_evap_out "Heat flow out of heat pipe wall [W]";
  output Real T_mean "Volume-weighted mean temperature [K]";

  parameter Real recoverable_power_fraction = 0.93703; //Near-field heating fraction as per "KRUSTY Reactor Design" paper
  
  parameter Real Q_loss_nominal = 350.0 "nominal heat loss rate, set to agree with experiment"; 
  
  //Determine what nominal nuclear heating is needed to get the final effective thermal power that we expect 
  parameter Real P_total_nom = (2350 + Q_loss_nominal)/recoverable_power_fraction "Nominal power from fission and decay [W]";
  
  
initial equation
  for i in 1:N loop 
    der(T[i]) = 0;
  end for;
  
equation
  if Q_loss_input > 1e-6 then 
    Q_loss = Q_loss_input;
  else
    Q_loss = Q_loss_nominal;
  end if;    
  
  //P_integral is the actual thermal power to compute the temperature field without explicitly modelling thermal loss 
  if P_integral_input > 1e-6 then
    P_integral = P_integral_input*recoverable_power_fraction - Q_loss;
  else
    P_integral = P_total_nom*recoverable_power_fraction - Q_loss;
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
  
  
  q_flow[1] = 0; // Inner Boundary: Adiabatic (Control Rod Hole) 
  
  // The temperature at the evap wall is fixed, received from HP solver. 
  if T_outer_wall_input > 1e-6 then
    T_outer_wall = T_outer_wall_input;
  else 
    T_outer_wall = T_HP_nominal;
  end if;
  
  // Energy Balance for each shell in terms of heat flux
  // q_flow is positive when heat flows towards positive r (negative T gradient)
  for i in 2:N loop    
    q_flow[i] = -k_correlation(T[i]) * (T[i] - T[i-1]) / dr; //Note this only looks like a backward difference, but it's actually forward difference due to indexing change
  end for;  
    
  // Outer boundary
  q_flow[N+1] = -k_correlation(T_outer_wall) * (T_outer_wall - T[N]) / (dr/2); //Note half delta r since T_outer_wall doesn't correspond to a node average
    
  for i in 1:N-1 loop    
    q_gen[i] = q_gen_prof[i]/power_profile_integral*P_integral; //update q_gen    
    
    //derived radial discretisation with spatially uniform conductivity:
    //rho_correlation(T[i]) * cp_correlation(T[i]) * der(T[i]) = q_flow[i]*(1/dr - 1/r_face[i]) - q_flow[i+1]/dr + q_gen[i];    
    
    //derivation with spatially-dependent conductivity:
    rho_correlation(T[i]) * cp_correlation(T[i]) * der(T[i]) = q_flow[i]*(1/dr - 1/r_face[i] + 1/k_correlation(T[i])*(k_correlation(T[i+1])-k_correlation(T[i]))/dr ) - q_flow[i+1]/dr + q_gen[i];    
  end for;  
  // handle outer boundary case
  q_gen[N] = q_gen_prof[N]/power_profile_integral*P_integral; //update q_gen
    rho_correlation(T[N]) * cp_correlation(T[N]) * der(T[N]) = q_flow[N]*(1/dr - 1/r_face[N] + 1/k_correlation(T[N])*(k_correlation(T_outer_wall)-k_correlation(T[N]))/(dr/2) ) - q_flow[N+1]/dr + q_gen[N];    
  
    
  verified_power_integral = sum(q_gen[i] * V_shell[i] for i in 1:N);   
  T_mean = sum(T[i] * V_shell[i] for i in 1:N) / fuel_volume;    
  Q_evap_out = q_flow[N+1]*(2*pi*r_face[N+1]*L); //Not modelling heat loss explicitly (else it would be subtracted here)
  
  annotation(uses(Modelica(version="4.0.0")));
end conductionKRUSTYRadialDisc;
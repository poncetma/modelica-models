model SimpleTHKRUSTY_ana
/*
 Control-volume-based implementation of radial heat conduction in KRUSTY with a uniform source. 
 This version includes fixed thermal losses out of the core (heat not going through the heat pipes).
 
 Inputs: integral power, T_wall
 Outputs: HP wall heat flux, T_avg
*/
import Modelica.Constants.pi;
import Modelica.Math; 

  // --- Physical Parameters (KRUSTY U-10Mo) ---
  final constant Integer N = 200 "Number of radial shells"; //tested up to 500
  parameter Real r_inner = 0.02 "Inner fuel radius [m]";
  //parameter Real r_outer = 0.055 "Outer fuel radius [m]";
  parameter Real L = 0.25 "Core height [m]";
  //parameter Real fuel_volume = pi*(r_outer^2. - r_inner^2.)*L;
  
  //The fuel volume should  match exactly: 0.00186454 m3
  parameter Real fuel_volume = 0.00186454; //computed in OpenFOAM from KRUSTY mesh. 
  //with fixed inner radius, compute the appropriate equivalent outer radius
  parameter Real r_outer = sqrt(r_inner^2 + fuel_volume/L/pi);
  
  //parameter Real rho = 17110 "Density [kg/m3]";
  //Real cp[N];// = 189; //140 "Specific heat capacity [J/kg.K]";
  //parameter Real k = 37.5; //20 "Thermal conductivity [W/m.K]";
  
  //correlations for mat properties
  function cp_correlation
    input Real T_local "Temperature (K)";
    output Real cp "isobaric specific heat (J/kg.K)";
    protected Real T_celcius_local;
  algorithm
    T_celcius_local := T_local - 273.15;
    //cp := (0.137E-3 + 5.12E-5*T_celcius_local + 1.99E-8*T_celcius_local^2.)*1000; 
    //cp := (0.137E-3 + (5.12E-5)*(1091.45-273.15) + (1.99E-8)*(1091.45-273.15)^2.)*1000;  //this correlation gives the wrong values...    
    cp := (0.00007333333*(T_celcius_local) + 0.134666667)*1000; //own correlation from Burkes data 
    //cp := (0.00007333333*(1091.45-273.5) + 0.134666667)*1000; 
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
    /*
    if T_celcius_local > 0 then
      k := 10.2 + 3.51E-2*T_celcius_local; 
    else 
      k:= 37.5;
    end if;    
    */
  end k_correlation;
  
  function rho_correlation
    input Real T_local "Temperature (K)";
    output Real rho "isobaric specific heat (J/kg.K)";
    protected Real T_celcius_local;
  algorithm
    T_celcius_local := T_local - 273.15;        
    rho := (17.15 - (8.63E-4)*(T_celcius_local + 20.))*1000.0; //in kg/m3
    //rho := (17.15 - (8.63E-4)*((1091.45-273.25) + 20.))*1000.0; //in kg/m3
    //rho := 17110; //Room temp value is actually 17320, 98.5% of theoretical 17580
    
  end rho_correlation;
  

// --- Heat Pipe / Boundary Parameters ---
  parameter Real T_HP_nominal = 1073.3074676166252; //1073.5; //873.15 "Fixed HP temperature for steady-state [K]";
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
  Real q_gen[N] "Required Q''' [W/m3]";
  //Define q_gen based on a prescribed radial profile
  Real q_gen_prof[N];
  
  // --- Variables ---
  //Real T[N](start=fill(800, N)) "Temperature at shell centers [K]";
  Real T[N] "Temperature at shell centers [K]";
  Real Q_flow[N+1] "Heat flow rate across shell faces [W]";
  Real V_shell[N] "Volume of each shell [m3]";
  
  Real T_outer_wall "Temperature at the heat pipe interface [K]";
  output Real Q_evap_out "Heat flow out of heat pipe wall [W]";
  output Real T_mean "Volume-weighted mean temperature [K]";

  // --- Geometry Constants ---
  final parameter Real dr = (r_outer - r_inner) / N;
  //final parameter Real r_face[N+1] = {r_inner + (i-1)*dr for i in 1:N+1};                               
  final parameter Real r_face[N+1] = {sqrt(r_inner^2 + (i-1)*(r_outer^2 - r_inner^2)/N) for i in 1:N+1}; //equal volume spacing
  
  // Node placement at logarithmic mean radius
  //final parameter Real r_n[N] = {(r_face[i+1] - r_face[i]) / Math.log(r_face[i+1]/r_face[i]) for i in 1:N};
  // Node placement at rms
  final parameter Real r_n[N] = {sqrt((r_face[i]^2 + r_face[i+1]^2)/2) for i in 1:N};
  
  
  Real power_profile_integral;
  Real verified_power_integral;

  //Real HTC_loss; //model temp-dependent losses out of MLI? 
  Real R_hp;
  Real R_constriction;
  
  final parameter Real G_corr_int[N-1] = {
    pi * L * (r_face[i+1]^2 - ( (r_n[i+1]^2 - r_n[i]^2) / (2 * log(r_n[i+1]/r_n[i])) ) )
    for i in 1:N-1
  };
  // Note: For the boundary flux (i=N+1), use r_n[N] and r_f[N+1] (the wall)
  final parameter Real G_corr_wall = pi * L * (r_face[N+1]^2 - ( (r_face[N+1]^2 - r_n[N]^2) / (2 * log(r_face[N+1]/r_n[N]))));
  

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

  //P_integral = P_total_nom; 
  //P_integral = P_integral_input;   
  if P_integral_input > 1e-6 then
    P_integral = P_integral_input;
  else
    P_integral = P_total_nom + Q_loss;
  end if;
  
  //q_flux_hp = P_total_nom/A_hp_eff; 
  //q_flux_hp = Q_hp_input/A_hp_eff;  

  for i in 1:N loop
    //q_gen[i] = P_integral / fuel_volume; //uniform heating
    //follow radial profile and renormalise to nominal integral power
    if i < 9./10.*N then
      q_gen_prof[i] = 1;
    else
      q_gen_prof[i] = 1; //1 reverts back to normal profile   
    end if;
    V_shell[i] = pi * (r_face[i+1]^2 - r_face[i]^2) * L;
  end for;   
  power_profile_integral = sum(q_gen_prof[i] * V_shell[i] for i in 1:N);

  

  // 1. Heat Flows at Faces
  Q_flow[1] = 0; // Inner Boundary: Adiabatic (Control Rod Hole) 
  
  // 1. Boundary Condition Logic
  // The temperature at the evap wall is fixed, received from HP solver. 
  
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
  //Now using logarithmic resistance to model curvature and achieve O(N^2) error
  for i in 2:N loop
    //Q_flow[i] = -k_correlation((T[i-1] + T[i])/2.) * (2 * pi * r_face[i] * L) * (T[i] - T[i-1]) / dr; 
    
    //Q_flow[i] = (T[i-1] - T[i]) / (Math.log(r_n[i]/r_n[i-1]) / (2 * pi * k_correlation((T[i-1]+T[i])/2) * L));
    //Q_flow[i] = (T[i-1] - T[i]) / (Math.log(r_n[i]/r_n[i-1]) / (2 * pi * k_correlation((T[i-1]+T[i])/2) * L));
    Q_flow[i] = (T[i-1] - T[i]) / (Math.log(r_n[i]/r_n[i-1]) / (2 * pi * k_correlation((T[i-1]+T[i])/2) * L)) + P_integral/fuel_volume * G_corr_int[i-1];
  end for;
  //Q_flow[N+1] = -k * (2 * pi * r_face[N+1] * L) * (T_outer_wall - T[N]) / (dr/2); //half-dr due to cell-centredness
  //T_outer_wall = T[N] - (Q_flow[N+1] / (2 * pi * r_face[N+1] * L)) * (dr/2) / k;
  //Q_flow[N+1] = -k * (A_hp_eff) * (T_outer_wall - T[N]) / (dr/2); //half-dr due to cell-centredness
  //Now there are two parallel heat flow paths at the outer wall. 
  //Assume Q_loss is constant (true for small variation of outer wall temp)
  //Q_flow[N+1] = -k_correlation((T[N]+T_outer_wall)/2.) * (A_hp_eff) * (T_outer_wall - T[N]) / (dr/2) + Q_loss; 
  
  //R_hp = (Math.log(r_outer/r_n[N]) * r_outer) / (k_correlation((T[N]+T_outer_wall)/2) * A_hp_eff);
  R_hp = (Math.log(r_outer/r_n[N]) * r_outer) / (k_correlation((T[N]+T_outer_wall)/2) * A_hp_eff);
  R_constriction = (Math.log(r_outer/r_n[N]) / (2*pi*k_correlation((T[N]+T_outer_wall)/2)*L)) * ((2*pi*r_outer*L) / A_hp_eff);
  //Q_flow[N+1] = (T[N] - T_outer_wall) / R_hp + Q_loss;
  //T[N] - T_outer_wall = ((P_integral/fuel_volume)/(4*k_correlation((T[N] + T_outer_wall)/2.))) * (r_outer^2 - r_n[N]^2)  + (Q_flow[N+1] - P_integral) * (Math.log(r_outer/r_n[N])/(2*pi*k_correlation((T[N] + T_outer_wall)/2.)*L)) + Q_flow[N+1] * (R_constriction - (Math.log(r_outer/r_n[N])/(2*pi*k_correlation((T[N] + T_outer_wall)/2.)*L)));
  Q_flow[N+1] = (T[N] - T_outer_wall) / (Math.log(r_face[N+1]/r_n[N]) / (2*pi*k_correlation((T[N]+T_outer_wall)/2)*L)) + P_integral/fuel_volume * G_corr_wall;
  
  for i in 1:N loop    
    q_gen[i] = q_gen_prof[i]/power_profile_integral*P_integral; //update q_gen
    rho_correlation(T[i]) * V_shell[i] * cp_correlation(T[i]) * der(T[i]) = Q_flow[i] - Q_flow[i+1] + q_gen[i] * V_shell[i];    
  end for;  
    verified_power_integral = sum(q_gen[i] * V_shell[i] for i in 1:N);
 

  T_mean = sum(T[i] * V_shell[i] for i in 1:N) / fuel_volume;
  Q_evap_out = Q_flow[N+1] - Q_loss;
  
  annotation(uses(Modelica(version="4.0.0")));
end SimpleTHKRUSTY_ana;
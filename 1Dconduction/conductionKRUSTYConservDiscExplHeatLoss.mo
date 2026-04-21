model conductionKRUSTYConservDiscExplHeatLoss
  /*
   Conservative, finite-volume implementation of 1D heat conduction in KRUSTY with a possibly non-uniform source. 
      
   To enable co-simulation as an FMU with a less robust CVODE setup, a semi-implicit/"lagged" scheme is used for the thermal conductivity: it is updated frequently but not at every time-step, and it is not solved simultaneously. 
   
   This version explicitly accounts thermal losses out of the core with a heat transfer coefficient tuned to match
   the nominal observed heat loss. The heat loss is no longer an input.
   
   Inputs: total integral power (heat deposition rate from fission + decay heat), evaporator wall temperature 
   Outputs: HP wall heat flux, Core average temperature
  */
  import Modelica.Constants.pi;
  input Real P_integral_input "Instantaneous integral power [W]";
  input Real T_outer_wall_input "Temperature of heat pipe wall [K]";
  input Real k_input "Conductivity (mean), controlled externally [W/m/K]"; 
  
  // --- Physical Parameters (KRUSTY U-10Mo) ---
  final constant Integer N = 25 "Number of radial shells";
  //tested up to several hundred nodes
  parameter Real r_inner = 0.02 "Inner fuel radius [m]";
  parameter Real L = 0.25 "Core height [m]";
  parameter Real fuel_volume = 0.00186454;
  //computed in OpenFOAM from KRUSTY mesh.
  parameter Real r_outer = sqrt(r_inner^2 + fuel_volume/L/pi);
  //with fixed inner radius, compute the appropriate equivalent outer radius
  //correlations for mat properties

  function cp_correlation
    input Real T_local "Temperature (K)";
    output Real cp "isobaric specific heat (J/kg.K)";
  protected
    Real T_celcius_local;
  algorithm
    T_celcius_local := T_local - 273.15;
    cp := (0.00007333333*(T_celcius_local) + 0.134666667)*1000;
//own correlation from Burkes data
//cp := 189;
  end cp_correlation;

  function k_correlation
    input Real T_local "Temperature (K)";
    output Real k "conductivity (W/m/K)";
  protected
    Real T_celcius_local;
  algorithm
    T_celcius_local := T_local - 273.15;
    k := 10.2 + (3.51E-2)*T_celcius_local;
    //k := 10.2 + (3.51E-2)*(800);
  end k_correlation;

  function rho_correlation
    input Real T_local "Temperature (K)";
    output Real rho "density (kg/m3)";
  protected
    Real T_celcius_local;
  algorithm
    T_celcius_local := T_local - 273.15;
    rho := (17.15 - (8.63E-4)*(T_celcius_local + 20.))*1000.0;
//in kg/m3
//rho := 17110; //Room temp value is actually 17320, 98.5% of theoretical 17580
  end rho_correlation;

  // --- Heat Pipe / Boundary Parameters ---
  parameter Real T_HP_nominal = 1073.15 "Fixed HP temperature for steady-state [K]";
  parameter Integer n_hp = 8 "Number of heat pipes";
  Real Q_loss "Heat lost through MLI [W], to be computed";
  Real P_integral;
  /*
    parameter Real w_contact = 0.015 "Effective contact width per heat pipe [m]";  
    parameter Real A_hp_eff = n_hp * w_contact * L "Total effective HP contact area [m2]";
    */
  Real q_gen[N] "Required Q''' [W/m3]";
  Real q_gen_profile[N];
  //Define q_gen based on a prescribed radial profile
  Real T[N](each start = 1073.15) "Temperature at shell centers [K]";
  Real T_lagged[N] (each start = 1073.15); //Introduce a lagged temperature for semi-implicit scheme
  
  Real q_flow[N + 1] "Heat flux across shell faces [W/m^2]";
  Real V_shell[N] "Volume of each shell [m3]";
  Real T_outer_wall "Temperature at the heat pipe interface [K]";
  Real T_outer_wall_lagged;
  //--Nodalisation
  parameter Real dr = (r_outer - r_inner)/N;
  parameter Real r_face[N + 1] = {r_inner + (i - 1)*dr for i in 1:N + 1};
  Real power_profile_integral;
  Real verified_power_integral;
  output Real Q_evap_out "Heat flow out of heat pipe wall [W]";
  output Real T_mean "Volume-weighted mean temperature [K]";
  parameter Real recoverable_power_fraction = 0.93703;
  //Near-field heating fraction as per "KRUSTY Reactor Design" paper
  parameter Real Q_loss_nominal = 350.0 "nominal heat loss rate through insulation + the rest of the surrounding components, set to agree with experiment";
  parameter Real T_outer_wall_nominal = T_HP_nominal;
  parameter Real T_ambient = 300.0 "ambient temperature";
  parameter Real HTC_loss = Q_loss_nominal/(T_outer_wall_nominal - T_ambient) "tuned HTC that gives the nominal heat loss [W/K]";
  parameter Real outer_wall_area = 2*pi*r_outer*L;
  Real Q_outerwall_out;
  //Determine what nominal nuclear heating (fission + decay power) is needed to get the final effective thermal power that we expect at nominal conditions [2350 W after heat losses]
  parameter Real P_total_nom = (2350 + Q_loss_nominal)/recoverable_power_fraction "Nominal power from fission and decay [W]";
  
  parameter Real dt_lag = 10.0 "frequency with which to update the conductivity within a given call to this Modelica model/FMU";    
  Real k_mean "mean conductivity, to ensure it's updating correctly";
  
initial equation
  for i in 1:N loop
    der(T[i]) = 0;
  end for;
equation

  //when sample(0, dt_lag) then
  when {initial()} then //trigger at the beginning and every dt_lag interval. dt_lag is chosen to be much greater than the typical FMU timestep so that it won't be called internally (encountered issues with that).
    T_lagged = pre(T);
    T_outer_wall_lagged = pre(T_outer_wall);
  end when;

  Q_loss = HTC_loss*(T_outer_wall - T_ambient); //At zero power, this should drive the temperature down to T_ambient
//  Q_loss = 350.0;
//P_integral is the actual thermal power to compute the temperature field while accounting for heat loss
  if P_integral_input > 1E-9 then
    P_integral = P_integral_input*recoverable_power_fraction; //no longer subtracting Q_loss here
  else
    P_integral = P_total_nom*recoverable_power_fraction;
  end if;
  
  for i in 1:N loop
//follow radial profile and renormalise to nominal integral power
    if i < 9./10.*N then
      q_gen_profile[i] = 1;
    else
      q_gen_profile[i] = 1;
//1 reverts back to normal profile
    end if;
    V_shell[i] = pi*(r_face[i + 1]^2 - r_face[i]^2)*L;
  end for;
  power_profile_integral = sum(q_gen_profile[i]*V_shell[i] for i in 1:N);
  q_flow[1] = 0;
// Inner Boundary: Adiabatic (Control Rod Hole)
// The temperature at the evap wall is fixed, received from HP solver.
  if T_outer_wall_input > 1E-9 then
    T_outer_wall = T_outer_wall_input;
  else
    T_outer_wall = T_HP_nominal;
  end if;
  
  if k_input > 1e-9 then
    k_mean = k_input;
    
    // Energy Balance for each shell in terms of heat flux
    // q_flow is positive when heat flows towards positive r (negative T gradient)
    for i in 2:N loop
      //use harmonic interpolation for k
      q_flow[i] = -1.*k_mean*(T[i] - T[i - 1])/dr;
      //Note this only looks like a backward difference, but it's actually forward difference due to indexing shift
    end for;
    // Outer boundary: heat flux based on half-delta-r to boundary and computed heat loss
    q_flow[N + 1] = -1.*k_mean*(T_outer_wall - T[N])/(dr/2) + Q_loss/outer_wall_area;
  else 
    k_mean = k_correlation(T_mean);
    
    // Energy Balance for each shell in terms of heat flux
    // q_flow is positive when heat flows towards positive r (negative T gradient)
    for i in 2:N loop
      //use harmonic interpolation for k
      q_flow[i] = -1.*(2/(1/k_correlation(T_lagged[i - 1]) + 1/k_correlation(T_lagged[i])))*(T[i] - T[i - 1])/dr;
      //Note this only looks like a backward difference, but it's actually forward difference due to indexing shift
    end for;
    // Outer boundary: heat flux based on half-delta-r to boundary and computed heat loss
    q_flow[N + 1] = -1.*(2/(1/k_correlation(T_lagged[N]) + 1/k_correlation(T_outer_wall_lagged)))*(T_outer_wall - T[N])/(dr/2) + Q_loss/outer_wall_area;
  end if;   
  

// Update heat source
  for i in 1:N - 1 loop
    q_gen[i] = q_gen_profile[i]/power_profile_integral*P_integral;
    rho_correlation(T[i])*cp_correlation(T[i])*der(T[i]) = ((2*pi*r_face[i]*L)*q_flow[i] - (2*pi*r_face[i + 1]*L)*q_flow[i + 1])/V_shell[i] + q_gen[i];
  end for;
// handle outer boundary case
  q_gen[N] = q_gen_profile[N]/power_profile_integral*P_integral; //update q_gen
  rho_correlation(T[N])*cp_correlation(T[N])*der(T[N]) = ((2*pi*r_face[N]*L)*q_flow[N] - (2*pi*r_face[N + 1]*L)*q_flow[N + 1])/V_shell[N] + q_gen[N];
  
  verified_power_integral = sum(q_gen[i]*V_shell[i] for i in 1:N);
  T_mean = sum(T[i]*V_shell[i] for i in 1:N)/fuel_volume;
  
  
  Q_outerwall_out = q_flow[N + 1]*(2*pi*r_outer*L);
//Integrated heat flux out of outer wall
  Q_evap_out = Q_outerwall_out - Q_loss;
//Outer wall heat flux minus heat loss contribution
  annotation(
    uses(Modelica(version = "4.0.0")),
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end conductionKRUSTYConservDiscExplHeatLoss;
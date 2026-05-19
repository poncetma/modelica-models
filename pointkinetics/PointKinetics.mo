model PointKinetics
"
A compact implemention of the point kinetics equations in their power form, 
Ref: https://www.osti.gov/servlets/purl/1484139
This model has been verified against an old scipy-based simple PK solver. 

Also added the ability to choose a power setpoint and insert external reactivity as needed to maintain that setpoint (as done in KRUSTY startup).

Now made exportable as an FMU. 
Inputs: Fuel temperature (T_Fuel), external reactivity, desired reactivity coefficient and reference temperature
Outputs: Power (P)
"
//forced to set these as inputs in order to control them within the FMU
input Real T_fuel_ref_input;
input Real T_fuel_input "Average fuel temperature, obtained from external solver";
input Real alpha_Tf_input; 
input Real rho_ext_input "external reactivity, outside input";
input Real P_0_input "initial fission power [W]";
Real P_0;
input Real P_setpoint_input "desired setpoint for power [W]";
input Boolean ACTIVE_POWER_CONTROL "Choose whether or not to actively control power using PI controller"; 

parameter Real rho_0 = 0.0e-5 "initial reactivity [pcm]"; 
parameter Real delta_rho = 9e-5 "amount to increase external reactivity by [pcm]"; 
parameter Real rho_max = 0.299*Beta; //~30 cent max

Real rho_ext ;//(start = 0, fixed = true); 
Real rho_cmd; // (start = 0, fixed = true);
Real T_fuel; 

output Real P "instantaneous fission power"; 
Real Cs[6] "instaneous group-wise precursor concentration"; 
Real rho "instantaneous net reactivity";
Real rho_fb "reactivity feedback";
Real P_setpoint; 
/*
Real P_set_latch (start=0, fixed=true); 
Boolean P_latch_alt (start=false, fixed=true); 
discrete Real rho_ext_input_d(start=0, fixed=true);
*/
Real P_max ;//(start=0, fixed=true);

parameter Real alpha_Tf_default = -0.1844*0.01*Beta "fuel TRC, Poston et al"; 
parameter Real T_fuel_ref_default = 1090 "reference temperature, from steady-state TH solve [K]"; //1091.173; //
Real alpha_Tf;
Real T_fuel_ref;

//This correlation gives the integral reactivity feedback (relative to *operating* (nominal) temp).
  //From "KRUSTY Design and Modelling" slide 88

function alpha_Tf_poly 
  input Real T_fuel "Temperature (K)";
  output Real alpha_Tf "Fuel temperature-reactivity coefficient (pcm/K)";
  algorithm
    alpha_Tf := (-1.6951E-11*T_fuel^3. + 5.0121E-8*T_fuel^2 - 1.4888E-4*T_fuel -7.9756E-2)*0.01*Beta; //in pcm/K
end alpha_Tf_poly;

//Values from KRUSTY papers - need to be parameters to be retrievable in Python
constant Real Beta = 0.00688; //from Stolte et al.
parameter Real betas[6] = {0.037, 0.211, 0.187, 0.407, 0.131, 0.027}*Beta; //from Grove et al.
parameter Real lambdas[6] = {0.01273, 0.03175, 0.116, 0.3118, 1.399, 3.876}; //from Grove et al.
constant Real Lambda = 5.20395e-6; //from Stolte et al.

//Values from MOOSE VTB (Serpent outputs)
/*
constant Real Beta = sum(betas); 
constant Real betas[6] = {23.06, 118.77, 115.32, 324.70, 96.07, 33.82}*1e-5;
constant Real lambdas[6] = {0.0125, 0.0318, 0.1094, 0.3173, 1.3529, 8.6655};
constant Real Lambda = 5.31603e-6; 
*/
parameter Real nu_bar = 2.54218 "mean fission neutron multiplicity";
parameter Real Sigma_f = 1e-2; //Don't know this value or n_speed a-priori since Serpent results are not condensed
parameter Real n_speed = 1e6 "one-group neutron speed";
//parameter Real Lambda = 5.5e-8; //seems very low, but set to match the Wang/Jiang paper. Little impact on gradual transients.

constant Real w_f = 200.0e6*1.602e-19 "Energy release per fission [J]";
  
//Conventional formulation:
Real n "normalised neutron density";
Real Cs_alt[6] "instaneous group-wise precursor concentration"; 
output Real P_alt;


//PI controller variables
input Real K_p_input;
input Real K_i_input;
Real K_p; 
Real K_i;
Boolean PI_active (start = false, fixed=true);
Real ctrl_out ;
Real ctrl_error ;
Real ctrl_error_integral (start = 0, fixed=true);

initial equation
//Power-formulation:
P = P_0;
nu_bar*betas/w_f*P = lambdas.*Cs; 

//Conventional formulation:
n = 1;
betas/Lambda*n = lambdas.*Cs_alt;
/*
P_max = P_alt;
ctrl_out = 0;
*/
equation

if P_0_input > 1e-12 then 
  P_0 = P_0_input;
else 
  P_0 = 1e-2;
end if;

if P_setpoint_input > 1e-9 then 
  P_setpoint = P_setpoint_input;
else
  P_setpoint = 3000;
end if;


//PI controller (proven approach)
if K_p_input > 1E-9 and K_i_input > 1E-9 then 
  K_p = K_p_input;
  K_i = K_i_input;
else
  K_p = 0.2; //Proportional gain
  K_i = K_p/100.; //Integral gain
end if; 

ctrl_error = if PI_active then (P_setpoint - P_alt) else 0;
der(ctrl_error_integral) = ctrl_error;
ctrl_out = K_p*ctrl_error + K_i*ctrl_error_integral;


P_max = max(P_max, P_alt); //update the maximum power yet reached

//For running the FMU, need to have a tolerance for these 'less than' conditions since they may otherwise be triggered prematurely.
when (not pre(PI_active)) and (P_alt - P_setpoint < -10) and (P_alt - P_max < -10) then 
      PI_active = true; //only activate when not previously set and past the setpoint and past previous max
      /*
      Modelica.Utilities.Streams.print("ACTIVATED PI CTRL");
      Modelica.Utilities.Streams.print("time: " + String(time));
      Modelica.Utilities.Streams.print("P_alt: " + String(P_alt));
      Modelica.Utilities.Streams.print("P_setpoint: " + String(P_setpoint));
      Modelica.Utilities.Streams.print("P_max: " + String(P_max));
      */
end when;   



if (ACTIVE_POWER_CONTROL) then 
  rho_ext = min(rho_max, ctrl_out + rho_ext_input);
else 
  rho_ext = rho_ext_input; 
end if;

rho_cmd =0;



if T_fuel_ref_input > 1e-9 then
  T_fuel_ref = T_fuel_ref_input;
else
  T_fuel_ref = T_fuel_ref_default;
end if;

if T_fuel_input > 1e-9 then
  T_fuel = T_fuel_input;
else
  T_fuel = T_fuel_ref; 
end if;

if alpha_Tf_input > 1e-9 then
    alpha_Tf = alpha_Tf_input;
else 
    alpha_Tf = alpha_Tf_poly(T_fuel);
end if;

rho_fb = alpha_Tf*(T_fuel - T_fuel_ref); 
rho = rho_0 + rho_ext + rho_fb; 

//Power formulation:
der(P) = (rho - Beta)/Lambda.*P + w_f/Lambda/nu_bar*sum(lambdas.*Cs);
der(Cs) = nu_bar*betas/w_f*P - lambdas.*Cs; //precursors all together

//Conventional formulation:
der(n) = (rho - Beta)/Lambda.*n + sum(lambdas.*Cs_alt);
der(Cs_alt) = betas/Lambda*n - lambdas.*Cs_alt;
//P_alt = n*P_0;
P_alt = max(1e-8, n*P_0);

annotation(
    __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian",
    __OpenModelica_simulationFlags(lv = "LOG_STDOUT,LOG_ASSERT,LOG_STATS", s = "dassl", variableFilter = ".*"));
end PointKinetics;
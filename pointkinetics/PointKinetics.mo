model PointKinetics
"
A compact implemention of the point kinetics equations in their power form, 
Ref: https://www.osti.gov/servlets/purl/1484139
This model has been verified against an old scipy-based simple PK solver. 

Now made exportable as an FMU. 
Inputs: Fuel temperature (T_Fuel), external reactivity
Outputs: Power (P)
"
Real rho_ext; 
Real T_fuel; 
input Real rho_ext_input "external reactivity, outside input";
input Real T_fuel_input "Average fuel temperature, obtained from external solver";
output Real P "instantaneous fission power"; 
Real Cs[6] "instaneous group-wise precursor concentration"; 
Real rho "instantaneous net reactivity";
Real rho_fb "reactivity feedback";

//forced to set these as inputs in order to control them within the FMU
input Real alpha_Tf_input; 
input Real T_fuel_ref_input;
parameter Real alpha_Tf_default = -0.1844*0.01*Beta "fuel TRC, Poston et al"; 
parameter Real T_fuel_ref_default = 1090 "reference temperature, from steady-state TH solve [K]"; //1091.173; //
Real alpha_Tf;
Real T_fuel_ref;

//This correlation gives the integral reactivity feedback (relative to *operating* (nominal) temp). 
//From "KRUSTY Design and Modelling" slide 88
function alpha_Tf_poly 
  input Real T_fuel "Temperature (K)";
  output Real alpha_Tf "isobaric specific heat (J/kg.K)";
  algorithm
    alpha_Tf := (-1.6951E-11*T_fuel^3. + 5.0121E-8*T_fuel^2 - 1.4888E-4*T_fuel -7.9756E-2)*0.01*Beta; //in pcm/K
end alpha_Tf_poly;

constant Real Beta = 0.00688; //value directly from KRUSTY papers
constant Real betas[6] = {0.037, 0.211, 0.187, 0.407, 0.131, 0.027}*Beta;
parameter Real lambdas[6] = {0.01273, 0.03175, 0.116, 0.3118, 1.399, 3.876}; //value from KRUSTY papers

parameter Real nu_bar = 2.54218 "mean fission neutron multiplicity";
parameter Real Sigma_f = 1e-2; //Don't know this value or n_speed a-priori since Serpent results are not condensed
parameter Real n_speed = 1e6 "one-group neutron speed";
parameter Real Lambda = 5.5e-8; //seems very low, but set to match the Wang/Jiang paper. Little impact on gradual transients.
constant Real w_f = 200.0e6*1.602e-19 "Energy release per fission [J]";

input Real P_0 = 5000.0 "initial fission power [W]";
parameter Real rho_0 = 0.0 "initial reactivity"; 
  
//Conventional formulation:
Real n "normalised neutron density";
Real Cs_alt[6] "instaneous group-wise precursor concentration"; 
output Real P_alt;

initial equation
//Power-formulation:
P = P_0;
nu_bar*betas/w_f*P = lambdas.*Cs; 

//Conventional formulation:
n = 1;
betas/Lambda*n = lambdas.*Cs_alt;

equation
if T_fuel_ref_input > 1e-6 then
  T_fuel_ref = T_fuel_ref_input;
else
  T_fuel_ref = T_fuel_ref_default;
end if;

if T_fuel_input > 1e-6 then
  T_fuel = T_fuel_input;
  alpha_Tf = alpha_Tf_poly(T_fuel);
else
  T_fuel = T_fuel_ref; 
  alpha_Tf = alpha_Tf_poly(T_fuel_ref);
end if;

rho_ext = rho_ext_input;

rho_fb = alpha_Tf*(T_fuel - T_fuel_ref); 
rho = rho_0 + rho_ext + rho_fb; 

//Power formulation:
der(P) = (rho - Beta)/Lambda.*P + w_f/Lambda/nu_bar*sum(lambdas.*Cs);
der(Cs) = nu_bar*betas/w_f*P - lambdas.*Cs; //precursors all together

//Conventional formulation:
der(n) = (rho - Beta)/Lambda.*n + sum(lambdas.*Cs_alt);
der(Cs_alt) = betas/Lambda*n - lambdas.*Cs_alt;
P_alt = n*P_0;

end PointKinetics;
model PointKineticsWithDecayHeat
"
Point kinetics with decay heat tracking, based on effective precursor groupings (similar to DNPs).
This is the same general approach used in FRINK, though the exact implementation there is not published. 
This model uses the TRACE implementation (a power formulation, see TRACE V5 theory manual) 
with input data suited for KRUSTY (fast-spectrum, U238 fission products) from the ANS Standard for Decay Heat (1979).

Inputs: Fuel temperature (T_Fuel), external reactivity
Outputs: Total thermal power (P) 
"
input Real T_fuel_ref_input;
input Real T_fuel_input "Average fuel temperature, obtained from external solver";
input Real alpha_Tf_input; 
input Real rho_ext_input "external reactivity, outside input";
input Real P_0 "initial thermal power [W]";

parameter Real alpha_Tf_default = -0.1844*0.01*Beta "fuel TRC, Poston et al"; 
parameter Real T_fuel_ref_default = 1090 "reference temperature, from steady-state TH solve [K]"; //1091.173; //

parameter Real Lambda = 5.5e-8; //seems very low, but set to match the Wang/Jiang paper. Little impact on gradual transients.
constant Real Beta = 0.00688; //value directly from KRUSTY papers
parameter Real betas[6] = {0.037, 0.211, 0.187, 0.407, 0.131, 0.027}*Beta;
parameter Real lambdas[6] = {0.01273, 0.03175, 0.116, 0.3118, 1.399, 3.876}; //values from KRUSTY papers

parameter Real E_f = 195.0 "Energy release per fission, U238"; 
parameter Real lambdas_dh[23] = {
3.29E+00,
9.38E-01,
3.71E-01,
1.11E-01,
3.61E-02,
1.33E-02,
5.01E-03,
1.37E-03,
5.52E-04,
1.79E-04,
4.90E-05,
1.71E-05,
7.05E-06,
2.32E-06,
6.45E-07,
1.26E-07,
2.55E-08,
8.48E-09,
7.51E-10,
2.42E-10,
2.27E-13,
9.05E-14,
5.61E-15
} "Decay heat group decay constants [s^-1]";
parameter Real EDs[23] = {
1.23E+00,
1.15E+00,
7.07E-01,
2.52E-01,
7.19E-02,
2.83E-02,
6.84E-03,
1.23E-03,
6.84E-04,
1.70E-04,
2.42E-05,
6.64E-06,
1.01E-06,
4.99E-07,
1.64E-07,
2.34E-08,
2.81E-09,
3.62E-11,
6.46E-11,
4.50E-14,
3.67E-16,
5.63E-17,
7.16E-17
} "Decay power release per fission [MeV/s]";
parameter Real Es[23] = EDs./lambdas_dh./E_f "Decay power fraction [-]";

parameter Real rho_0 = 0.0 "initial reactivity"; 

//This correlation gives the integral reactivity feedback (relative to *operating* (nominal) temp). 
//From "KRUSTY Design and Modelling" slide 88
function alpha_Tf_poly 
  input Real T_fuel;
  output Real alpha_Tf;
  algorithm
    alpha_Tf := (-1.6951E-11*T_fuel^3. + 5.0121E-8*T_fuel^2 - 1.4888E-4*T_fuel -7.9756E-2)*0.01*Beta; //in pcm/K
end alpha_Tf_poly;

Real rho_ext; 
Real T_fuel; 
Real rho "instantaneous net reactivity";
Real rho_fb "reactivity feedback";
Real alpha_Tf;
Real T_fuel_ref;

Real Cs[6] "instaneous group-wise precursor power [W] "; 
Real Hs[23] "decay heat precursor energy [J]";
Real decay_heat_fraction;

output Real P_fiss "instantaneous fission power"; 
output Real P_dec "instantaneous decay heat power";
output Real P_tot "total instantaneous thermal power";

initial equation
if P_0 > 0 then
  //P = P_0;
  P_tot = P_0;
else
  //P = 5000.0;
  P_tot = 5000.0;
end if; 

betas/Lambda*P_fiss = lambdas.*Cs;

Es*P_fiss = lambdas_dh.*Hs "Inititial condition for infinite irradiation";

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

if rho_ext_input > 0 then
  rho_ext = rho_ext_input;
else 
  //rho_ext = 0.0; 
  if time < 100 then
    rho_ext= min(100, time)*1e-5 ;
  else
    rho_ext = 0.0;
  end if;
  
end if;

rho_fb = alpha_Tf*(T_fuel - T_fuel_ref); 
rho = rho_0 + rho_ext + rho_fb; 

der(P_fiss) = (rho - Beta)/Lambda.*P_fiss + sum(lambdas.*Cs);
der(Cs) = betas/Lambda*P_fiss - lambdas.*Cs;

der(Hs) = Es*P_fiss - lambdas_dh.*Hs;
P_dec = sum(lambdas_dh.*Hs);
P_tot = P_fiss + P_dec;
decay_heat_fraction = P_dec/P_tot;




end PointKineticsWithDecayHeat;
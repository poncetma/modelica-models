/*
A Stirling engine model which uses the "FTB" engine correlation (similar to KRUSTY's ASC) in tandem with a single lumped thermal mass.

Inputs: Heat pipe condenser temperature, taken as T_hot 
Outputs: Lumped stirling temperature (T_stirlin) 
*/
model StirlingOneNodeWithCorrelation

parameter Real C_cond_wall = 17.780965 "Heat pipe condenser wall capacitance [J/K]";
parameter Real C_stirling = 15*C_cond_wall "Stirling capacitance [J/K]"; 

parameter Real Q_draw_nominal = 2350/8 "nominal power draw [W]";
parameter Real T_stirling_activation = 650 + 273.15 "Temperature at which Stirlings are turned on in start-up run of KRUSTY"; 
parameter Real T_stirling_nominal = 650 + 273.15 "Stirling hot-side temperature, Poston et al., Fig. 10 [K]"; //630 + 273.15
parameter Real T_stirling_cold_nominal  = 65 + 273.15 "Cold sink temperature [K]"; 
parameter Real HTC_cold = Q_draw_nominal/(T_stirling_nominal - T_stirling_cold_nominal);

input Real Q_cs_input "incoming heat transfer rate from the heat pipe condenser [W]";
output Real T_stirling; 

Real Q_cs "instantaneous heat transfer rate from the heat pipe condenser to the Stirling engine [W]";
Real Q_stirling_bc "instantaneous heat transfer from the Stirling engine to the cold sink"; 

input Boolean COLD_START; 

initial equation

if COLD_START then 
  T_stirling = 15 + 273.15;
else 
  T_stirling = T_stirling_nominal;
end if;   

equation

//Heat flows
/*
if Q_cs_input > 1E-9 then 
  Q_cs = Q_cs_input;
else 
  Q_cs = Q_draw_nominal;
end if;
*/

//Receive the hot-side temperature as the external boundary condition and use that to compute Q_cs from the FTB correlation. 

Q_stirling_bc = HTC_cold*(T_stirling - T_stirling_cold_nominal);

//Energy balance
C_stirling*der(T_stirling) = Q_cs - Q_stirling_bc;


end StirlingOneNodeWithCorrelation;
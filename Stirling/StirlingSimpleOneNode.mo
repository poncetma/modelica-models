/*
A simple lumped capacitance model of the Stirling engine with a single node. Intended to be a modular version of the one built into the HeatPipeWithVapourCore model. 

Inputs: Heat transfer rate from heat pipe condenser (Q_cs_input)
Outputs: Lumped stirling temperature (T_stirling) 
*/
model StirlingSimpleOneNode

parameter Real C_cond_wall = 17.780965 "Heat pipe condenser wall capacitance [J/K]";
parameter Real C_stirling = 15*C_cond_wall "Stirling capacitance [J/K]"; 

parameter Real Q_draw_nominal = 2350/8 "nominal power draw [W]";
parameter Real T_stirling_activation = 650 + 273.15 "Temperature at which Stirlings are turned on in start-up run of KRUSTY"; 
parameter Real T_stirling_nominal = 650 + 273.15 "Stirling hot-side temperature, Poston et al., Fig. 10 [K]"; //630 + 273.15
parameter Real T_stirling_cold_nominal  = 65 + 273.15 "Cold sink temperature [K]"; 
parameter Real HTC_cold = Q_draw_nominal/(T_stirling_nominal - T_stirling_cold_nominal);

input Boolean COLD_START; 
input Real Q_cs_input "incoming heat transfer rate from the heat pipe condenser [W]";
output Real T_s; 

Real Q_cs "instantaneous heat transfer rate from the heat pipe condenser to the Stirling engine [W]";
Real Q_bc "instantaneous heat transfer from the Stirling engine to the cold sink"; 



parameter Boolean MODEL_STIRLING_ACTIVATION = true; 


Boolean STIRLING_ACTIVATED(start=false,fixed=true); //activated with 'when' statement

initial equation

if COLD_START then 
  T_s = 15 + 273.15;
else 
  T_s = 15 + 273.15; //T_stirling_nominal;
end if;   

equation

if MODEL_STIRLING_ACTIVATION then 
  when T_s > T_stirling_activation then 
    STIRLING_ACTIVATED = true;
    Modelica.Utilities.Streams.print("Stirling engine activated!");
  end when;
else 
  STIRLING_ACTIVATED = true;
end if; 

//Heat flows
if Q_cs_input > 1E-9 then 
  Q_cs = Q_cs_input;
else 
  Q_cs = Q_draw_nominal;
end if;

if STIRLING_ACTIVATED then 
  Q_bc = HTC_cold*(T_s - T_stirling_cold_nominal);
else
  Q_bc = 0;
end if;

//Energy balance
C_stirling*der(T_s) = Q_cs - Q_bc;


end StirlingSimpleOneNode;
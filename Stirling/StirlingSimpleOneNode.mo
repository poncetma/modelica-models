/*
A simple lumped capacitance model of the Stirling engine with a single node. Intended to be a modular version of the one built into the HeatPipeWithVapourCore model. 

Inputs: Heat transfer rate from heat pipe condenser (Q_cs_input)
Outputs: Lumped stirling temperature (T_s) 
*/
model StirlingSimpleOneNode

parameter Real C_adiabcond_wall = 17.780965 "Heat pipe adiabatic + condenser wall capacitance [J/K]";
parameter Real C_stirling = 15*C_adiabcond_wall "Stirling capacitance [J/K]"; 

parameter Real T_stirling_activation = 650 + 273.15 "Temperature at which Stirlings are turned on in start-up run of KRUSTY"; 
parameter Real T_stirling_nominal = 650 + 273.15 "Stirling hot-side temperature, Poston et al., Fig. 10 [K]"; //630 + 273.15
parameter Real T_stirling_cold_nominal  = 65 + 273.15 "Cold sink temperature [K]"; 


input Boolean COLD_START; 
input Real Q_cs_input "incoming heat transfer rate from the heat pipe condenser [W]";
input Real Q_bc_input "outgoing heat transfer rate to cold sink [W]";
input Real T_s_init_input; 
input Real Q_draw_nominal_fullpower_input; 
output Real T_s; 

Real Q_cs "instantaneous heat transfer rate from the heat pipe condenser to the Stirling engine [W]";
Real Q_bc "instantaneous heat transfer from the Stirling engine to the cold sink"; 

Real Q_draw_nominal_fullpower "nominal power draw at full power [W]";
Real HTC;

parameter Boolean MODEL_STIRLING_ACTIVATION = false "Require the Stirling engine to reach a setpoint temperature before drawing heat";
parameter Boolean FIX_CONDENSER_HEATFLUX = false;

Boolean STIRLING_ACTIVATED (start = false, fixed=true); //activated with 'when' statement

initial equation

if COLD_START then 
  T_s = 15 + 273.15;
else 
  if T_s_init_input > 1E-9 then 
    T_s = T_s_init_input;
  else     
    der(T_s) = 0.;
  end if; 
end if;   


equation

if Q_draw_nominal_fullpower_input > 1E-9 then 
  Q_draw_nominal_fullpower = Q_draw_nominal_fullpower_input;
else 
  Q_draw_nominal_fullpower = 2250/8;
end if; 

HTC = Q_draw_nominal_fullpower/(T_stirling_nominal - T_stirling_cold_nominal);

if MODEL_STIRLING_ACTIVATION then 
  when T_s > T_stirling_activation then 
    STIRLING_ACTIVATED = true;
    Modelica.Utilities.Streams.print("Stirling engine activated!");
  end when;
else 
  STIRLING_ACTIVATED = true;
end if; 

//Heat flows
if FIX_CONDENSER_HEATFLUX then
  Q_cs = 2250/8;  
else 
  Q_cs = Q_cs_input;
end if;


if Q_bc_input > 1E-9 then //override
  Q_bc = Q_bc_input;
else
  if MODEL_STIRLING_ACTIVATION then 
    if STIRLING_ACTIVATED then 
      Q_bc = HTC*(T_s - T_stirling_cold_nominal);
    else
      Q_bc = 0;
    end if;
  else 
    Q_bc = HTC*(T_s - T_stirling_cold_nominal); 
  end if; 
end if;

//Energy balance
C_stirling*der(T_s) = Q_cs - Q_bc;


end StirlingSimpleOneNode;
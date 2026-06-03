/*
A simple lumped capacitance model of the Stirling engine with a single node. Intended to be a modular version of the one built into the HeatPipeWithVapourCore model. 

Inputs: Heat transfer rate from heat pipe condenser (Q_cs_input)
Outputs: Lumped stirling temperature (T_stirling) 
*/
model StirlingSimpleOneNode

parameter Real C_cond_wall = 17.780965 "Heat pipe condenser wall capacitance [J/K]";
parameter Real C_stirling = 15*C_cond_wall "Stirling capacitance [J/K]"; 

parameter Real T_stirling_activation = 650 + 273.15 "Temperature at which Stirlings are turned on in start-up run of KRUSTY"; 
parameter Real T_stirling_nominal = 650 + 273.15 "Stirling hot-side temperature, Poston et al., Fig. 10 [K]"; //630 + 273.15
parameter Real T_stirling_cold_nominal  = 65 + 273.15 "Cold sink temperature [K]"; 


input Boolean COLD_START; 
input Real Q_cs_input "incoming heat transfer rate from the heat pipe condenser [W]";
input Real Q_bc_input "outgoing heat transfer rate to cold sink [W]";
input Real T_s_init_input; 
input Real Q_draw_nominal_input; 
output Real T_s; 

Real Q_cs "instantaneous heat transfer rate from the heat pipe condenser to the Stirling engine [W]";
Real Q_bc "instantaneous heat transfer from the Stirling engine to the cold sink"; 

Real Q_draw_nominal "nominal power draw [W]";
Real HTC_cold;

parameter Boolean MODEL_STIRLING_ACTIVATION = false "Require the Stirling engine to reach a setpoint temperature before drawing heat";

//input Boolean MODEL_STIRLING_ACTIVATION_input; 
//Boolean MODEL_STIRLING_ACTIVATION;

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

if Q_draw_nominal_input > 1E-9 then 
  Q_draw_nominal = Q_draw_nominal_input;
else 
  Q_draw_nominal = 2350/8;
end if; 

HTC_cold = Q_draw_nominal/(T_stirling_nominal - T_stirling_cold_nominal);

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

if Q_bc_input > 1E-9 then //override
  Q_bc = Q_bc_input;
else
  if MODEL_STIRLING_ACTIVATION then 
    if STIRLING_ACTIVATED then 
      Q_bc = HTC_cold*(T_s - T_stirling_cold_nominal);
    else
      Q_bc = 0;
    end if;
  else 
    Q_bc = HTC_cold*(T_s - T_stirling_cold_nominal); 
  end if; 
end if;

//Energy balance
C_stirling*der(T_s) = Q_cs - Q_bc;


end StirlingSimpleOneNode;
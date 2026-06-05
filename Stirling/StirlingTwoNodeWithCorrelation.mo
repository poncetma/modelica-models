/*
A two-node Stirling engine thermal mass model with the FTB correlation used to get the heat transfer between hot/cold nodes. 

Inputs: Heat transfer rate from heat pipe condenser (Q_cs_input)
Outputs: Hot-side Stirling temperature (T_s, cold-side heat rejection, electric power).  
*/
model StirlingTwoNodeWithCorrelation

parameter Real C_adiabcond_wall = 17.780965 "Heat pipe adiabatic + condenser wall capacitance [J/K]";
parameter Real C_stirling_onenode = 15*C_adiabcond_wall "Stirling capacitance [J/K]"; 
//Both hot and cold sides have thermal masses, but note that the hot-side one is the one that really matters for the system start-up and it has a greater impact on the power throughput.
parameter Real C_stirling_h = 1.5*C_stirling_onenode;
parameter Real C_stirling_c = 1.5*C_stirling_onenode; 


//parameter Real T_stirling_activation = 650 + 273.15 "Temperature at which Stirlings are turned on in start-up run of KRUSTY"; 
parameter Real T_stirling_activation = 950 "Temperature at which Stirlings are turned on in start-up run of KRUSTY"; 
parameter Real T_stirling_nominal = 650 + 273.15 "Stirling hot-side temperature, Poston et al., Fig. 10 [K]"; //630 + 273.15
parameter Real T_stirling_cold_nominal  = 65 + 273.15 "Cold-side temperature [K]"; 
parameter Real T_sink = 15 + 273.15 "Chiller temperature (T_inf) [K]"; 

input Boolean COLD_START; 
input Real Q_cs_input "incoming heat transfer rate from the heat pipe condenser [W]";
input Real Q_bc_input "desired boundary condition, here applied to internal heat throughput (Q_internal) [W]";
input Real T_s_init_input; 
input Real Q_draw_nominal_input; 

output Real T_s(start = T_stirling_nominal, fixed=false) "Stirling hot side (acceptor) temperature [K]"; 
Real T_sc (start = T_stirling_cold_nominal, fixed=false) "Stirling cold side temperature [K]";


Real Q_cs "instantaneous heat transfer rate from the heat pipe condenser to the Stirling engine [W]";
Real Q_internal "instantaneous internal heat transfer from hot to cold sides ('throughput') [W]";
Real Q_bc "instantaneous heat transfer from the Stirling engine to the cold sink"; 

Real Q_draw_nominal "nominal power draw [W]";
output Real Q_rejected_th "thermal power rejected [W]";
output Real Q_electric "electric power generated [W]";

parameter Boolean MODEL_STIRLING_ACTIVATION = true "Require the Stirling engine to reach a setpoint temperature before drawing heat";

Boolean STIRLING_ACTIVATED (start = false, fixed=true); //activated with 'when' statement

//FTB engine performance correlation parameters
parameter Real T_h_min = 400 + 273.15;
parameter Real T_h_max = 650 + 273.15;
parameter Real T_c_max = 70 + 273.15;
parameter Real T_c_min = 30 + 273.15;
parameter Real Q_Th_coeff = (80-52)/200.0; //approximate slope from plot
parameter Real Q_Tc_coeff = (68-86)/40.0; //opposite dependence
parameter Real eta_Th_coeff = (0.35-0.25)/250;
parameter Real eta_Tc_coeff = (0.30-0.35)/40;

//parameter Real eta = 0.325 "thermal efficiency @T_h = 650 C, T_c = 60 c"; 

Real eta "thermal efficiency";

Real eta_nominal;
Real Q_corr_nominal;
Real cf;

Real HTC_cold;


initial equation

if COLD_START then 
  T_s = 15 + 273.15;
  T_sc = 15 + 273.15;
else 
  if T_s_init_input > 1E-9 then 
    T_s = T_s_init_input;
    T_sc = T_s_init_input;
  else     
    der(T_s) = 0.;
    der(T_sc) = 0.;
  end if; 
end if;   


equation

if Q_draw_nominal_input > 1E-9 then 
  Q_draw_nominal = Q_draw_nominal_input;
else 
  Q_draw_nominal = 2350/8;
end if; 

//Compute correction factor to give the correct power draw at nominal T_h, T_c
eta_nominal = (0.25 + eta_Th_coeff*(T_stirling_nominal-T_h_min))*(0.35 + (eta_Tc_coeff*(T_stirling_cold_nominal - T_c_min)))/0.35;
Q_corr_nominal = ( (52 + Q_Th_coeff*(T_stirling_nominal-T_h_min))*( 86 + (Q_Tc_coeff*(T_stirling_cold_nominal-T_c_min)) )/86 )/eta_nominal;
cf = Q_draw_nominal/Q_corr_nominal; 

//Compute HTC for cold side -> Give the correct T_sc at steady-state with presumed cooler/sink temperature
//HTC_cold = (Q_draw_nominal*(1-eta_nominal))/(T_stirling_cold_nominal - T_sink);
HTC_cold = (Q_draw_nominal)/(T_stirling_cold_nominal - T_sink);


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


eta =  max( 0.15,(0.25 + eta_Th_coeff*(T_s-T_h_min)) )*( max( 0.15, 0.35 + (eta_Tc_coeff*(T_stirling_cold_nominal - T_c_min)) )/0.35);
//eta = 0.325; 

if Q_bc_input > 1E-9 then //override
  Q_internal = Q_bc_input;
else
  if MODEL_STIRLING_ACTIVATION then 
    if STIRLING_ACTIVATED then 
      //Q_bc = HTC_cold*(T_s - T_stirling_cold_nominal);
      //Now use the correlation 
      //Q_bc = ( (52 + Q_Th_coeff*(T_s-T_h_min))*( 1 + (Q_Tc_coeff*(T_stirling_cold_nominal-T_c_min))/86 ) )/eta;       
      Q_internal = max(0, cf*( ( (52 + Q_Th_coeff*(T_s-T_h_min))*( 86 + (Q_Tc_coeff*(T_stirling_cold_nominal-T_c_min)) )/86 )/eta ) );       
    else
      Q_internal = 0;
    end if;
  else 
    //Q_bc = HTC_cold*(T_s - T_stirling_cold_nominal); 
    Q_internal = max(0, cf* (( (52 + Q_Th_coeff*(T_s-T_h_min))*( 86 + (Q_Tc_coeff*(T_stirling_cold_nominal-T_c_min)) )/86 )/eta ) );       
    
  end if; 
end if;

Q_bc = HTC_cold*(T_sc - T_sink); 

//Energy balance
//C_stirling*der(T_s) = Q_cs - Q_bc;
C_stirling_h*der(T_s) = Q_cs - Q_internal;
C_stirling_c*der(T_sc) = Q_internal - Q_bc;

//other outputs
Q_electric = eta*Q_internal;
Q_rejected_th = (1-eta)*Q_internal;

end StirlingTwoNodeWithCorrelation;
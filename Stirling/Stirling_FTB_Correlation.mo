model Stirling_FTB_Correlation
"Simple correlation for power output of the 'FTB' (ASC-like) Stirling engine as a function of 
hot- and cold-side temperatures. "

//Approximately linear dependence on both hot- and cold-side temperatures (Wood et al.), with one side held constant
//Ultimately, it depends on the total delta T. 

parameter Real T_h_min = 400 + 273.15;
parameter Real T_h_max = 650 + 273.15;
parameter Real T_c_max = 70 + 273.15;
parameter Real T_c_min = 30 + 273.15;

parameter Real T_h_coeff = (80-52)/200.0; //approximate slope from plot
parameter Real T_c_coeff = (68-86)/40.0; //opposite dependence
Real Q;
Real T_h;
Real T_c;
equation
//base dependence on T_h with a modifier for T_c
T_h = 500 + 273.15;
T_c = 50 + 273.15;
Q = (52 + T_h_coeff*(T_h-T_h_min))*( 1 + (T_c_coeff*(T_c-T_c_min))/86 ); 
end Stirling_FTB_Correlation;

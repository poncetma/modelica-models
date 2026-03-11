model PointKineticsWithDecayHeat
"
Point kinetics with decay heat tracking, based on effective precursor groupings (similar to DNPs).
This is the same general approach used in FRINK, though the exact implementation there is not published. 
This model uses the TRACE implementation (a power formulation, see TRACE V5 theory manual) 
with input data suited for KRUSTY (fast-spectrum, U238 fission products) from the ANS Standard for Decay Heat (1979).

Inputs: Fuel temperature (average, T_Fuel), external reactivity (rho_ext)
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

parameter Real E_f = 195.0 "Energy release per fission, U238 [MeV]"; 
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
} "Decay heat group decay constants for U238 [s^-1]";
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
} "Decay power release per fission for U238 [MeV/s]";
parameter Real E_fracs[23] = EDs./lambdas_dh./E_f "Decay heat fraction [-]";

parameter Real rho_0 = 0.0 "initial reactivity"; 

//hardcoded experimental data from KRUSTY (for power history)
constant Real exp_times[138] = {
0,
1116,
1435,
1754,
3429,
3748,
4147,
5503,
5981,
6619,
7097,
21213,
21372,
21851,
21930,
22409,
22887,
23206,
23525,
24004,
24642,
24961,
25280,
25599,
26157,
26635,
27194,
28789,
29427,
30065,
30623,
31261,
32776,
33334,
33813,
34291,
34929,
35328,
35806,
36125,
36524,
37003,
37481,
37800,
38358,
38837,
39634,
40192,
40671,
41229,
41867,
42346,
42665,
43223,
43781,
44419,
44897,
45695,
46014,
46732,
47370,
47928,
48566,
49204,
49603,
50320,
50878,
51357,
51835,
52314,
53032,
53989,
54547,
55025,
55663,
56062,
56780,
57577,
57896,
58454,
59013,
59491,
60049,
61086,
61485,
62043,
62601,
63159,
63718,
64356,
64675,
65153,
65791,
66270,
66828,
67227,
67705,
68263,
68662,
69141,
69699,
70177,
70656,
71852,
72649,
73367,
74005,
74803,
75600,
76238,
76876,
77514,
78072,
78551,
79109,
81422,
82139,
83016,
83814,
86206,
86924,
87482,
88120,
89795,
90433,
90911,
91390,
91709,
92427,
96813,
96972,
97132,
97929,
98727,
99444,
100162,
100641,
100880
} "linearised experimental data, times [s]";

constant Real exp_powers[138] = {
0,
0,
3675,
3072,
3072,
3603,
3072,
3072,
2921,
2986,
2921,
2950,
0,
0,
2921,
2756,
2971,
2914,
2555,
2993,
2813,
2203,
2935,
2218,
2885,
2612,
2684,
2699,
1866,
2146,
1995,
2053,
2038,
2871,
2591,
2763,
2677,
2749,
2706,
2742,
4163,
3861,
4048,
4012,
4077,
4033,
4041,
2541,
2871,
2684,
2742,
2663,
2713,
2656,
2519,
2598,
2569,
2892,
2828,
2907,
2605,
2713,
2699,
2978,
2957,
2993,
3940,
3660,
3797,
3718,
3754,
3754,
2447,
2770,
2548,
2620,
2548,
2555,
1313,
3014,
2182,
2498,
2390,
2397,
1313,
2813,
2067,
2361,
2239,
2318,
2275,
5160,
2160,
3014,
2627,
2785,
2699,
2734,
1522,
3093,
2368,
2648,
2541,
2577,
1234,
1665,
1450,
1522,
1464,
1493,
1457,
1486,
1443,
1493,
1457,
1464,
1313,
1364,
1328,
1335,
1952,
1737,
1830,
1780,
2699,
2490,
2634,
2605,
2670,
2742,
0,
2742,
1148,
1550,
1335,
1421,
1378,
0
} "linearised experimental data, powers [W]";

//input Real t_exp_onset "time at which the transient began during the 28-h run of KRUSTY [s]";
parameter Real t_exp_onset = 8*3600; //10.02*3600 "time at which the transient began during the 28-h run of KRUSTY [s]"; 
parameter Integer nearest_exp_index = findNearestIndex(exp_times, t_exp_onset);

Real[nearest_exp_index] fissionpower_history = Modelica.Math.Vectors.reverse(exp_powers[1:nearest_exp_index]);
Real[nearest_exp_index] time_history = { exp_times[nearest_exp_index] - t for t in Modelica.Math.Vectors.reverse(exp_times[1:nearest_exp_index])} ;

/*
  This correlation gives the integral reactivity feedback (relative to *operating* (nominal) temp).
  From "KRUSTY Design and Modelling" slide 88
*/
function alpha_Tf_poly 
  input Real T_fuel;
  output Real alpha_Tf;
  algorithm
    alpha_Tf := (-1.6951E-11*T_fuel^3. + 5.0121E-8*T_fuel^2 - 1.4888E-4*T_fuel -7.9756E-2)*0.01*Beta; //in pcm/K
end alpha_Tf_poly;

/*
  Finds the nearest value in an array and returns the corresponding index
*/
function findNearestIndex
  input Real array[:];
  input Real target;
  output Integer index;
  protected 
    Real minDiff;
    Real currentDiff;
  algorithm
    index := 1;
    minDiff := abs(array[1] - target);
    
    for i in 2:size(array, 1) loop
      currentDiff := abs(array[i] - target);
      if currentDiff < minDiff then
        minDiff := currentDiff;
        index := i;
      end if;
    end for;
end findNearestIndex;


/*
  Computes initial conditions based on fission power history (excluding decay heat). 
  Adapted from the Modelica TRANSFORM lib  which applies the TRACE methodology. 
  For KRUSTY analysis, the exp. fission power curve gives this type of power history.
  Since the experimental data is constant, it is hardcoded (allowing it to be packaged in the FMU).
*/
function computeICsFromPowerHistory
  input Real[:] time_hist;
  input Real[:] fisspower_hist; //defined as only the fission power (no decay heat in the history)  
  input Real[:] lambdas;
  input Real[:] betas;
  input Real Lambda;
  input Real[:] lambdas_dh;
  input Real[:] E_fracs "Decay-heat fraction of fission power";
  
  //input Real[size(lambdas, 1)] Cs_0=fill(0, size(lambdas, 1)) "Precursor concentration at history time = 0";
  //input Real[size(lambdas_dh, 1)] Es_0=fill(0, size(lambdas_dh, 1)) "Decay-heat concentration at history time = 0";
  output Real Cs_0[6]; //Real[size(lambdas, 1)] Cs_0;
  output Real Hs_0[23];//Real[size(lambdas_dh, 1)] Hs_0;
  protected
    //Integer nT = nearest_index; //nT=138 "# of time history points to include (integrate up to)";
    Integer nT = size(time_hist, 1); //nT=138 "# size of time history";
    Integer nK=6 ;//size(lambdas, 1) "# of delayed-neutron precursors groups";
    Integer nDH=23; //size(lambdas_dh, 1)"# of decay-heat groups";
    Real dt;
    Real a; //y-intercept
    Real b; //slope
    Real elamdt;
    Real elamdt_dh;  
  
  algorithm //manually implemented from TRACE manual
  
    for k in 1:nK  loop       
      Cs_0[k] := 0.0; //assumes that the full history is provided      
      for i in 1:nT - 1 loop 

        //slope and y-intercept computed based on time history indexing 
        b := (fisspower_hist[i] - fisspower_hist[i+1])/(time_hist[i] - time_hist[i+1]); //slope
        a := fisspower_hist[i] + time_hist[i]*b; //y-intercept        
        dt := time_hist[i+1] - time_hist[i];
        Cs_0[k] := Cs_0[k] 
                + ( (a - b/lambdas[k])*(1 - exp(-lambdas[k]*dt)) + b*(time_hist[i+1]*exp(-lambdas[k]*dt) - time_hist[i]) )*exp(-lambdas[k]*time_hist[i]);               
      end for;
      
      Cs_0[k] := betas[k]/lambdas[k]/Lambda*Cs_0[k];     

    end for;
    
    for j in 1:nDH  loop       
      Hs_0[j] := 0.0; //assumes that the full history is provided
      for l in 1:nT - 1 loop 
        //slope and y-intercept computed based on time history indexing 
        b := (fisspower_hist[l] - fisspower_hist[l+1])/(time_hist[l] - time_hist[l+1]); //slope
        a := fisspower_hist[l] + time_hist[l]*b; //y-intercept        
        dt := time_hist[l+1] - time_hist[l];        
        
        Hs_0[j] := Hs_0[j] 
                + ( (a - b/lambdas_dh[j])*(1 - exp(-lambdas_dh[j]*dt)) + b*(time_hist[l+1]*exp(-lambdas_dh[j]*dt)  - time_hist[l] ) ) 
                * exp(-lambdas_dh[j]*time_hist[l]);               
      end for;
      
      Hs_0[j] := E_fracs[j]/lambdas_dh[j]*Hs_0[j];     

    end for;        
  
end computeICsFromPowerHistory;

Real rho_ext; 
Real T_fuel; 
Real rho "instantaneous net reactivity";
Real rho_fb "reactivity feedback";
Real alpha_Tf;
Real T_fuel_ref;

Real Cs[6] "instaneous group-wise precursor power [W] "; 
Real Hs[23] "decay heat precursor energy [J]";
//Real Cs_0[6];
//Real Hs_0[23];
//Real decay_heat_fraction;

output Real P_fiss "instantaneous fission power"; 
output Real P_dec "instantaneous decay heat power";
output Real P_tot "total instantaneous thermal power";

/*
initial equation

if P_0 > 0 then  
  P_tot = P_0;
else  
  P_tot = 5000.0;
end if; 


if (t_exp_onset < 0) then //any negative input will trigger the infinite power history IC
  betas/Lambda*P_fiss = lambdas.*Cs;
  Es*P_fiss = lambdas_dh.*Hs ;
else
  
  (Cs, Hs) = computeICsFromPowerHistory(powerhist_times,powerhist_powers,lambdas, betas, Lambda, lambdas_dh,Es);   
  
end if;
*/

/*
initial algorithm
  (Cs, Hs) := computeICsFromPowerHistory(powerhist_times,powerhist_powers,lambdas, betas, Lambda, lambdas_dh,Es);   
  P_fiss := 0;
*/

  
initial equation
//  (time_history, fissionpower_history) = genFissionPowerHistory(exp_times, exp_powers, nearest_exp_index);

  (Cs, Hs) = computeICsFromPowerHistory(time_history, fissionpower_history, lambdas, betas, Lambda, lambdas_dh, E_fracs);   
  P_fiss = exp_powers[nearest_exp_index]; //The fission power has to match the chosen point in the history
  

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
else //rho_ext = 0.0;
  if time < 100 then
    rho_ext = 0.0;
    //rho_ext= min(100, time)*1e-5 ;
  else
    rho_ext = 0.0;
  end if;
  
end if;

rho_fb = alpha_Tf*(T_fuel - T_fuel_ref); 
rho = rho_0 + rho_ext + rho_fb; 

der(P_fiss) = (rho - Beta)/Lambda.*P_fiss + sum(lambdas.*Cs);
der(Cs) = betas/Lambda*P_fiss - lambdas.*Cs;

der(Hs) = E_fracs*P_fiss - lambdas_dh.*Hs;
P_dec = sum(lambdas_dh.*Hs);
P_tot = P_fiss + P_dec;
//decay_heat_fraction = P_dec/P_tot;


end PointKineticsWithDecayHeat;
model FPSE_Nonlinear_Majidniya
 "Implementation of Majidniya et al. (2019) Eqs. (15)/(16) for the RE-1000 case.
   Uses Table 2 parameter values (converted to SI).
   Implements EXACT algebraic structure of Eq.(15) and Eq.(16) with named intermediates."
  import Modelica.Constants.pi;

  // ------------------- Paper constants (Table 2, Majidniya) -------------------
  // Temperatures (Table 2)
  //parameter Real T_h = 814.3 "K (heater temperature from Table 2)";
  parameter Real T_h_nominal = 814.3; 
  parameter Real T_k = 322.8 "K (cooler/regenerator reference temp from Table 2)";
  //parameter Real T_r = (T_h_nominal-T_k)/log(T_h_nominal/T_k) "K (regenerator log-mean used in paper)";
  parameter Real T_b = T_k "K (buffer assumed same as cold side in paper)";

  // Mean pressure (Table 2)
  // Note: Table 2 lists "Pmean 71 bars" -> 71 bar = 7.1e6 Pa
  parameter Real p_mean = 70e5; //71e5 "Pa (71 bar) - RE-1000 case from Table 2";

  // Volumes (Table 2)
  parameter Real V_D = 37.97e-6 "m^3 (VD from table: 37.97 cm^3)";
  parameter Real V_r = 56.37e-6 "m^3 (Vr from table: 56.37 cm^3)";
  //parameter Real V_k = A_r*L_r; //computed from total XS area //V_r //assume 
  parameter Real V_k = A_k*L_k; //computed from total XS area //V_r //assume 
  parameter Real V_h = A_h*L_h; //computed from total XS area //V_r; //assume
  parameter Real V_B = 2615.0e-6 "m^3 (VB from table: 2615 cm^3 -> buffer large)";
  // mean compression/expansion clearance lengths (Cc,Ce) from Table 2:  
  parameter Real C_e = 1.861e-2 "m expansion clearance Ce (Table2)";
  parameter Real C_c = 1.83e-2  "m compression clearance Cc (one side) - paper uses two?";
  // The paper uses geometric amplitudes: Xp and Xd (Table 2)  
  //parameter Real X_p = 1.145e-2 "m (Xp = 1.145 cm - from Table 2)";
  //parameter Real X_d = 1.233e-2 "m (Xd = 1.233 cm - from Table 2)";
  parameter Real X_p = 9.0e-3 "m, from PhD thesis";
  parameter Real X_d = 8.3e-3 "m, from PhD thesis";

  // Cross-sectional areas from Table 2 (convert cm^2 to m^2)
  parameter Real A_h = 1.4898e-4 "m^2 (Ah from Table2)";
  parameter Real A_k = 2.6163e-4 "m^2 (Ak from Table2)";
  //parameter Real A_k = 135*
  parameter Real A_r = 8.745e-4  "m^2 (Ar from Table2)";

  // Diameters  from Table 2 (convert cm to m)
  parameter Real d_p = 5.718e-2 "m piston diameter (dp)";
  parameter Real d_d = 5.67e-2  "m displacer diameter (dd)";
  parameter Real d_rod = 1.663e-2 "m rod diameter (drod)";

  // Cross-sectional areas used in geometry (derived)
  parameter Real A_p = pi*(d_p/2)^2 "m^2 power piston area (derived from dp)";
  parameter Real A_d = pi*(d_d/2)^2 "m^2 displacer area (derived from dd)";
  parameter Real A_rod = pi*(d_rod/2)^2 "m^2 rod area (derived from drod)";

  // Masses & mechanical params (Table 2)
  parameter Real m_d = 0.426 "kg (displacer mass from Table 2)";
  parameter Real m_p = 6.2   "kg (power piston mass from Table 2)";
  parameter Real K_p = 0 "N/m (piston centering spring, not present in Majidniya model))";
  parameter Real K_d = 0.0 ;//5000 "N/m (displacer spring), not in the Majidniya model but allows final adjustment";
  // Viscous damping (paper uses load damping = 580 N s/m)
  parameter Real C_load  = 580.0;//580.0 "N*s/m (load damping from Table 2)";
  parameter Real C_p = 0.0 "additional (viscous) piston damping, assume negligible"; 
  parameter Real C_d = 0.0 "displacer mechanical viscous damping (small / not explicit in table)";

  // Gas constants
  parameter Real R_gas = 2077.0 "J/(kg K) helium specific gas constant (paper uses He)";
  parameter Real gamma = 5195.0/3117.0; // "polytropic exponent ~1.66666
  parameter Real T_mean = (T_h_nominal + T_k)/2; //Don't have enough data to improve the T_avg estimate
  parameter Real rho_mean = p_mean / ( R_gas * T_mean ); 
  parameter Real mu = 30e-6 "dynamic viscosity, approximate (src: 'Stirling Cycle Machine Analysis') "; //19.85e-6

  // Pressure drop geometry & friction (Table 2)
  parameter Real L_h = 0.1834 "m (Lh 18.34 cm)";
  parameter Real L_k = 0.0792 "m (Lk 7.92 cm)";
  parameter Real L_r = 0.0644 "m (Lr 6.44 cm)";

  parameter Real d_h_h = 0.2362e-2; //4.*V_h/A_w_h; //0.2362e-2 "m (dh 0.2362 cm -> 0.002362 m)";
  parameter Real d_h_k = 4.*V_k/A_w_k; //d_h_h;
  parameter Real d_w = 0.00889e-2 "Wire mesh diameter, m (dw from table 0.00889 cm -> 8.89e-5 m)";
  parameter Real porosity = 0.759; //Also given in table 2
  parameter Real d_h_r = (d_w * porosity)/(1 - porosity);
  // Wetted flow areas (approx - Table 2 gives Awk 115.2 cm2 -> 0.01152 m2)
  parameter Real A_w_k = 115.2e-4 "m^2 (Awk from table 115.2 cm^2)";
  //parameter Real A_w_h = A_w_k "use same order as Aw_k if not separately given"; //Actually, don't need it since we have d_h_h (ostensibly)
  //parameter Real A_w_r = 4*V_r/d_h_r; //A_r "approx use regenerator area Ar from table (converted above)"; //At no point do we actually need the wetted area of the regenerator (we get d_h_r for free)

  
  
  // Numerical safe values
  //parameter Real eps = 1e-9;
  
  // Expected solution
  parameter Real omega = 2*pi*30.0; // 30 Hz
  parameter Real phi = -42.5*pi/180; //phase shift
  //parameter Real Vdot_max = 0.0080674; //graphically solved max from analytical expression
                                         //analytically solved max in closed form 
  parameter Real Vdot_max = abs(A_p*X_p*omega*sin(phi)/sin(atan((A_p*X_p*omega*sin(phi))/(A_p*X_p*omega*cos(phi)-(2*A_d-A_rod)*X_d*omega))));

  // ------------------- STATES -------------------
  Real x_p(start=X_p*sin(phi)) "power piston displacement (m)"; //Must start with correct phase, X_p*sin(omega*t + phi), t=0  //
  //Real x_p(start=X_p*sin(-pi/2)) "power piston displacement (m)"; //Try with "perfect" phase
  //Real v_p(start=X_p*omega*cos(phi)) "power piston velocity (m/s)"; //Need a little kick? //Can use the same trick, analytical derivative
  Real v_p(start=0.0); //just keep it at 0, better matches the pressure drop from Majidniya et al.

  Real x_d(start=0.0) "displacer displacement (m)";
  //Real v_d(start=X_d*omega) "displacer velocity (m/s)";
  Real v_d(start=0.0);

  // ------------------- Intermediates (named to match paper) -------------------
  Real Vc "compression instantaneous volume";
  Real Ve "expansion instantaneous volume";
  Real Vd "displacer gas-spring instantaneous volume";
  Real Vb "buffer instantaneous volume";

  //Real V_avT "sum term defined in Eq.(7) (units m^3/K)";
  //Real MR "intermediate MR = m_total * R (paper notation)";
  Real p;
  Real p_b;
  Real p_d;
  
  Real deltaP(start=0.0) "pressure drop driving displacer (from heat exchanger/regenerator pressure drops)";

  // For pressure-drop velocity proxies
  //Real Vdot_c, Vdot_e;
  Real Vdot_flow;
  Real u_h, u_k, u_r;  
  Real u_h_max, u_k_max, u_r_max;
   
  Real Re_h; 
  Real Re_k; 
  Real Re_r; 
  Real Cf_h "heater Darcy friction factor (tunable)";
  Real Cf_k "cooler Darcy friction factor (tunable)";
  Real Cf_r  "regenerator Darcy friction factor (tunable)";
  
  Real dP_h;   
  Real dP_k;   
  Real dP_r; 

  Real F_load; 
  Real P_out_mech, E_out_mech, P_mech_avg; 

  //Dynamically calculated densities and viscosities following Majidniya's thesis
  Real rho_h;
  Real rho_k;
  Real rho_r; 
  Real mu_h;
  Real mu_k;
  Real mu_r;
  
  
  Real T_h; //Try with time-varying signal
  Real T_r; //If T_h varies then so does T_r.
  Real V_avT;     
  
equation
  T_h = 814.3 + 0*time; // + 100*sin(omega*time);
  T_r = (T_h_nominal-T_k)/log(T_h_nominal/T_k) "K (regenerator log-mean used in paper)";
  V_avT =  A_p*C_c/T_k + A_d*C_e/T_h + V_h/T_h + V_k/T_k + V_r/T_r;

  // ---------- instantaneous geometric volumes (Eqs. 3,4,13,14 in paper) ----------
  // Clearances: Ce, Cc used to set mean volumes (paper notation uses Ce, Cc)
  // For RE-1000, we treat Ve mean and Vc mean implicitly via the Table2 amplitude & clearances:
  Vc = A_p*(x_p + C_c) - (A_d - A_rod)*x_d;
  Ve = A_d*(x_d + C_e);
  
  Vb = V_B - A_p*x_p;   // buffer instantaneous volume (Eq.13 style)
  Vd = V_D - A_rod*x_d; // displacer gas spring instantaneous volume (Eq.14 style)

  // ---------- VavT and M definitions (Eq.6 & 7 mapping) ----------
  //MR = p_mean*(Vc/T_k + Ve/T_h + V_h/T_h + V_k/T_k + V_r/T_r);
  //MR = p_mean*((A_p*C_c)/T_k + (A_d*C_e)/T_h + V_h/T_h + V_k/T_k + V_r/T_r); 
  
  //p is taken tto be the pressure in the compression volume
  p = p_mean*(1 + (A_p*x_p - (A_d - A_rod)*x_d)/(T_k*V_avT) + (A_d*x_d)/(T_h*V_avT) )^(-1) ;
  
  p_b = p_mean*(V_B/Vb)^gamma; //Gas springs, may need slight adjustment
  p_d = p_mean*(V_D/Vd)^gamma;
  
  // ---------- Pressure drop DeltaP (Eq.17 & Eqs.18-19) ----------
  // volumetric rates from Eqs.(18) (paper): Vdot_c = A_p*v_p - A_d*v_d ; Vdot_e = A_d*v_d  
  Vdot_flow = A_p*v_p - (2*A_d - A_rod)*v_d; //Vdot_c - Vdot_e;   
  //Vdot_flow = der(Vc) - der(Ve);

  // instantaneous velocities (numerical)    
  u_h = Vdot_flow/A_h;
  u_k = Vdot_flow/A_k;
  u_r = Vdot_flow/A_r;
  
  // max velocities (analytical)
  //Areas don't seem to be consistent with hydraulic diameter estimates? That's because there isn't just a single pipe for the heater/cooler!
  //With these proper areas, the velocities come down to much more realistic levels
  /*
  u_h_max = omega*sqrt( (A_p*X_p)^2 + ((2*A_d-A_rod)*X_d)^2 - 2*A_p*X_p*(2*A_d - A_rod)*X_d*sin(phi) )/A_h;
  u_k_max = omega*sqrt( (A_p*X_p)^2 + ((2*A_d-A_rod)*X_d)^2 - 2*A_p*X_p*(2*A_d - A_rod)*X_d*sin(phi) )/A_k;
  u_r_max = omega*sqrt( (A_p*X_p)^2 + ((2*A_d-A_rod)*X_d)^2 - 2*A_p*X_p*(2*A_d - A_rod)*X_d*sin(phi) )/A_r;
  */
  //solve using predicted Vdot_max calculated previously
  
  u_h_max = Vdot_max/A_h;
  u_k_max = Vdot_max/A_k;
  u_r_max = Vdot_max/A_r;
  
  
  //friction factors (Cf)     
  rho_k = (48.14*(p + dP_k/2.0)*1e-5)/(T_k*(1 + 0.4446*(p + dP_k/2.0)*1e-5*T_k^(-1.2)));
  rho_r = (48.14*(p + deltaP - dP_r/2.0)*1e-5)/(T_r*(1 + 0.4446*(p + deltaP - dP_r/2.0)*1e-5*T_r^(-1.2)));
  rho_h = (48.14*(p + dP_k + dP_r/2.0)*1e-5)/(T_h*(1 + 0.4446*(p + dP_k + dP_r/2.0)*1e-5*T_h^(-1.2)));
  mu_k = 3.674*1e-7*T_k^0.7;
  mu_r = 3.674*1e-7*T_r^0.7;
  mu_h = 3.674*1e-7*T_h^0.7;
  
  //compute from predicted max flow velocity based on sinusoidal motion   
  /*
  Re_h = rho_mean*u_h_max*d_h_h/mu;
  Re_k = rho_mean*u_k_max*d_h_k/mu;
  Re_r = rho_mean*u_r_max*d_h_r/mu;    
  */
  Re_h = rho_h*abs(u_h_max)*d_h_h/mu_h;
  Re_k = rho_k*abs(u_k_max)*d_h_k/mu_k;
  Re_r = rho_r*abs(u_r_max)*d_h_r/mu_r;
  
  //compute dynamically according to instantaneous local velocity--this approach seems cause damping to zero, for some reason.
  /*
  Re_h = max(rho_mean*u_h*d_h_h/mu, eps);//again guard against division by 0
  Re_k = max(rho_mean*u_k*d_h_k/mu, eps);
  Re_r = max(rho_mean*u_r*d_h_r/mu, eps);
  */
  if Re_h < 2000 then
    Cf_h = 64/Re_h;
  else    
    Cf_h = 0.316*Re_h^(-0.25); //Assume turbulent conditions for now
  end if;
  
  if Re_k < 2000 then
    Cf_k = 64/Re_k;
  else     
    Cf_k = 0.316*Re_k^(-0.25); 
  end if;
  
  if Re_r < 60 then 
    Cf_r = 4*10^(1.73 - 0.93 * log10(Re_r));    
  elseif Re_r < 1000 then
    Cf_r = 4*10^(0.714 - 0.365 * log10(Re_r));
  else
    Cf_r = 4*10^(0.015 - 0.125 * log10(Re_r));
  end if;
  

  // nonlinear Darcy pressure drops in the "bypass" region (heater, regenerator, cooler) (Eq.17)
  // DeltaP_i = Cf_i * (L_i/dh_i) * 0.5 * rho * u_i * |u_i|
  /*
  dP_h = 0.5 * rho_mean * Cf_h * (L_h/d_h_h) * u_h * abs(u_h);
  dP_k = 0.5 * rho_mean * Cf_k * (L_k/d_h_k) * u_k * abs(u_k);
  dP_r = 0.5 * rho_mean * Cf_r * (L_r/d_h_r) * u_r * abs(u_r);
  */
  dP_h = 0.5 * rho_h * Cf_h * (L_h/d_h_h) * u_h * abs(u_h);
  dP_k = 0.5 * rho_k * Cf_k * (L_k/d_h_k) * u_k * abs(u_k);
  dP_r = 0.5 * rho_r * Cf_r * (L_r/d_h_r) * u_r * abs(u_r);
  // The paper sums these into a DeltaP that appears in Eq (16) as DeltaP: total drop between hot and cold
  deltaP = dP_h + dP_r + dP_k;
    
  F_load = C_load*v_p; //damping force due to power conversion 
  P_out_mech = F_load*v_p; 
  der(E_out_mech) = P_out_mech;  
  // average power computed externally (user can post-process or use E_elec_int/time)
  //P_mech_avg = E_out_mech/ max(time, 1e-5); 
  P_mech_avg = ( E_out_mech - delay(E_out_mech, 0.1) ) / 0.1;//moving average with 0.1s window

  //Mechanics: introduce fudge factors (spring consts, damping coeffs) to tune frequency, amplitude.
  der(x_p) = v_p;  
  der(v_p) = ( A_p*(p - p_b)  - F_load - C_p*v_p - K_p*x_p )/m_p; 

  der(x_d) = v_d;  
  der(v_d) = ( A_d*deltaP + A_rod*(p - p_d) - C_d*v_d - K_d*x_d )/ m_d; 
  //der(v_d) = ( A_d*deltaP + (A_d - A_rod)*p - A_rod*p_d  )/ m_d; //wrong!

  // Simple assertions to help catch blow-ups
  //assert(V_B - A_p * x_p > -0.9*V_B, "Compression volume (V_B - A_p x_p) approaching zero or negative");
  //assert(V_D - A_rod * x_d > -0.9*V_D, "Displacer-spring volume (V_D - A_rod x_d) approaching zero or negative");

  //annotation (experiment(StartTime=0, StopTime=0.15, Tolerance=1e-5, Interval=1e-5));
end FPSE_Nonlinear_Majidniya;
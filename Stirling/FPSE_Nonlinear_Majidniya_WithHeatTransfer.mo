model FPSE_Nonlinear_Majidniya_WithHeatTransfer
 "Implementation of Majidniya et al. (2019) Eqs. (15)/(16) for the RE-1000 case.
   Uses Table 2 parameter values (converted to SI).
   "
  import Modelica.Constants.pi;

  // ------------------- Paper constants (Table 2, Majidniya) -------------------
  // Temperatures (Table 2)
  //parameter Real T_h = 814.3 "K (heater temperature from Table 2)";
  parameter Real T_h_nominal = 814.3; 
  parameter Real T_k_nominal = 322.8 "K (cooler/regenerator reference temp from Table 2)";
  //parameter Real T_r = (T_h_nominal-T_k)/log(T_h_nominal/T_k) "K (regenerator log-mean used in paper)";
  parameter Real T_b = T_k_nominal "K (buffer assumed same as cold side in paper)";

  // Mean pressure (Table 2)
  // Note: Table 2 lists "Pmean 71 bars" -> 71 bar = 7.1e6 Pa
  parameter Real p_mean = 71e5; //71e5 "Pa (71 bar) - RE-1000 case from Table 2, equal to MR/V_avT";

  // Volumes (Table 2)
  parameter Real V_D = 37.97e-6 "m^3 (VD from table: 37.97 cm^3)";
  parameter Real V_r = 56.37e-6 "m^3 (Vr from table: 56.37 cm^3)";
  parameter Real V_r1 = V_r/2; 
  parameter Real V_r2 = V_r/2;
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
  parameter Real B_load = 10 "damping coefficient"; 
  parameter Real C_p = 0.0 "additional (viscous) piston damping, assume negligible"; 
  parameter Real C_d = 0.0 "displacer mechanical viscous damping (small / not explicit in table)";

  // Gas constants
  parameter Real R_gas = 2077.0 "J/(kg K) helium specific gas constant (paper uses He)";
  parameter Real gamma = 5195.0/3117.0; // "polytropic exponent ~1.66666
  parameter Real T_mean_approx = (T_h_nominal + T_k_nominal)/2; //Don't have enough data to improve the T_avg estimate
  //parameter Real rho_mean = p_mean / ( R_gas * T_mean_approx ); 
  parameter Real m_total_gas = 0.0014 "total gas content in engine [kg]"; //From thesis 
  parameter Real mu = 30e-6 "dynamic viscosity, approximate (src: 'Stirling Cycle Machine Analysis') "; //19.85e-6

  // Pressure drop geometry & friction (Table 2)
  parameter Real L_h = 0.1834 "m (Lh 18.34 cm)";
  parameter Real L_k = 0.0792 "m (Lk 7.92 cm)";
  parameter Real L_r = 0.0644 "m (Lr 6.44 cm)";

  parameter Real d_h_h = 0.2362e-2; //4.*V_h/A_w_h; //0.2362e-2 "m (dh 0.2362 cm -> 0.002362 m)";
  parameter Real d_h_k = d_h_h; //4.*V_k/A_w_k; //d_h_h;
  //parameter Real d_h_k = A_w_k/pi/L_k; //equivalent formulation of the above
  parameter Real d_w = 0.00889e-2 "Wire mesh diameter, m (dw from table 0.00889 cm -> 8.89e-5 m)";
  parameter Real porosity = 0.759; //Also given in table 2
  parameter Real d_h_r = (d_w * porosity)/(1 - porosity);
  // Wetted flow areas (approx - Table 2 gives Awk 115.2 cm2 -> 0.01152 m2)
  parameter Real A_w_k = 115.2e-4 "m^2 (Awk from table 115.2 cm^2)";
  //parameter Real A_w_h = A_w_k "use same order as Aw_k if not separately given"; //Actually, don't need it since we have d_h_h (ostensibly)
  
  parameter Real A_w_h = A_w_k/L_k*L_h;  
  parameter Real A_w_r1_wall = A_w_k/L_k*(L_r/2); //wetted area of wall, NOT wire mesh
  parameter Real A_w_r2_wall = A_w_k/L_k*(L_r/2); 
  parameter Real A_w_r1_mesh = 0.2019/2; //4*porosity*V_r1/d_h_r;
  parameter Real A_w_r2_mesh = 0.2019/2; //4*porosity*V_r2/d_h_r; //doesn't really affect the solution 

  
  // Expected solution (from linearised approximation)
  parameter Real omega = 2*pi*30.0; // 30 Hz
  parameter Real phi = -42.5*pi/180; //phase shift
  //parameter Real Vdot_max = 0.0080674; //graphically solved max from analytical expression                                          
  //parameter Real Vdot_max = abs(A_p*X_p*omega*sin(phi)/sin(atan((A_p*X_p*omega*sin(phi))/(A_p*X_p*omega*cos(phi)-(2*A_d-A_rod)*X_d*omega)))); 
  parameter Real Vdot_max = omega*sqrt( (A_p*X_p)^2. + ((2*A_d - A_rod)*X_d)^2 - 2*A_p*X_p*(2*A_d - A_rod)*X_d*sin(phi) ); //Equation from 1st Majidniya paper
  parameter Real u_h_max = Vdot_max/A_h;
  parameter Real u_k_max = Vdot_max/A_k;
  parameter Real u_r_max = Vdot_max/A_r;
  
  /*********/
  /*Thermal*/
  /*********/
  parameter Real cp = 5/2*R_gas "heat capacity of helium at const pressure";
  parameter Real cv = 3/2*R_gas "heat capacity of helium at const volume";
  

  // ------------------- STATES -------------------
  Real x_p(start=X_p*sin(phi), fixed=true) "power piston displacement (m)"; //Must start with correct phase, X_p*sin(omega*t + phi), t=0  //
  //Real x_p(start=X_p*sin(-pi/2)) "power piston displacement (m)"; //Try with "perfect" phase
  
  //Real v_p(start=0.0, fixed=false); 
  Real v_p(fixed=false); 
  //Real v_p(start=X_p*omega*cos(phi)) "power piston velocity (m/s)"; //Need a little kick? //Can use the same trick, analytical derivative

  //Real x_d(start=0.0, fixed=false) "displacer displacement (m)";
  Real x_d(fixed=false) "displacer displacement (m)";
  //Real v_d(start=X_d*omega) "displacer velocity (m/s)";
  //Real v_d(start=0.0, fixed=false);
  Real v_d(fixed=false);

  // ------------------- Intermediates (named to match paper) -------------------
  Real Vc "compression instantaneous volume";
  Real Ve "expansion instantaneous volume";
  Real Vd "displacer gas-spring instantaneous volume";
  Real Vb "buffer instantaneous volume";
  
  //Real p;  
  Real p_b;
  Real p_d;
  
  Real p_c(start=p_mean);
  Real p_k(start=p_mean);
  Real p_r1(start=p_mean);
  Real p_r2(start=p_mean); 
  Real p_h(start=p_mean);
  Real p_e(start=p_mean);
  
  Real dP_h(start=1.0, fixed=false);   
  Real dP_k(start=1.0, fixed=false);   
  Real dP_r1(start=1.0, fixed=false); 
  Real dP_r2(start=1.0, fixed=false); 
  Real deltaP(start=1.0, fixed=false) "pressure drop driving displacer (from heat exchanger/regenerator pressure drops)";

  // For pressure-drop velocity proxies
  //Real Vdot_c, Vdot_e;
  Real Vdot_flow(start = 0.01);
  Real u_h, u_k, u_r;    
   
  Real Re_h(start = 1000); 
  Real Re_k(start = 1000); 
  
  Real Re_r1(start = 1000);
  Real Re_r2(start = 1000);
  
  Real Cf_h "heater Darcy friction factor";
  Real Cf_k "cooler Darcy friction factor";
  Real Cf_r1  "regenerator Darcy friction factor";
  Real Cf_r2 ;  
    
  Real F_load; 
  Real P_out_mech;
  Real E_out_mech (start=0, fixed=true); //this is the integral of power, should definitely start at 0. 
  Real P_mech_avg; 

  //Dynamically calculated densities and viscosities following Majidniya's thesis
  Real rho_h(start=1,fixed=false);
  Real rho_k(start=1,fixed=false);
  //Real rho_r; 
  Real rho_r1(start=1,fixed=false); 
  Real rho_r2(start=1,fixed=false);
  
  Real mu_h(start=1e-6);
  Real mu_k(start=1e-6);
  //Real mu_r(start=1e-6);  
  Real mu_r1(start=1e-6);  
  Real mu_r2(start=1e-6);   

  
  Real V_avT (start=4e-7, fixed=false);     
  //parameter Real V_avT = A_p*C_c/T_k + A_d*C_e/T_h + V_h/T_h + V_k/T_k + V_r/T_r;     
  
  
  parameter Real u_eps = 1e-3;
    
  /* Heat transfer */   
  parameter Real T_wall_h = T_h_nominal;
  parameter Real T_wall_k = T_k_nominal;  
  parameter Real T_r_nominal_mean = (T_h_nominal-T_k_nominal)/log(T_h_nominal/T_k_nominal);
  parameter Real T_r1_nominal = 418.4;
  parameter Real T_r2_nominal = 662.7;
  parameter Real T_wall_r1 = T_r1_nominal;
  parameter Real T_wall_r2 = T_r2_nominal; 
  parameter Real T_mesh_r1_bc = T_r1_nominal;
  parameter Real T_mesh_r2_bc = T_r2_nominal;
  
  Real T_h(start = T_h_nominal);
  Real T_k(start = T_k_nominal);
  Real T_r1(start = T_r1_nominal);
  Real T_r2(start = T_r2_nominal);
  Real T_c(start = T_r_nominal_mean);
  Real T_e(start = T_r_nominal_mean);
  
  Real rho_c;
  Real rho_e;
  Real m_h;
  Real m_k;
  Real m_r1;
  Real m_r2;
  Real m_c; 
  Real m_e;
  
  Real mflow_kc;
  Real mflow_r1k;
  Real mflow_r2r1;
  Real mflow_he;
  Real mflow_r2h;
  Real mflow_r1r2;
  
  Real Pr_k; 
  Real Pr_h;
  Real Pr_r1;
  Real Pr_r2;
  
  Real HTC_k;
  Real HTC_h;
  Real HTC_r1;
  Real HTC_r2;
  
  Real Q_k;
  Real Q_h;
  Real Q_r1;
  Real Q_r2;
  
  Real W_cp "rate of work done on the power piston by the gas in the compression volume"; 
  Real W_ed "rate of work done on the displacer piston by the gas in the expansion volume";
  
  /* //Don't need explicit work terms for the buffer and displacer gas spring volumes since their thermodynamics are implicit
  Real W_bp;
  Real W_dd;
  */
  
  //Modelling the wire mesh explicitly 
  parameter Real rho_steel = 8000;
  parameter Real cp_steel = 500;
  parameter Real V_mesh = V_r*(1-porosity); 
  parameter Real m_mesh_1 = V_mesh*rho_steel/2; 
  parameter Real m_mesh_2 = V_mesh*rho_steel/2;
  Real T_mesh_r1(start=T_mesh_r1_bc); 
  Real T_mesh_r2(start=T_mesh_r2_bc);
  
//Now following the equation order proposed in Majidniya's follow-up 'performance' paper 
equation
    /************/
   /*Mechanical*/
  /************/
  
  // Net volumetric flow rate in the system 
  Vdot_flow = A_p*v_p - (2*A_d - A_rod)*v_d; //exact
  //Vdot_flow = der(Vc) - der(Ve); //numerical
  
  //Instantaneous velocities depending on XS area 
  u_h = Vdot_flow/A_h;
  u_k = Vdot_flow/A_k;
  u_r = Vdot_flow/A_r;
  
  //Viscosities 
  mu_k = 3.674*1e-7*T_k^0.7;
  //mu_r = 3.674*1e-7*T_r^0.7;
  mu_r1 = 3.674*1e-7*T_r1^0.7;
  mu_r2 = 3.674*1e-7*T_r2^0.7;
  mu_h = 3.674*1e-7*T_h^0.7;
  
  //Instantaneous volumes   
  Vc = A_p*(x_p + C_c) - (A_d - A_rod)*x_d;
  Ve = A_d*(x_d + C_e);
  
  Vb = V_B - A_p*x_p;   // buffer instantaneous volume 
  Vd = V_D - A_rod*x_d; // displacer gas spring instantaneous volume   
  
  
  //V_avT =  A_p*C_c/T_k + A_d*C_e/T_h + V_h/T_h + V_k/T_k + (V_r1+V_r2)/((T_r1 + T_r2)/2.);
  //V_avT =  A_p*C_c/T_k + A_d*C_e/T_h + V_h/T_h + V_k/T_k + V_r1/T_r1 + V_r2/T_r2;
  
  /*Static V_avT: only works without thermal coupling*/
  //V_avT =  A_p*C_c/T_k_nominal + A_d*C_e/T_h_nominal + V_h/T_h_nominal + V_k/T_k_nominal + V_r1/T_r1_nominal + V_r2/T_r2_nominal;
  //p_c = p_mean*(1 + (A_p*x_p - (A_d - A_rod)*x_d)/(T_k*V_avT) + (A_d*x_d)/(T_h*V_avT) )^(-1) ;  
  //p_c = m_total_gas*R_gas/V_avT*(1 + (A_p*x_p - (A_d - A_rod)*x_d)/(T_k*V_avT) + (A_d*x_d)/(T_h*V_avT) )^(-1) ;  
  
  /*/Dynamic V_avT*/
  V_avT = Vc/T_c + Ve/T_e + V_h/T_h + V_k/T_k + V_r1/T_r1 + V_r2/T_r2;  
  p_c = m_total_gas*R_gas/V_avT;   
  
  
  // Darcy pressure drops in the "bypass" region (heater, regenerator, cooler)    
  dP_h = 0.5 * rho_h * Cf_h * (L_h/d_h_h) * u_h * sqrt(u_h^2 + u_eps^2);
  dP_k = 0.5 * rho_k * Cf_k * (L_k/d_h_k) * u_k * sqrt(u_k^2 + u_eps^2);  
  dP_r1 = 0.5 * rho_r1 * Cf_r1 * ((L_r/2)/d_h_r) * u_r * sqrt(u_r^2 + u_eps^2);
  dP_r2 = 0.5 * rho_r2 * Cf_r2 * ((L_r/2)/d_h_r) * u_r * sqrt(u_r^2 + u_eps^2); //still assume the same friction factor for now
  
  
  //Friction factors 
  //Laminar, transition, turbulent
  
  if (Re_h < 2000) then
    Cf_h = 64/Re_h;    
  else    
    Cf_h = 0.316*Re_h^(-0.25); //Assume turbulent conditions for now
  end if;
  
  if (Re_k < 2000) then
    Cf_k = 64/Re_k;
  else     
    Cf_k = 0.316*Re_k^(-0.25); 
  end if;
  
  if (Re_r1 < 60) then 
    Cf_r1 = 4*10^(1.73 - 0.93 * log10(Re_r1));    
  elseif (Re_r1 < 1000) then
    Cf_r1 = 4*10^(0.714 - 0.365 * log10(Re_r1));
  else
    Cf_r1 = 4*10^(0.015 - 0.125 * log10(Re_r1));
  end if; 
  
  if (Re_r2 < 60) then 
    Cf_r2 = 4*10^(1.73 - 0.93 * log10(Re_r2));    
  elseif (Re_r1 < 1000) then
    Cf_r2 = 4*10^(0.714 - 0.365 * log10(Re_r2));
  else
    Cf_r2 = 4*10^(0.015 - 0.125 * log10(Re_r2));
  end if; 
  
  
  //Only turbulent/transition, avoids instability
  /*
  Cf_h = 0.316*Re_h^(-0.25); //Assume turbulent conditions for now
  Cf_k = 0.316*Re_k^(-0.25); 
  Cf_r1 = 4*10^(0.714 - 0.365 * log10(Re_r1));
  Cf_r2 = 4*10^(0.714 - 0.365 * log10(Re_r2));
  */
  
  //Reynolds numbers      
  //--compute according to estimated max velocity, as per first Majidniya paper  
  
  Re_h = rho_h*abs(u_h_max)*d_h_h/mu_h;
  Re_k = rho_k*abs(u_k_max)*d_h_k/mu_k;
  Re_r1 = rho_r1*abs(u_r_max)*d_h_r/mu_r1;
  Re_r2 = rho_r2*abs(u_r_max)*d_h_r/mu_r2;
  
  
  //--compute according to instantaneous local velocity, as per Majidniya's thesis and 2nd paper 
  // Can lead to instability when u_h, u_k, u_r are close to 0    
  // Can also cause overdamping 
  /*
  Re_h = max(1,rho_h*sqrt(u_h^2 + u_eps^2)*d_h_h/mu_h);
  Re_k = max(1,rho_k*sqrt(u_k^2 + u_eps^2)*d_h_k/mu_k);
  Re_r1 = max(1,rho_r1*sqrt(u_r^2 + u_eps^2)*d_h_r/mu_r1);
  Re_r2 = max(1,rho_r2*sqrt(u_r^2 + u_eps^2)*d_h_r/mu_r2);
  */
  
  //Pressure states
  p_k = p_c + dP_k/2.;
  p_r1 = p_k + dP_k/2. + dP_r1/2.;
  p_r2 = p_r1 + dP_r1/2. + dP_r2/2.;
  p_h = p_r2 + dP_r2/2. + dP_h/2;
  p_e = p_h + dP_h/2.;  
  
  //Densities   
  rho_k = 48.14*p_k*1e-5/(T_k*(1 + 0.4446*p_k*1e-5*T_k^(-1.2)));
  rho_r1 = 48.14*p_r1*1e-5/(T_r1*(1 + 0.4446*p_r1*1e-5*T_r1^(-1.2)));
  rho_r2 = 48.14*p_r2*1e-5/(T_r2*(1 + 0.4446*p_r2*1e-5*T_r2^(-1.2)));
  rho_h = 48.14*p_h*1e-5/(T_h*(1 + 0.4446*p_h*1e-5*T_h^(-1.2)));
  rho_c = 48.14*p_c*1e-5/(T_c*(1 + 0.4446*p_c*1e-5*T_c^(-1.2))); 
  rho_e = 48.14*p_e*1e-5/(T_e*(1 + 0.4446*p_e*1e-5*T_e^(-1.2)));   
  
  //Gas springs  
  p_b = p_mean*(V_B/Vb)^gamma; 
  p_d = p_mean*(V_D/Vd)^gamma;
  
  //Total pressure drop between hot/cold sides (compression and expansion spaces)  
  deltaP = dP_h + dP_r1 + dP_r2 + dP_k;
    
  F_load = B_load*v_p; //value from second Majidniya pub, with thermal coupling
  //F_load = C_load*v_p; //value from first Majidniya pub, no thermal coupling
  
  P_out_mech = F_load*v_p; 
  der(E_out_mech) = P_out_mech;  
  // average power computed externally (user can post-process or use E_elec_int/time)
  //P_mech_avg = E_out_mech/ max(time, 1e-5); 
  P_mech_avg = ( E_out_mech - delay(E_out_mech, 0.1) ) / 0.1;//moving average with 0.1s window

  //Piston mechanics. Can introduce fudge factors (spring consts, damping coeffs) to tune frequency & amplitude if need be.
  der(x_p) = v_p;  
  der(v_p) = ( A_p*(p_c - p_b)  - F_load )/m_p; //- C_p*v_p - K_p*x_p

  der(x_d) = v_d;  
  der(v_d) = ( A_d*deltaP + A_rod*(p_c - p_d) - C_d*v_d - K_d*x_d )/ m_d; 
  
    /*********/
   /*Thermal*/
  /*********/
  
  //masses
  m_h = rho_h*V_h; 
  m_k = rho_k*V_k;
  m_r1 = rho_r1*V_r1;
  m_r2 = rho_r2*V_r2;       
  m_c = rho_c*Vc;  
  m_e = rho_e*Ve;
  
  //Mass flow rates
  //original formulation from Majidiya's second paper & thesis
  mflow_kc = der(m_c);
  mflow_r1k = mflow_kc + der(m_k);
  mflow_r2r1 = mflow_r1k + der(m_r1);
  mflow_he = der(m_e); 
  mflow_r2h = mflow_he + der(m_h);
  mflow_r1r2 = mflow_r2h + der(m_r2);
  
  // Version that enforces flow continuity
  /*
  mflow_kc   = der(m_c);
  mflow_r1k  = -(-mflow_kc  - der(m_k));
  mflow_r1r2 = -mflow_r1k - der(m_r1);
  mflow_r2h  = mflow_r1r2 - der(m_r2);
  mflow_he   = mflow_r2h - der(m_h);
  mflow_r2r1 = -mflow_r1r2;
  */
  
  //Wall heat transfer coeffs
  Pr_h = 0.7117*T_h^(-1.*(0.01 - 1.42E-4*p_h*1e-5))/(1 + 1.123E-3*p_h*1e-5);
  Pr_k = 0.7117*T_k^(-1.*(0.01 - 1.42E-4*p_k*1e-5))/(1 + 1.123E-3*p_k*1e-5);
  Pr_r1 = 0.7117*T_r1^(-1.*(0.01 - 1.42E-4*p_r1*1e-5))/(1 + 1.123E-3*p_r1*1e-5);
  Pr_r2 = 0.7117*T_r2^(-1.*(0.01 - 1.42E-4*p_r2*1e-5))/(1 + 1.123E-3*p_r2*1e-5);  
  HTC_h = Cf_h/4*(Re_h*mu_h*cp)/(2*d_h_h*Pr_h);
  HTC_k = Cf_k/4*(Re_k*mu_k*cp)/(2*d_h_k*Pr_k);  
  HTC_r1 = Cf_r1/4*(Re_r1*mu_r1*cp)/(2*d_h_r*Pr_r1);
  HTC_r2 = Cf_r2/4*(Re_r2*mu_r2*cp)/(2*d_h_r*Pr_r2);
  
  //Wall heat transfer
  Q_k = HTC_k*A_w_k*(T_wall_k - T_k);
  Q_h = HTC_h*A_w_h*(T_wall_h - T_h);
  
  Q_r1 = HTC_r1*A_w_r1_wall*(T_wall_r1 - T_r1);
  Q_r2 = HTC_r2*A_w_r2_wall*(T_wall_r2 - T_r2);
  
  //Regenerative heat transfer (to wire mesh--hence the use of d_h_r)
  /*
  Q_r1 = HTC_r1*A_w_r1_mesh*(T_mesh_r1_bc - T_r1); 
  Q_r2 = HTC_r2*A_w_r2_mesh*(T_mesh_r2_bc - T_r2);
  */
  //Q_r1 = HTC_r1*A_w_r1_mesh*(T_mesh_r1 - T_r1); 
  //Q_r2 = HTC_r2*A_w_r2_mesh*(T_mesh_r2 - T_r2);
  
  //NOTE: The instantaneous mesh temperature is not resolved, hence this model cannot be used for cold start-up.
  
  //Work done by the gas. There is no work done in const-volume spaces. (h,k,r1,r2)
  W_cp = A_p*(p_c - p_b)*v_p;   
  //W_ed = A_d*(p_e - p_c)*v_d;
  //W_cp = A_p*p_c*v_p; 
  W_ed = A_d*p_e*v_d;
  
  //W_cp = p_c*der(Vc);  //numerical
  //W_ed = p_e*der(Ve);  
  //W_cp = p_c*(A_p*v_p - (A_d - A_rod)*v_d); //analytic
  //W_ed = p_e*A_d*v_d;
  
  //Energy balance. Energy advection depends on the instantaneous flow direction   
  
  m_c*cv*der(T_c) = mflow_kc*cp*(if mflow_kc > 0 then T_k else T_c) - W_cp - der(m_c)*cv*T_c;
  
  m_k*cv*der(T_k) = Q_k + (-mflow_kc)*cp*(if (-mflow_kc) > 0 then T_c else T_k) + mflow_r1k*cp*(if mflow_r1k > 0 then T_r1 else T_k) - der(m_k)*cv*T_k ;
  
  m_r1*cv*der(T_r1) = Q_r1 + (-mflow_r1k)*cp*(if (-mflow_r1k) > 0 then T_k else T_r1) + mflow_r2r1*cp*(if mflow_r2r1 > 0 then T_r2 else T_r1) - der(m_r1)*cv*T_r1 ;
  
  m_r2*cv*der(T_r2) = Q_r2 + mflow_r1r2*cp*(if mflow_r1r2 > 0 then T_r1 else T_r2) + (-mflow_r2h)*cp*(if (-mflow_r2h) > 0 then T_h else T_r2) - der(m_r2)*cv*T_r2 ;
  
  m_h*cv*der(T_h) = Q_h + mflow_r2h*cp*(if mflow_r2h > 0 then T_r2 else T_h)  + (-mflow_he)*cp*(if (-mflow_he) > 0 then T_e else T_h) - der(m_h)*cv*T_h ;
  
  m_e*cv*der(T_e) = mflow_he*cp*(if mflow_he > 0 then T_h else T_e) - W_ed - der(m_e)*cv*T_e; 
  
  
  //New equation for mesh temperature and modified regenerator equations to exchange with its explicit representation. 
  m_mesh_1*cp_steel*der(T_mesh_r1) = -Q_r1;
  m_mesh_2*cp_steel*der(T_mesh_r2) = -Q_r2;
   
  //--------------------
  //Alt fixed values in case no thermal coupling is desired (must also comment out heat balance equations). 
  /*
  T_h = T_h_nominal;
  T_r1 = T_r1_nominal;
  T_r2 = T_r2_nominal;
  T_k = T_k_nominal;
  T_c = T_k_nominal; 
  T_e = T_h_nominal; 
  */
  
  //annotation (experiment(StartTime=0, StopTime=0.15, Tolerance=1e-5, Interval=1e-5));
end FPSE_Nonlinear_Majidniya_WithHeatTransfer;
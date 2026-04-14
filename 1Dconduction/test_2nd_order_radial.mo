model test_2nd_order_radial
// -------------------------------
  // Geometry and discretization
  // -------------------------------
  parameter Integer N = 100 "Number of radial cells";
  parameter Real ri = 0.005 "Inner radius [m]";
  parameter Real ro = 0.015 "Outer radius [m]";
  parameter Real L  = 1.0   "Axial length [m]";
  parameter Real dr = (ro - ri)/N;

  // -------------------------------
  // Material properties (U-10Mo)
  // -------------------------------
  parameter Real rho = 17000 "Density [kg/m3]";
  parameter Real cp  = 150   "Specific heat [J/kg-K]";
  parameter Real k   = 15    "Thermal conductivity [W/m-K]";

  // -------------------------------
  // Heat generation and BCs
  // -------------------------------
  parameter Real qppp = 1e7 "Volumetric heat generation [W/m3]";
  parameter Real T_outer = 600 "Outer wall temperature [K]";

  // -------------------------------
  // State variables
  // -------------------------------
  Real T[N](start=fill(T_outer, N));

  // -------------------------------
  // Geometry (PARAMETERS, not variables)
  // -------------------------------
  parameter Real r_f[N+1] = {ri + (i-1)*dr for i in 1:N+1};
  parameter Real r_c[N]   = {ri + (i-0.5)*dr for i in 1:N};

  parameter Real V[N] =
    {Modelica.Constants.pi*L*(r_f[i+1]^2 - r_f[i]^2) for i in 1:N};

  parameter Real A_f[N+1] =
    {2*Modelica.Constants.pi*L*r_f[i] for i in 1:N+1};

initial equation
  for i in 1:N loop 
    der(T[i]) = 0;
  end for;
equation
  // -------------------------------
  // Inner boundary: adiabatic
  // -------------------------------
  
  der(T[1]) =
    ( k*A_f[2]*(T[2]-T[1])/dr
    + qppp*V[1] )
    / (rho*cp*V[1]);
  /*
  der(T[1]) =
  ( k*A_f[2]*(T[2]-T[1])/dr
  - k*A_f[1]*(T[1]-T[1])/dr
  + qppp*V[1] )
  / (rho*cp*V[1]);
  */
  // -------------------------------
  // Interior cells
  // -------------------------------
  
  for i in 2:N-1 loop
    der(T[i]) =
      ( k*A_f[i+1]*(T[i+1]-T[i])/dr
      - k*A_f[i]*(T[i]-T[i-1])/dr
      + qppp*V[i] )
      / (rho*cp*V[i]);
  end for;

  // -------------------------------
  // Outer boundary: prescribed T
  // -------------------------------
  /*
  der(T[N]) =
    ( k*A_f[N]*(T[N-1]-T[N])/dr
    + k*A_f[N+1]*(T_outer-T[N])/dr
    + qppp*V[N] )
    / (rho*cp*V[N]);
    */
  der(T[N]) =
  ( k*A_f[N]*(T[N-1]-T[N])/dr
  + k*A_f[N+1]*(T_outer-T[N])/(dr/2.)
  + qppp*V[N] )
  / (rho*cp*V[N]);  

end test_2nd_order_radial;
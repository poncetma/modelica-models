model FreePistonOneNode_test
  "Single-node, KRUSTY-flavoured thermoacoustic toy model.
   - Phenomenological pressure drive p = p0 + beta * sin(omegaDrive*time + phi) * x_p
     (beta and phi tuneable to produce a positive average p·V work)
   - Single piston mass with mechanical damping and an electrical viscous damping B_elec
   - Built-in RMS estimator for v_p (first-order filter) so you can compute suggested B_elec
   - Parameters chosen to be in the ASC/ KRUSTY ballpark (see comments)
   - Use this to tune frequency/amplitude/power quickly; once tuned you can
     move to the more detailed 2-node model and map parameters across.
  "

  // ---------------- PARAMETERS (KRUSTY-flavoured defaults) ----------------
  // Phenomenological single-node volume (sum of typical ASC dead + swept volumes)
  parameter Real V0 = 1.30e-4 "Nominal total gas volume (m^3) ~ ASC dead vol sum (starting guess)";

  // Physical reference values (for documentation / mapping)
  parameter Real p_ref = 3.52e6 "Mean charge pressure (Pa) used as context (not used directly)";
  parameter Real T_hot = 1123.0 "Hot boundary temperature (K)";
  parameter Real T_cold = 363.0 "Cold boundary temperature (K)";

  // Piston geometry / mechanics (ASC-style numbers)
  parameter Real A_p = 1.539e-4 "Power piston area (m^2) (ASC estimate)";
  // Choose piston mass so mechanical natural frequency is in the target range
  // For a target f ≈ 100 Hz with a gas+spring effective k ≈ 2e5 N/m -> m ≈ k/ω^2
  // I set a starting piston mass somewhat heavy to match earlier tuning suggestions.
  parameter Real m_p = 0.49 "Piston mass (kg) - tune this to move natural freq";

  // mechanical linear spring (represents mechanical centring + stiffness)
  // choose k_p to give a nominal undamped mechanical natural frequency near the target
  parameter Real k_p = 1.0e5 "Piston spring stiffness (N/m) - rough starting value";

  // mechanical viscous damping (non-electrical losses) - keep small initially
  parameter Real c_mech = 2.0 "Mechanical viscous damping (N·s/m)";

  // ---------------- THERMOACOUSTIC DRIVE (phenomenological) ----------------
  // beta controls the amplitude of the pressure drive due to timed heat input: larger -> stronger drive
  // Start with a modest value and increase if oscillations die (try 1e4..1e6)
  parameter Real beta = 6.0e4 "Thermoacoustic drive coefficient (Pa) - tune this";

  // drive frequency (rad/s). For KRUSTY target f ≈ 100 Hz -> omega ≈ 2π·100
  parameter Real f_target = 100.0 "Target mechanical frequency (Hz)";
  parameter Real omegaDrive = 2*Modelica.Constants.pi * f_target "Drive angular frequency (rad/s)";

  // drive phase (radians) — set so that pressure leads or lags displacement to produce positive p·V
  // typical useful range for thermoacoustic effective phasing is ~ 45–90° (0.785..1.57 rad)
  parameter Real phi = 1 "Drive phase (rad) - tune this for best p·V work";

  // ---------------- ALTERNATOR (simple viscous) ----------------
  // Use viscous electrical damping (easy tuning)
  parameter Real B_elec(start=0.0) "Electrical equivalent damping (N·s/m) - tune from v_rms & P_target";

  // ---------------- RMS estimator (filter) ----------------
  // Implements a first-order low-pass for v^2 to produce an RMS estimate:
  // d(s)/dt = (v_p^2 - s)/tau_rms  => v_rms = sqrt(s)
  parameter Real tau_rms = 0.01 "Filter time constant for RMS estimator (s) — choose a few cycles length (0.01s => ~1 cycle at 100Hz)";

  // ---------------- STATES ----------------
  Real x_p(start=0.0) "Piston displacement (m). Positive direction arbitrary but consistent";
  Real v_p(start=0.1) "Piston velocity (m/s) - nonzero start helps get motion";
  Real s_v2(start=( (2.0)^2 )) "Filter state approximating mean(v^2) (initialized to expected v_rms^2)";

  // bookkeeping
  Real v_rms "estimated RMS of piston velocity (m/s)";
  Real p "instantaneous pressure (Pa)";
  Real V "instantaneous single-node volume (m^3)";
  Real P_pv_inst "instantaneous p·Vdot power (W) (for diagnostics)";
  Real P_elec_inst "instantaneous electrical power extracted (W)";
  Real P_elec_avg "time-averaged electrical power (J/s) approximated via E_int/time";
  Real E_elec_int(start=0.0) "integral of electrical energy (J)";

  // optional target power convenience parameter (set by you)
  parameter Real P_target = 80.0 "Target electrical power (W) - used to compute suggested B if desired";

equation
  // ---------- volume & kinematics ----------
  // Single-node model: choose a mapping from piston displacement to effective gas volume.
  // We keep it linear and simple: V = V0 - A_eff * x_p.
  // Choose A_eff to scale realistic volume swings (we use piston area A_p as A_eff).
  V = V0 - A_p * x_p;
  // instantaneous volume rate
  // Vdot = - A_p * v_p  (not needed explicitly but used in P_pv_inst)
  P_pv_inst = p * ( - A_p * v_p ); // p * Vdot  (power delivered by gas to piston; positive if gas does work on piston)

  // ---------- phenomenological pressure drive ----------
  // p = p_mean + beta * sin(omegaDrive * t + phi) * x_p
  // - This is a compact way to introduce a timed pressure term whose amplitude scales with displacement.
  // - The product sin(...) * x_p produces the p-x loop area needed for net p·V work.
  p = p_ref + beta * sin(omegaDrive * time + phi) * x_p;

  // ---------- piston dynamics (including simple electrical damping) ----------
  der(x_p) = v_p;
  der(v_p) = ( p * A_p  // gas force on piston
               - c_mech * v_p
               - k_p * x_p
               - B_elec * v_p  // electrical reaction treated as viscous damping
             ) / m_p;

  // ---------- electrical power bookkeeping ----------
  P_elec_inst = B_elec * v_p * v_p; // instantaneous mechanical power removed
  der(E_elec_int) = P_elec_inst;
  P_elec_avg = E_elec_int / max(time, 1e-6);

  // ---------- RMS estimator (first-order filter of v^2) ----------
  der(s_v2) = ( v_p*v_p - s_v2 ) / tau_rms;
  v_rms = sqrt( max(s_v2, 1e-12) );

  // ---------- convenience computed suggestion for B (not used by model unless you set it) --
  // Recommended B to achieve P_target assuming P_elec_avg ≈ B * v_rms^2 (approx)
  // (user may read this and then set B_elec accordingly)
  // suggested_B = P_target / (v_rms^2)  -- leave as algebraic comment/guide only
  // (we do not assign it into B_elec automatically so you can ramp it manually)

  annotation (experiment(StartTime=0.0, StopTime=10, Interval=1e-4, Tolerance=1e-6));
end FreePistonOneNode_test;
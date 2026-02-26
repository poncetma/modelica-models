model Stirling_controlledheatremoval

input Real power_setpoint;
parameter Real cycle_efficiency = 0.2;

Real Q_elec;
Real Q_waste;

equation
Q_elec = power_setpoint*cycle_efficiency;
Q_waste = power_setpoint*(1-cycle_efficiency);


end Stirling_controlledheatremoval;
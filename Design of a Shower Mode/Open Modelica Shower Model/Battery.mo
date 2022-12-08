model Battery
  // Battery Parameters
  parameter Real eff_bat = 0.95;                 //One trip Battery efficiency
  parameter Real battery_cap = 13500*3.6e3;      //J, Capacity of battery
  parameter Real max_power = 1200;               //W, Maximum power transfer rate
  parameter Real max_SOC = 0.9;                  //%, Maximum state of charge
  parameter Real min_SOC = 0.1;                  //%, Minimum state of charge
  
  // Battery Variables
  Real battery_energy(start=0.5*battery_cap);
  Real SOC;

  // Inputs - Outpus
  Modelica.Blocks.Interfaces.RealInput power_in annotation(
    Placement(visible = true, transformation(origin = {-80, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-72, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput power_out annotation(
    Placement(visible = true, transformation(origin = {70, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {86, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput discharge_command annotation(
    Placement(visible = true, transformation(origin = {-6, -74}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {2, 62}, extent = {{20, -20}, {-20, 20}}, rotation = 90)));

  // Auxiliary variables (To avoid Open Modelica and FMU crashing)
  Real power_out_var;
  Real power_in_var;

  // Auxiliary variables (To avoid Open Modelica and FMU crashing)
  Real time_counter_up(start=0);
  Real time_counter_down(start=0);
  Real time_counter_value_up(start=0);
  Real time_counter_value_down(start=0);
  parameter Real cool_off_period = 1*60; 
  Boolean is_cooling_off_period_happening_up(start=false);
  Boolean is_cooling_off_period_happening_down(start=false);

algorithm

  // Battery overcharge and undercharge protection 
  when SOC <= min_SOC then
    // When battery is empty, raise flag
    is_cooling_off_period_happening_down:=true;
  end when;
  when SOC >= max_SOC then
    // When battery is full, raise flag
    is_cooling_off_period_happening_up:=true;
  end when;
  
  if is_cooling_off_period_happening_down then 
    // If flag is raised because battery is empty, start counter and stop battery from draining 
    time_counter_value_down:=1;
    power_out_var:=0;
  else 
    // Stop counter and allow battery to drain
    time_counter_value_down:=0;
    if SOC >= min_SOC then
      power_out_var:=discharge_command*max_power;
    end if;
  end if;
  
  if is_cooling_off_period_happening_up then 
    // If flag is raised because battery is full, start counter and stop battery from charging
    time_counter_value_up:=1;
    power_in_var:=0;
  else 
    // Stop conter and allow battery to charge
    time_counter_value_up:=0;
    if SOC <= max_SOC then
      power_in_var:=power_in;
    end if;
  end if; 
  
  when (time_counter_down > cool_off_period) then
    // When the counter is above treshold, lower flag
    is_cooling_off_period_happening_down:=false;
  end when; 
  when (time_counter_up > cool_off_period) then
    // When the counter is above treshold, lower flag
    is_cooling_off_period_happening_up:=false;
  end when;  
  
equation

  // Energy balance
  der(battery_energy) = eff_bat *  power_in_var - power_out_var / eff_bat;
  // SOC calculation
  SOC=battery_energy/battery_cap;
  
  // Battery overcharge and undercharge protection 
  der(time_counter_up)=time_counter_value_up;
  when not is_cooling_off_period_happening_up then
    reinit(time_counter_up,0);
  end when;  
  der(time_counter_down)=time_counter_value_down;
  when not is_cooling_off_period_happening_down then
    reinit(time_counter_down,0);
  end when;  

algorithm
  // Output
  power_out:=power_out_var;
  annotation(
    uses(Modelica(version = "4.0.0")),
    Icon(graphics = {Rectangle(lineColor = {96, 160, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.HorizontalCylinder, extent = {{-78, 60}, {78, -60}})}));
end Battery;

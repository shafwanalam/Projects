model tank

  //#################### Tank parameters
  parameter Real tank_mass = 2000;        //kg
  parameter Real max_capacity = 0.9;      //%
  parameter Real min_capacity = 0.1;      //%
  parameter Real specific_heat = 4182;    // J/Kg.degC
  parameter Real T_ini = 20;              // degC
  parameter Real T_fin = 60;
  parameter Real T_op = 32;
  
  //#################### WATER LEVEL CONTROL VARIABLES
  Real m_in_var, m_out_var;
  Real water_mass(start = tank_mass*0.5);
  Real water_temp(start = 32);
  Boolean is_waiting_time_water_happening_min(start = false);
  Boolean is_waiting_time_water_happening_max(start = false);
  Real time_counter_value_water_min(start = 0);
  Real time_counter_value_water_max(start = 0);
  Real time_counter_water_min(start = 0);
  Real time_counter_water_max(start = 0);
  parameter Real waiting_period = 10;
  
  //#################### INPUT AND OUTPUT
  Modelica.Blocks.Interfaces.RealInput m_in annotation(
    Placement(visible = true, transformation(origin = {-62, 50}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-88, 74}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput m_out annotation(
    Placement(visible = true, transformation(origin = {-62, -28}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-94, -68}, extent = {{20, -20}, {-20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput heat_in annotation(
    Placement(visible = true, transformation(origin = {32, 48}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {24, 84}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealInput desired_water_flow annotation(
    Placement(visible = true, transformation(origin = {54, -32}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {82, -82}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));

algorithm

  //#################### WATER LEVEL PROTECTION
  when water_mass <= tank_mass * min_capacity then
  // When tank is empty, raise flag
    is_waiting_time_water_happening_min := true;
  end when;
  when water_mass >= tank_mass * max_capacity then
  // When tank is full, raise flag
    is_waiting_time_water_happening_max := true;
  end when;
  if is_waiting_time_water_happening_min then
  // If flag is raised, start counter and stop output to impede tank from further draining
    time_counter_value_water_min := 1;
    m_out_var := 0;
  else
  // If flag is down, stop counter and allow tank to drain
    time_counter_value_water_min := 0;
    if water_mass >= tank_mass * min_capacity then
      m_out_var := desired_water_flow;
    end if;
  end if;
  if is_waiting_time_water_happening_max then
  // If flag is raised, start counter and stop input to impede tank from overflowing
    time_counter_value_water_max := 1;
    m_in_var := 0;
  else
  // If flag is down, stop counter and allow tank to fill
    time_counter_value_water_max := 0;
    if water_mass <= tank_mass * max_capacity then
      m_in_var := m_in;
    end if;
  end if;
  when time_counter_water_min > waiting_period then
  // When the counter is above treshold, lower flag.
    is_waiting_time_water_happening_min := false;
  end when;
  when time_counter_water_max > waiting_period then
  // When the counter is above treshold, lower flag.
    is_waiting_time_water_happening_max := false;
  end when;

equation

  //#################### MASS BALANCE
  der(water_mass) = m_in_var - m_out_var;
  //#################### ENERGY BALANCE
  specific_heat * der(water_temp * water_mass) = (m_in_var * T_ini - m_out_var * water_temp) * specific_heat + heat_in;
  //#################### WATER LEVEL PROTECTION
  der(time_counter_water_min) = time_counter_value_water_min;
  when not is_waiting_time_water_happening_min then
    reinit(time_counter_water_min, 0);
  end when;
  der(time_counter_water_max) = time_counter_value_water_max;
  when not is_waiting_time_water_happening_max then
    reinit(time_counter_water_max, 0);
  end when;

algorithm
  // Output
  m_out := m_out_var;
  annotation(
    uses(Modelica(version = "4.0.0")),
    Icon(graphics = {Rectangle( fillColor = {170, 170, 255}, fillPattern = FillPattern.Solid, extent = {{-96, 96}, {96, -96}}), Text(origin = {-59, 74}, extent = {{-55, 8}, {55, -8}}, textString = "IN"), Text(origin = {23, 56}, extent = {{-19, 10}, {19, -10}}, textString = "Heat"), Polygon(origin = {0, -41}, fillColor = {0, 255, 255}, fillPattern = FillPattern.Backward, points = {{-96, 41}, {-62, 31}, {-32, 41}, {2, 33}, {32, 43}, {64, 33}, {96, 45}, {96, -55}, {-96, -55}, {-96, 1}, {-96, 41}}), Text(origin = {-56, -68}, extent = {{-16, 8}, {16, -8}}, textString = "OUT")}));
end tank;

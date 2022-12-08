model pv_bat_comb
  pv_power pv_power1 annotation(
    Placement(visible = true, transformation(origin = {-32, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Battery battery1 annotation(
    Placement(visible = true, transformation(origin = {20, 48}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Solar123 solar123 annotation(
    Placement(visible = true, transformation(origin = {-74, 74}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  heat_pump heat_pump1 annotation(
    Placement(visible = true, transformation(origin = {70, 18}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  tank tank1 annotation(
    Placement(visible = true, transformation(origin = {36, -12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Step Battery_control(startTime = 10000) annotation(
    Placement(visible = true, transformation(origin = {68, 68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Step Tank_control(height = 8 * 50 / 600, startTime = 300) annotation(
    Placement(visible = true, transformation(origin = {-50, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  shower shower1 annotation(
    Placement(visible = true, transformation(origin = {-34, -86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.Constant Shower_control(k = 1)  annotation(
    Placement(visible = true, transformation(origin = {-50, -44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(pv_power1.power_out, battery1.power_in) annotation(
    Line(points = {{-22, 48}, {13, 48}}, color = {255, 0, 0}));
  connect(solar123.out_value, pv_power1.meteo_data) annotation(
    Line(points = {{-66, 75}, {-54, 75}, {-54, 48}, {-42, 48}}, color = {255, 0, 0}));
  connect(battery1.power_out, heat_pump1.input_power_heatpump) annotation(
    Line(points = {{29, 48}, {61, 48}, {61, 18}}, color = {255, 0, 0}));
  connect(heat_pump1.output_power_heatpump, tank1.heat_in) annotation(
    Line(points = {{80, 18}, {38, 18}, {38, -4}}, color = {255, 0, 0}));
/*connect(pulse1.y, tank1.m_in) annotation(
    Line(points = {{-61, 20}, {54, 20}}, color = {0, 0, 127}));*/
  connect(Battery_control.y, battery1.discharge_command) annotation(
    Line(points = {{79, 68}, {79, 68.5}, {20, 68.5}, {20, 54}}, color = {0, 0, 127}, pattern = LinePattern.Dot));
  connect(Tank_control.y, tank1.m_in) annotation(
    Line(points = {{-39, 12}, {-0.5, 12}, {-0.5, -5}, {27, -5}}, color = {0, 0, 127}));
  connect(tank1.m_out, shower1.flow_water) annotation(
    Line(points = {{27, -19}, {-70, -19}, {-70, -86}, {-42, -86}}, color = {0, 0, 127}));
  connect(shower1.desired_water_flow, tank1.desired_water_flow) annotation(
    Line(points = {{-24.5, -86}, {92, -86}, {92, -20}, {44, -20}}, color = {0, 0, 127}, pattern = LinePattern.Dot));
  connect(Shower_control.y, shower1.shower_command) annotation(
    Line(points = {{-39, -44}, {-34, -44}, {-34, -78}}, color = {0, 0, 127}, pattern = LinePattern.Dash));
protected
  annotation(
    uses(Modelica(version = "4.0.0")),
    Diagram(graphics = {Polygon(points = {{56, 52}, {56, 52}, {56, 52}})}));
end pv_bat_comb;

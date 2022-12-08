model heat_pump

  parameter Real efficiency_heater = 4;

  Modelica.Blocks.Interfaces.RealInput input_power_heatpump annotation(
    Placement(visible = true, transformation(origin = {-48, 32}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-94, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput output_power_heatpump annotation(
    Placement(visible = true, transformation(origin = {70, 36}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

equation
  
  output_power_heatpump = efficiency_heater * input_power_heatpump;

annotation(
    uses(Modelica(version = "3.2.3")),
    Icon(graphics = {Ellipse(origin = {0, 1}, fillColor = {226, 167, 218}, fillPattern = FillPattern.Solid, extent = {{-96, 95}, {96, -95}}), Text(origin = {-50, -1}, extent = {{-20, 11}, {20, -11}}, textString = "Power in"), Text(origin = {65, -3}, extent = {{-19, 11}, {19, -11}}, textString = "Heat out")}));
end heat_pump;

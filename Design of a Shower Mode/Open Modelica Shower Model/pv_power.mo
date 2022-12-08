model pv_power

  parameter Real efficiency = 0.2;
  parameter Real area = 6;
  Modelica.Blocks.Interfaces.RealOutput power_out annotation(
    Placement(visible = true, transformation(origin = {70, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput meteo_data annotation(
    Placement(visible = true, transformation(origin = {-62, 6}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-96, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));

equation
power_out = efficiency * area * meteo_data;

annotation(
    uses(Modelica(version = "4.0.0")),
    Icon(graphics = {Polygon(origin = {-64, -48}, fillColor = {51, 130, 247}, fillPattern = FillPattern.Cross, points = {{-42, -48}, {-12, 48}, {42, 48}, {12, -48}, {-6, -48}, {-42, -48}}), Polygon(origin = {-4, -48}, fillColor = {51, 130, 247}, fillPattern = FillPattern.Cross, points = {{-42, -48}, {-12, 48}, {42, 48}, {12, -48}, {-6, -48}, {-42, -48}}), Polygon(origin = {56, -48}, fillColor = {51, 130, 247}, fillPattern = FillPattern.Cross, points = {{-42, -48}, {-12, 48}, {42, 48}, {12, -48}, {-6, -48}, {-42, -48}}), Polygon(origin = {-62, 50}, fillColor = {51, 130, 247}, fillPattern = FillPattern.Cross, points = {{-42, -48}, {-12, 48}, {42, 48}, {12, -48}, {-6, -48}, {-42, -48}}), Polygon(origin = {-2, 50}, fillColor = {51, 130, 247}, fillPattern = FillPattern.Cross, points = {{-42, -48}, {-12, 48}, {42, 48}, {12, -48}, {-6, -48}, {-42, -48}}), Polygon(origin = {58, 50}, fillColor = {51, 130, 247}, fillPattern = FillPattern.Cross, points = {{-42, -48}, {-12, 48}, {42, 48}, {12, -48}, {-6, -48}, {-42, -48}})}));
end pv_power;

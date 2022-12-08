model shower
  //Input: flow of water at 32 C
  //Input: command to start a shower
  //Output: number of shower availables (showers remaining)
  //Output: flow of water required
  //Script will always prioritize shower from number 1 to number n_shower
  parameter Real n_shower = 5;
  parameter Real water_flow_per_shower = 50 / (10 * 60);
  parameter Real time_per_shower = 10 * 60;
  parameter Real t_p = 50;
 
  Boolean is_shower_active_1(start=false);
  Boolean is_shower_active_2(start=false);
  Boolean is_shower_active_3(start=false);
  Boolean is_shower_active_4(start=false);
  Boolean is_shower_active_5(start=false);

  Real chattering_counter(start=0);
  Real chattering_counter_value(start=0);

  Real shower_time_counter_1(start=0);
  Real shower_time_counter_2(start=0);
  Real shower_time_counter_3(start=0);
  Real shower_time_counter_4(start=0);
  Real shower_time_counter_5(start=0);

  Real shower_time_counter_value_1(start=0);
  Real shower_time_counter_value_2(start=0);
  Real shower_time_counter_value_3(start=0);
  Real shower_time_counter_value_4(start=0);
  Real shower_time_counter_value_5(start=0);

//Real available_flow;
  Real count(start=0);
  
  Boolean auxiliar(start=false);
  
  Modelica.Blocks.Interfaces.RealInput flow_water annotation(
    Placement(visible = true, transformation(origin = {-74, 14}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {-84, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealOutput desired_water_flow(start=0) annotation(
    Placement(visible = true, transformation(origin = {80, 40}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {95, 1}, extent = {{-21, -21}, {21, 21}}, rotation = 0)));
  Modelica.Blocks.Interfaces.RealInput shower_command annotation(
    Placement(visible = true, transformation(origin = {-8, 70}, extent = {{-20, -20}, {20, 20}}, rotation = 0), iconTransformation(origin = {0, 82}, extent = {{-20, -20}, {20, 20}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealOutput showers_remaining(start=n_shower) annotation(
    Placement(visible = true, transformation(origin = {6, -78}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {0, -92}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
  Modelica.Blocks.Interfaces.RealOutput showers_completed(start=0) annotation(
    Placement(visible = true, transformation(origin = {80, -14}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {100, -92}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));

algorithm

//#################### RESPONSE TO START SHOWER COMMAND
  chattering_counter_value := shower_command;
  when chattering_counter >= 0.5 * t_p then
    if showers_remaining > 0 then
      if flow_water >= desired_water_flow or true then
        auxiliar := true;
        if auxiliar and not is_shower_active_1 then
          auxiliar := false;
          shower_time_counter_value_1 := 1;
          is_shower_active_1 := true;
        end if;
        if auxiliar and not is_shower_active_2 then
          auxiliar := false;
          shower_time_counter_value_2 := 1;
          is_shower_active_2 := true;
        end if;
        if auxiliar and not is_shower_active_3 then
          auxiliar := false;
          shower_time_counter_value_3 := 1;
          is_shower_active_3 := true;
        end if;
        if auxiliar and not is_shower_active_4 then
          auxiliar := false;
          shower_time_counter_value_4 := 1;
          is_shower_active_4 := true;
        end if;
        if auxiliar and not is_shower_active_5 then
          auxiliar := false;
          shower_time_counter_value_5 := 1;
          is_shower_active_5 := true;
        end if;
      end if;
    end if;
  end when;
//#################### LOOPING TO SET PARAMETERS
  showers_remaining := 0;
// Reset value of showers remaining to 0 and then do the count later in FOR loop
// Stop a shower and count it
  when shower_time_counter_1 >= time_per_shower then
    is_shower_active_1 := false;
    count := count + 1;
    shower_time_counter_value_1 := 0;
  end when;  
    when shower_time_counter_2 >= time_per_shower then
      is_shower_active_2 := false;
      count := count + 1;
      shower_time_counter_value_2 := 0;
    end when;  
    when shower_time_counter_3 >= time_per_shower then
      is_shower_active_3 := false;
      count := count + 1;
      shower_time_counter_value_3 := 0;
    end when;  
    when shower_time_counter_4 >= time_per_shower then
      is_shower_active_4 := false;
      count := count + 1;
      shower_time_counter_value_4 := 0;
    end when;    
    when shower_time_counter_5 >= time_per_shower then
      is_shower_active_5 := false;
      count := count + 1;
      shower_time_counter_value_5 := 0;
    end when;
// Count not active showers
  if not is_shower_active_1 then
    showers_remaining := showers_remaining + 1;
  end if;
    if not is_shower_active_2 then
      showers_remaining := showers_remaining + 1;    
    end if;
    if not is_shower_active_3 then
      showers_remaining := showers_remaining + 1;    
    end if;    
    if not is_shower_active_4 then
      showers_remaining := showers_remaining + 1;    
    end if;
    if not is_shower_active_5 then
      showers_remaining := showers_remaining + 1;    
    end if;

equation
//#################### EQUATIONS
    der(shower_time_counter_1) = shower_time_counter_value_1; 
    der(shower_time_counter_2) = shower_time_counter_value_2;     
    der(shower_time_counter_3) = shower_time_counter_value_3;   
    der(shower_time_counter_4) = shower_time_counter_value_4;   
    der(shower_time_counter_5) = shower_time_counter_value_5; 
    der(chattering_counter)=chattering_counter_value; 
    
  when chattering_counter > t_p then
    reinit(chattering_counter, 0);
  end when;
  when chattering_counter_value == 0 then
    reinit(chattering_counter, 0);
  end when;


    when not is_shower_active_1 then
      reinit(shower_time_counter_1, 0);
    end when;   
    when not is_shower_active_2 then
      reinit(shower_time_counter_2, 0);
    end when;  
    when not is_shower_active_3 then
      reinit(shower_time_counter_3, 0);
    end when;  
    when not is_shower_active_4 then
      reinit(shower_time_counter_4, 0);
    end when;  
    when not is_shower_active_5 then
      reinit(shower_time_counter_5, 0);
    end when; 


algorithm
//#################### OUTPUT ALGORITHM

showers_completed := count;
desired_water_flow := (n_shower-showers_remaining)*water_flow_per_shower;

  annotation(
    uses(Modelica(version = "4.0.0")),
    Icon(graphics = {Rectangle(origin = {0, 1}, fillColor = {0, 255, 127}, fillPattern = FillPattern.Solid, extent = {{-94, 93}, {94, -93}}), Text(origin = {49, 1}, extent = {{-23, 19}, {23, -19}}, textString = "Water needed"), Text(origin = {-43, 1}, extent = {{-21, 15}, {21, -15}}, textString = "Mass flow in"), Text(origin = {0, 55}, extent = {{-26, 7}, {26, -7}}, textString = "Request command"), Text(origin = {1, -70}, extent = {{-35, 14}, {35, -14}}, textString = "Showers available"), Text(origin = {89, -77}, extent = {{-37, 7}, {37, -7}}, textString = "Showers completed")}));
end shower;

clear all;
timestep = 1/20;
global graphstep; %Time step in hours
graphstep = 1/20;
%% Material Properties
K_rw = 0.044; %rock_wool_conductivity %w/m.
K
K_pes = 0.4; %polyethelene_sheet_conductivity %w/m.
K
K_pss = 45.2736; %prefabricated_steel_sheet_conductivity %w/m.
K
K_css = 45.2736; %corrugated_steel_sheet_conductivity
assuming same as prefabricated %w/m.K
alpha_steel = 0.354; %absorptance_steel
h_in = 8.29; %internal_air_h %
W/m^2.K
K_lsg = 1; %Laminated_safety_glass_conductivity %
W/m.K did not find the data so assuming the properties of glass
K_leg = 1; %Low_emissivity_glass_conductivity %W/m.
K did not find the data so assuming the properties of glass
beta_lsg = 0.35; %Low_emissivity_glass_transmittance %
solar transmittance
h_af = 6.25; %air_facade_h % W/m^2.K
K_conc = 1.4; %concrete_conductivity % W/m.K
Cp_air = 1005; %cp of air % J/Kg.K
rho_air = 1.225; %air density % kg/m^3
%% Computation
T_room = zeros(24/timestep,1);
T_comfort = zeros(1,288/timestep); % 12*24 = 288 taking 24 hours for each of the
month
T_room = zeros(1,288/timestep);
T_external = zeros(1,288/timestep);
T_room_original = zeros(1,288/timestep);
HVAC = zeros(1,288/timestep);
Solar_Irradiance = zeros(1,288/timestep);
delta_t = timestep*3600; % seconds
t_12 = [timestep:timestep:288];
counter = 0;
for j = 1:12
 if j==1 || j==2 || j==3 || j==4 || j==11 || j==12
 h_ext = 34; %external_air_h_winter %
W/m^2.K
 f_f = 0.3; % reduction factor due to moveable curtains
 elseif j==5 || j==6 || j==7 || j==8 || j==9 || j==10
 h_ext = 22.7; %external_air_h_summer %
W/m^2.K
 f_f = 0.8; % reduction factor due to moveable curtains
 else
 fprintf("Invalid Month");
 end
 %% Areas and Volume of Different Places/Walls
 Areas_L1_L2 = 129 ; % m^2
 Area_T = 1092.55; % m^2
 Area_B1_B2_without_door = 53.6; % m^2
 Area_B1_B2_door = 29.6; % m^2
 Volume_total = 3828.93; % m^3
 %% Resistance of sides of the Building named as L1 and L2
 t_g1 = 0.020; % m
 t_ag1= 1; % m
 t_g2 = 0.008; % m
 t_ag2 = 0.012; % m
 t_g3 = 0.012; % m
 R_L=((1/h_ext)+t_g1/(K_lsg)+(1/h_af)+t_g2/(K_leg)+(1/h_af)+t_g3/(K_leg)+(1/h_in));
 %% Resistance of the top
 external_steel_thickness = 0.001; %m
 polythene_sheet_thickness = 0.001; %m
 rock_wool_thickness = 0.080; %m
 internal_steel_thickness = 0.0006; %mm
 R_T= (1/h_ext + external_steel_thickness/K_pss + polythene_sheet_thickness/(K_pes) +
rock_wool_thickness/(K_rw) + internal_steel_thickness/(K_pss) + (1/h_in));
 %% Resistance of the front and Back without door
 R_fb = (1/h_ext + external_steel_thickness/K_pss + polythene_sheet_thickness/(K_pes)
+ rock_wool_thickness/(K_rw) + internal_steel_thickness/(K_pss) + (1/h_in));
 %% Resistance of B1 Door
 t_d_g1=0.008; % m
 t_d_g2=0.012; % m
 R_D = ((1/h_ext) + t_d_g1/(K_lsg) + 1/h_af + t_d_g2/K_lsg + 1/h_in);
 %% UA of the building
 UA = 2*Areas_L1_L2/R_L + Area_T/R_T + 2*Area_B1_B2_without_door/R_fb +
2*Area_B1_B2_door/R_D ;
 %% Ventilation/Convection Heat Transder
 n = 0.7/3600; % number of air shifts per second
 %% Heat Generation Inside the Building
 number_of_people = 25;
 number_of_PCs = 25;
 heat_per_person = 130; % watts Reference of Cenegel book
 heat_per_PC = 50; % watts Reference of Cenegel book
 total_heat_generation = (number_of_people*heat_per_person) +
(number_of_PCs*heat_per_PC);
 %% Solar Irradiation Heat Addition
 A_S_L1_L2 = beta_lsg*Areas_L1_L2*(1-f_f); % Asun for a sides L1 and L2 same for the
both
 A_S_B1_B2_door = beta_lsg*Area_B1_B2_door; % Asun for doors in the sides B1 and B2
 A_S_B1_B2 = (alpha_steel*Area_B1_B2_without_door)/ R_fb; % Asun for B1 and B2
without doors
 A_S_top = (alpha_steel*Area_T)/R_T;
 h_rad = 6;
 %% Temperature Profile
 A = readmatrix('SOLAR DATA SHAFWAN\Roof_Yearly.xlsx');
 t = A(1:24,25);
 T_ext_int=zeros(24/timestep,1);

 T_ext = A(1:24,2*j);
 t_int=[0:timestep:24];
 T_ext_int=interp1(t,T_ext,t_int,"spline");
 %% Solar Irradiation Data
 G_t = readmatrix('SOLAR DATA SHAFWAN\Roof_Yearly.xlsx');
 G_top = G_t(1:24,(2*j)-1);
 G_top_int = interp1(t,G_top,t_int, "spline");
 G_Le1 = readmatrix('SOLAR DATA SHAFWAN\Wall_L1_Yearly.xlsx');
 G_L1 = G_Le1(1:24,j);
 G_L1_int = interp1(t,G_L1,t_int,"spline");
 G_Le2 = readmatrix('SOLAR DATA SHAFWAN\Wall_L2_Yearly.xlsx');
 G_L2 = G_Le2(1:24,j);
 G_L2_int = interp1(t,G_L2,t_int,"spline");
 G_Be1 = readmatrix('SOLAR DATA SHAFWAN\Wall_B1_Yearly.xlsx');
 G_B1 = G_t(1:24,j);
 G_B1_int = interp1(t,G_B1,t_int,"spline");
 G_Be1 = readmatrix('SOLAR DATA SHAFWAN\Wall_B2_Yearly.xlsx');
 G_B2 = G_t(1:24,j);
 G_B2_int = interp1(t,G_B2,t_int,"spline");

 T = T_ext(1);
 T_orig = T_ext(1);
 term1 = Cp_air*rho_air*Volume_total;
 term = (delta_t/(Cp_air*rho_air*Volume_total))*(UA);
 T_comfort_local = zeros(24/timestep,1);
 T_room_local = zeros(24/timestep,1);
 T_room_original_local = zeros(24/timestep,1);
 T_room_original_local = zeros(24/timestep,1);
 HVAC_local = zeros((24/timestep)+1,1);
 HVAC_local(1) = 0;
 for i = 1:(24/timestep)
 counter = counter + 1;
 if i<(18/timestep) && i>(8/timestep) % This condition counts the heat
generation only during the working time
 y = 1;
 else
 y = 0;
 end
 T_comfort_local(i) = 23.9 + 0.295*(T_ext_int(i)-22)*exp(-((T_ext_int(i)-22)/
(22*1.414))^2);
 T_comfort(counter) = T_comfort_local(i);
 T_external(counter) = T_ext_int(i);
 T_room_local(i) = T;
 T_room_original_local(i) = T_orig;
 T_orig = T_orig + term*(T_ext_int(i)-T_orig) + y*n*(T_ext_int(i)-T_orig)*delta_t+
y*(total_heat_generation)*delta_t/term1 + (A_S_top*G_top_int(i))/((h_rad+h_ext)*term1)
*delta_t + (A_S_L1_L2*G_L1_int(i)+A_S_L1_L2*G_L2_int(i)+A_S_B1_B2*G_B1_int(i)/
(h_rad+h_ext)+(A_S_B1_B2*G_B2_int(i))/(h_rad+h_ext)+ A_S_B1_B2_door*G_B1_int(i)
+A_S_B1_B2_door*G_B2_int(i) )*delta_t/term1;
 T = T+ ((HVAC_local(i)*delta_t)/term1)+ term*(T_ext_int(i)-T) + y*n*(T_ext_int(i)
-T)*delta_t+ y*(total_heat_generation)*delta_t/term1 + (A_S_top*G_top_int(i))/
((h_rad+h_ext)*term1)*delta_t + (A_S_L1_L2*G_L1_int(i)+A_S_L1_L2*G_L2_int(i)
+A_S_B1_B2*G_B1_int(i)/(h_rad+h_ext)+(A_S_B1_B2*G_B2_int(i))/(h_rad+h_ext)+
A_S_B1_B2_door*G_B1_int(i)+A_S_B1_B2_door*G_B2_int(i) )*delta_t/term1;
 T_room(counter) = T_room_local(i);
 T_room_original(counter) = T_room_original_local(i);
 if j == 1 || j == 2 || j==3 || j==4 || j== 11 || j == 12
 if i<(18/timestep+1) && i>(8/timestep+1)
 if T < T_comfort_local(i)
 HVAC_local(i+1) = HVAC_local(i) + 250*(T_comfort_local(i)-T);
 elseif T>T_comfort_local(i)
 HVAC_local(i+1) = HVAC_local(i) - 250*(T-T_comfort_local(i));
 else
 HVAC_local(i+1) = HVAC_local(i);
 end
 HVAC(counter) = HVAC_local(i);
 end
 else
 if i<(18/timestep+1) && i>(8/timestep+1)
 if T < T_comfort_local(i)
 HVAC_local(i+1) = HVAC_local(i) + 400*(T_comfort_local(i)-T);
 elseif T>T_comfort_local(i)
 HVAC_local(i+1) = HVAC_local(i) - 400*(T-T_comfort_local(i));
 else
 HVAC_local(i+1) = HVAC_local(i);
 end
 HVAC(counter) = HVAC_local(i);
 end
 end
 Solar_Irradiance(counter) = G_top_int(i);
 end
end
HVAC=HVAC*1.1;
%% Plotting of Graphs
figure(1);
plot(t_12,T_comfort)
hold on
plot(t_12,T_room)
hold on
plot(t_12,T_external)
hold on
plot(t_12,T_room_original)
legend('Comfort Temperature','Troom with HVAC','Ambinet Temperature','Troom without
HVAC')
xlabel("Time (Hours)")
ylabel("Temperature (C)")
figure(2);
plot(t_12,HVAC)
xlabel("Time")
ylabel("HVAC Load")
csvwrite('HVAC.csv', HVAC)
%% Energy Required from Storage Tank
COP_chiller = 0.7;
COP_radiator = 0.9;
Load = zeros(1,288/timestep);
for a = 1:(288/timestep)
 if HVAC(a) < 0
 Load(a) = -(HVAC(a)/COP_chiller);
 else
 Load(a) = (HVAC(a)/COP_radiator);
 end
end
figure(3);
plot(t_12,Load)
xlabel("Time");
ylabel("Energy Load from Storage Tank");
%% Energy Produced by One Solar Collector on Each Average Day of 12 Months
tiledlayout(2,2)
E=zeros(1,12);
E_load = cumsum(Load)*timestep;
E_panel = zeros(1,5760);
E_panel = cumsum(2.33*0.5*Solar_Irradiance)*timestep;
x=(length(t_12)/12);
y=x;
for i=1:12
 E(i)=E_panel(y);
 y=y+x;
end
nexttile
bar(E,'y')
title('commulative energy over 12 months')
for i=12:-1:2
 E(i)=E(i)-E(i-1);
end
nexttile
bar(E,'r')
title('Monthly Energy Production')
month = { 'Jan' , 'Feb' , 'Mar' , 'Apr' , 'May' , 'Jun' , 'Jul' , 'Aug' , 'Sep' , 'Oct' ,
'Nov' , 'Dec'};
xticks(1:1:numel(month)) % give an increment for ticks
xticklabels(month(1:1:end))
writematrix(E,'Monthy Energy Production.csv')
E=cumsum(E);
nexttile
plot(E,'g')
title('summation check')
nexttile
bar(E(12))
title('1st year energy production')
%% Optmization Starts Here
global collector_efficiency
collector_efficiency = 0.5;
global density Cp
density = 1000;
global area
area = 2.33;
Cp = 4186;
% Time definition
% A day will be defined by 24 hours by steps of 30 minutes :
global deltaTT deltaT Time
deltaTT = timestep; % half an hour
Time = deltaTT:deltaTT:288;
deltaT = deltaTT*3600;
%% PV power
global GSR % global solar radiation
global P_load % load required from solar collectors
GSR = Solar_Irradiance; % W/m^2
T_PV = T_external; % celsius degrees
P_load = Load;
%%
%
%
% Battery and system specifications including investment costs
global eta_tank eta_sc_tank eta_tank_ac Cinv_bat Cinv_PV Cinv_tank Cinv_collector
eta_tank = 0.92; % charging efficiency
eta_sc_tank = 0.92; % DC/DC converter efficiency
eta_tank_ac = 0.92; % AC DC converter efficiency
Cinv_tank = 473; % price in euros for 1m^3 tank volume
Cinv_collector = 701; % Euros/solar collector
E_load = cumsum(P_load)*deltaT; % one day energy of the load
E_collector = cumsum(comp_P_SC(1))*deltaT; % one day energy of one single panel
Nsc_min = 4; % minimal number of panels
% bounds for SOC constraints
global temp_min temp_max
temp_min = 343; % min temperature
temp_max = 353; % max temperature
Nsc_max = 500; % max collectors
Ntank_min = 1; % min volume
Ntank_max = 10; % max volume
temp_ini_min = 343; % min initial temperature
temp_ini_max = 353; % max initial temperature
% lower and upper bounds constraints
Xl = [Nsc_min;Ntank_min;temp_ini_min];
Xu = [Nsc_max;Ntank_max;temp_ini_max];
obj = @(X) X(1)*Cinv_collector + X(2)*Cinv_tank; % investment cost
% no linear constraints
% non-linear constraints are defined by a specific function : "constraints"
options = optimoptions(@ga);
[Xopt,cost] = ga(obj,3,[],[],[],[],Xl,Xu,@constraints,[1 2],options)
[P_pv,P_bat,SOC,P_lost,P_unsup] = comp_Power(Xopt);
figure(4)
plot(Time,P_load,Time,P_pv,Time,P_bat,Time,P_lost,Time,P_unsup)
legend('Load','SC','Tank','Lost','Unsupplied')
%% Functions definition
% SC power
function P_SC = comp_P_SC(Nsc)
 % compute the power produced by Npv panel during a day
 global GSR collector_efficiency area
 P_SC = Nsc*area*collector_efficiency*GSR;
end
%%
% Powers computation
function [P_sc,P_tank,temp,P_lost,P_unsup] = comp_Power(X)
 % compute the powers in the different elements
 global deltaT P_load eta_sc_tank eta_tank_ac eta_tank temp_min temp_max density Cp
 Nsc = X(1);
 Ntank = X(2);
 temp_ini = X(3);

 P_sc = comp_P_SC(Nsc);
 n = length(P_sc);
 P_unsup = zeros(1,n);
 P_tank = zeros(1,n);
 P_lost = zeros(1,n);
 temp = zeros(1,n+1);

 temp(1) = temp_ini;
 for k = 1:n
 if P_sc(k) >= P_load(k)
 if (temp(k) < temp_max)
 P_tank(k) = -eta_sc_tank*(P_sc(k)-P_load(k));
 temp(k+1) = min(temp(k)-eta_tank*P_tank(k)*deltaT/(Ntank*Cp*density),
temp_max);
 else
 P_lost(k) = P_sc(k)-P_load(k);
 temp(k+1) = temp(k);
 end
 else
 if (temp(k) > temp_min)
 P_tank(k) = -eta_tank*(P_sc(k)-P_load(k));
 temp(k+1) = max(temp(k)-eta_tank_ac*P_tank(k)*deltaT/(Ntank*density*Cp),
temp_min);
 else
 P_unsup(k) = -(P_sc(k)-P_load(k));
 temp(k+1) = temp(k);
 end
 end
 end
end
% Optimization problem constraints
function [c,ceq] = constraints(X)
 % constraints of our optimization problem
 global temp_min temp_max

 [~,~,temp,~,P_unsup] = comp_Power(X);


 c1 = temp_min - min(temp);

 c2 = max(temp) - temp_max;
 c3 = max(P_unsup);

 c4 = temp(1)-temp(length(temp)) ;

 c = [c1 c2 c3 c4];
 ceq = [];
end
%% Optimization of PV Plant (System A)
% A day will be defined by 24 hours by steps of 30 minutes
timestep= 1/20;
global deltaT Time
deltaT = 1/20; % half an hour
Time = deltaT:deltaT:288;
timestep=1/20;
t_12 = timestep:timestep:288;
%% For Calculation of Heat Pump Load
% Making the whole load positive and converting heating/cooling load to electrical load
global GSR T_PV
Solar_Irradiance= readmatrix('Solar_Irradiance.csv');
T_external= readmatrix('T_external.csv');
GSR=Solar_Irradiance;
T_PV=T_external;
global P_load
B = readmatrix('HVAC.csv');
figure(13)
plot(B)
Q_summer = 0.27; %taking into account DC-DC, DC-AC converter efficiencies and
heat pump COP for summer
Q_winter = 0.23; %taking into account DC-DC, DC-AC converter efficiencies and
heat pump EER for winter
for i = 1:(length(B))
 if B(i)>0
 B(i)=B(i)*Q_summer;
 else
 B(i)=B(i)*(-1)*Q_winter;
 end
end
tiledlayout(3,3)
nexttile
P_load=B;
plot(P_load)
%% Solar Irradiance for PV Panels
nexttile
plot(t_12,Solar_Irradiance)
xlabel("Time")
ylabel("Solar_Irradiance")
%% PV power
global P_PV_STC % PV array rated power in standard test conditions
25/06/22 16:26 plost_elim.m 2 of 5
(STC)
global G_STC % solar radiation under STC
global CT % temperature coefficient
global T_STC % temperature STC
global eta_PV % panels efficiency (electrical)
 % temperature
P_PV_STC = 410; % Capacity of each PV Panel %Watts
G_STC = 1000; % Irradiance at standard test conditions W/m^2
CT = -0.36e-2; % efficiency change in panel per degree temperature
difference
T_STC = 25; % Standard temperature in celsius
eta_PV = 0.95; % electrical efficiency of PV
% Plot of power produced by one panel in one day
nexttile
plot(Time,comp_P_PV(1),'LineWidth',2)
axis([0 24 0 inf],'auto y')
xticks(0:2:24)
xlabel('Day hours')
ylabel('One panel power (W)')
grid
%% Battery and system specifications including investment costs
global C_cell eta_bat eta_DCDC eta_ACDC Cinv_bat Cinv_PV Cinv_cell Cinv_panel
C_cell = 230*12; % cell capacity (W.h) 230 Ah / 12 V
eta_bat = 0.95; % charging efficiency
eta_DCDC = 0.98; % DC/DC converter efficiency
eta_ACDC = 0.95; % AC DC converter efficiency
Cinv_bat = 330; % € / kWh
Cinv_PV = 900; % € / kW
Cinv_cell = Cinv_bat*1e-3*C_cell;
Cinv_panel = Cinv_PV*1e-3*P_PV_STC;
E=rand(12,1);
E_load = cumsum(P_load)*deltaT;% one day energy of the load
E_panel = cumsum(comp_P_PV(1))*deltaT;
x=(length(t_12)/12);
y=x;
for i=1:12
 E(i)=E_panel(y);
 y=y+x;
end
nexttile
25/06/22 16:26 plost_elim.m 3 of 5
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
% one day energy of one single pane
plot(E_panel)
Npv_min = nearest(E_load/E_panel); % minimal number of panels
%% Bounds on constraints
% Bounds for SOC constraints
global SOC_min SOC_max
SOC_min = 0.05;
SOC_max = 0.95;
Npv_max = 500;
Ncell_min = 25;
Ncell_max = 5000;
SOC_ini_min = SOC_min*1e2;
SOC_ini_max = SOC_max*1e2;
% Lower and upper bound constraints
Xl = [Npv_min;Ncell_min;SOC_ini_min];
Xu = [Npv_max;Ncell_max;SOC_ini_max];
%% Objective Function
obj = @(X) X(1)*Cinv_panel+X(2)*Cinv_cell; % investment cost
options = optimoptions(@ga);
[Xopt,cost] = ga(obj,3,[],[],[],[],Xl,Xu,@constraints,[1 2],options)
[P_pv,P_bat,SOC,P_lost,P_unsup] = comp_Power(Xopt);
25/06/22 16:26 plost_elim.m 4 of 5
nexttile
plot(Time,P_load,Time,P_pv,Time,P_bat,Time,P_lost,Time,P_unsup)
legend('Load','PV','Battery','Lost','Unsupplied','Location',"bestoutside")
%% Functions definition
% PV power
function P_PV = comp_P_PV(Npv)
 % compute the power produced by Npv panel during a day
 global P_PV_STC G_STC CT T_STC eta_PV GSR T_PV
 P_PV = Npv*eta_PV*P_PV_STC/G_STC*(1-CT*(T_PV-T_STC)).*GSR;
end
% Powers computation
function [P_pv,P_bat,SOC,P_lost,P_unsup] = comp_Power(X)
 % compute the powers in the different elements
 global deltaT P_load eta_DCDC eta_ACDC eta_bat C_cell SOC_min SOC_max
 Npv = X(1);
 Ncell = X(2);
 SOC_ini = X(3)*1e-2;
 P_pv = comp_P_PV(Npv);
 n = length(P_pv);
 P_unsup = zeros(1,n);
 P_bat = zeros(1,n);
 P_lost = zeros(1,n);
 SOC = zeros(1,n+1);

 SOC(1) = SOC_ini;
 for k = 1:n
 if ( (eta_DCDC*P_pv(k)) >= (P_load(k))/eta_ACDC)
 if (SOC(k) < SOC_max)
 P_bat(k) = -eta_DCDC*(eta_DCDC*P_pv(k)-P_load(k)/eta_ACDC);
 SOC(k+1) = min(SOC(k)-eta_bat*P_bat(k)*deltaT/(Ncell*C_cell),SOC_max);
 else
 P_lost(k) = eta_DCDC*P_pv(k)-P_load(k)/eta_ACDC;
 SOC(k+1) = SOC(k);
 end
 else
 if (SOC(k) > SOC_min)
 P_bat(k) = -1/eta_DCDC*(eta_DCDC*P_pv(k)-P_load(k)/eta_ACDC);
 SOC(k+1) = max(SOC(k)-P_bat(k)*deltaT/(Ncell*C_cell),SOC_min);
 else
 P_unsup(k) = -(eta_DCDC*P_pv(k)-P_load(k)/eta_ACDC);
 SOC(k+1) = SOC(k);
 end
 end
 end
25/06/22 16:26 plost_elim.m 5 of 5
end
function [c,ceq] = constraints(X)
 % constraints of our optimization problem
 global SOC_min SOC_max

 [~,~,SOC,~,P_unsup] = comp_Power(X);


 c1 = SOC_min - min(SOC);

 c2 = max(SOC) - SOC_max;
 c3 = max(P_unsup);

 c4 = SOC(1)-SOC(length(SOC)) ;

 c = [c1 c2 c3 c4];
 ceq = [];
end
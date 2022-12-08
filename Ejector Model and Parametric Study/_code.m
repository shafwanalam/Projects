%% *Code by Md Shafwan Ul Alam and Lisa Zumsteeg* 

clear all;
close all;
clc;


pyversion
[v,e] = pyversion; system([e,' -m pip uninstall -y CoolProp']);
[v,e] = pyversion; system([e,' -m pip install --user -U CoolProp']);

import py.CoolProp.CoolProp.*

%% Given Values

fluid = 'R141b';
m = PropsSI ('M', fluid ); % molar mass [kg/mol]


% Geometry data and efficiencies
Dt = 2.64*10^-3; % diameter at throat(t) position
D1 = 4.5* 10^-3; % diameter at point 1
D2 = 6.7* 10^-3; % diameter at point 2

eta_p = 0.95;
eta_py = 0.88;
eta_s = 0.85;
Phi_m = 0.82; %for frictional losses

%Pressure [Pa] and temperature [K] of primary(p) and secondary(s) fluids
K = 273.15; 
Tp = 95 + K; 
Ts = 8+K;
Pp = 6.04 * 10^5;
Ps = 0.4 * 10^5;

%% Calculation of the primary flow, mp

Mt = 1; % Mach number at position t, primary flow supposed to be choked
R = 8.314;
Rg = R/m; %[J/(kg.K)]
Cp = PropsSI ('C', 'T', Tp, 'P', Pp, fluid ); %[J/(kg.K)]
Cv = Cp - Rg; %[J/(kg.K)]
gamma = (Cp / Cv);
At = pi * (Dt/2).^2; %Area at position t given by Dt

mp = sqrt(gamma/Rg)*(Pp/sqrt(Tp)) * Mt * (1+ ((gamma-1)/2)*Mt^2).^(-(gamma+1)/(2*(gamma-1)))* At* sqrt(eta_p);

%% Determination of conditions for primary fluid at point 1

% Mach number, Mp1
Ap1 = pi * (D1/2).^2; %Area at point 1 for primary fluid given by D1
fun1 = @(Mp1) (-Ap1/At + ((Mt/Mp1) * ((1+(gamma-1)/2 * Mt.^2)/(1+(gamma-1)/2 * Mp1.^2)).^(-(gamma+1)/(2*(gamma-1)))));
Mp1 = fzero(fun1,1.5); %resolution of Mp1

% Pressure, Pp1
Pp1 = Pp * (1+ ((gamma-1)/2) * Mp1.^2).^(-gamma/(gamma-1));

%% Determination of conditions at position y (start of mixing fluids)

% M and P conditions for secondary fluid
Msy = 1; % Mach number at position y, secondary flow supposed to be choked
Psy = Ps * (1+ ((gamma-1)/2)*Msy.^2).^(-gamma/(gamma-1));

% M and P conditions for primary fluid (deduced from secondary fluid)
Ppy = Psy; % same pressure for p than s at the start of mixing fluids
fun2 = @(Mpy) (-Ppy/Pp1 + ((1+(gamma-1)/2 * Mp1.^2)/(1+(gamma-1)/2 * Mpy.^2)).^(-(gamma+1)/(2*(gamma-1))));
Mpy = fzero(fun2,1.5); %resolution of Mpy

% Section conditions for primary and secondary fluid
Apy = Ap1 * eta_py * (Mp1/Mpy) * (((1+(gamma-1)/2 * Mp1.^2)/(1+(gamma-1)/2 * Mpy.^2)).^(-(gamma+1)/(2*(gamma-1))));
A2 = pi * (D2/2).^2; % where A2=Ay because positions are in the constant area
Asy = A2 - Apy;

%% Calculation of the secondary flow, ms
ms = sqrt(gamma/Rg)* Ps/ sqrt(Ts) * Msy * (1+(gamma-1)/2 * Msy.^2).^(-(gamma+1)/(2*(gamma-1))) * Asy * sqrt(eta_s);

%% Determination of conditions for the mixing section at m position

% Velocity condition, Vm by deducing the velocities Vpy and Vsy at y position
Tpy = Tp * (1+ (gamma-1)*Mpy.^2/2).^-1;
Tsy = Tp * (1+ (gamma-1)*Msy.^2/2).^-1;
Vpy = Mpy * sqrt(gamma* Rg * Tpy);
Vsy = Msy * sqrt(gamma * Rg * Tsy);

Vm = Phi_m*(mp*Vpy+ms*Vsy)/(mp+ms); %knowing the mass flow of the 2 fluids

% Mach number condition, Mm by deducing Tm
C1 = mp*(Cp*Tpy+ Vpy.^2/2) + ms*(Cp*Tsy+ Vsy.^2/2);
fun3 = @(Tm) (-C1 + (mp+ms)*(Cp*Tm+ Vm.^2/2));
Tm = fzero(fun3,2);%Resolution of Tm

Mm = Vm/sqrt(gamma* Rg* Tm);

% Pressure condition, Pm
Pm = Ppy; % same pressure in the mixing constant area

%% Determination of conditions at position 2 (separation of constant area and diffuser)
% Pressure condition, P2
P2 = Pm * (1+(2*gamma/(gamma+1))*(Mm.^2-1));

% Mach number, M2
M2 = sqrt((1+(gamma-1)*Mm.^2/2)/(gamma*Mm.^2 - (gamma-1)/2));

%% Determination of the condition at position c (stagnation)
% Pressure Pc
Pc = P2 / (1+(gamma-1)*M2.^2/2).^(-gamma/(gamma-1));

%% Performance parameter of the ejector
Ent_ratio = ms/mp; % Entrainment Ratio


%% 2/ Representation of the Pressure and Velocity conditions progress

%Hypothetical Distance along the ejector
dis_prim_fluid = [0 1 1.5 2 4.5 6 6.01 8];
dis_sec_fluid  = [0.5 2 4.5 6 6.01 8];

%Values of the primary fluid
PrimaryPressure = [Pp 0.5*(-Pp1+Pp) Pp1 Ppy Pm Pm P2 Pc];
PrimaryVelocity = [0 Mt Mp1 Mpy Mm Mm M2 0];

%Values of the secondary fluid
SecondaryPressure = [Ps Psy Pm Pm P2 Pc];
SecondaryVelocity = [0 Msy Mm Mm M2 0];

%% Plot
%{we can see that the velocity of the primary fluid increases rapidly once
%it enters the ejector, for the secondary fluid the velocity also increases
%but not as rapidly as the primary fluid. For the primary fluid the
%velocity is at its peak just before entering the mixing chamber, and once
%the mixing starts the velocity of the primary fluid keeps on decreasing
%and the secondary fluid keeps on increasing, until they are at an
%equilibrium about half way down the mixing chamber. The moment the mixed fluid
%leaves the mixing chamber and enters the diffuser, there is a sudden drop
%of velocity, and after that the valocity keeps on decreasing slowly as it
%travels through the diffuser.}

figure(1)
plot(dis_sec_fluid, SecondaryVelocity)
hold on
plot(dis_prim_fluid, PrimaryVelocity)
xlabel('Hypothetical Distance along the ejector');
ylabel('Velocity variation');
title('Velocity Variation in the Ejector')
legend('Secondary fluid','Primary fluid')

%{we can see that the pressure of the primary fluid decreases rapidly once
%it enters the ejector, for the secondary fluid the pressure also decreases
%but not as rapidly as the primary fluid. The pressure of both the fluids 
% reaches eequilibrium just before entering the mixing chamber, and remains
% constant throughout the mixing chamber. The moment the mixed fluid
%leaves the mixing chamber and enters the diffuser, there is a sudden
%increase of pressure, and after that the pressure keeps on increasing slowly
% as it travels through the diffuser.}

figure(2)
plot(dis_sec_fluid, SecondaryPressure)
hold on
plot(dis_prim_fluid, PrimaryPressure)
xlabel('Hypothetical Distance along the ejector');
ylabel('Pressure variation');
title('Pressure Variation in the Ejector')
legend('Secondary fluid','Primary fluid')
disp(Ent_ratio)

%{
Description:
For the nominal point, we can see the progress of conditions inside the ejector gas-gas.
The variation velocity is represented with the Mach number, indeed when
M=1, we know the fluid is at the sonic velocity.
The fluid draws in a secondary fluid and produces a mixed fluid at the ejector outlet. 
The primary fluid has a higher total pressure than the ejector outlet, and higher than the secondary fluid.
So, the primary fluid undergoes acceleration at the throat and reaches its sonic velocity.
The p fluid becomes supersonic (M>1) until the constant area.
After the primary fluid reached the throat, the pressure and the velocity start to decrease. 
Once the primary and secondary fluids are mixed, they enter a diffuser where the velocity is diffused that means converted into
pressure. It's why the velocity decrease and the pressure increase in the diffuser.
%} 
%% 1. For Changing Dt from 0.00264 to 0.002 and 0.003 and Ent_ratio changes from 0.3571 to 0.4973 and 0.2873
% 
% <<C:\Users\shafw\Downloads\Graph\1st Dt 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\1st Dt 2.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd Dt 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd Dt 2.jpg>>
% 
%% Parametric Study Result Discussion 1
% 1. For Changing the values of Dt we notice that by increasing the value
% we decrease the Value of Entrainment ratio and vice versa.
%% 2. For Changing D1 from 0.0045 to 0.004 and 0.005 and Ent_ratio changes from 0.3571 to 0.3896 and 0.3069
% 
% <<C:\Users\shafw\Downloads\Graph\1st D1 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\1st D1 2.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd D1 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\1st D1 2.jpg>>
% 
%% Parametric Study Result Discussion 2
% 2. For Changing the values of D1 we notice that by increasing the value
% we decrease the Value of Entrainment ratio and vice versa.
%% 3. For Changing D2 from 0.0067 to 0.006 and 0.007 and Ent_ratio changes from 0.3571 to 0.2656 and 0.3993
%%
% 
% <<C:\Users\shafw\Downloads\Graph\1st D2 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\1st D2 2.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd D2 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd D2 2.jpg>>
% 
%% Parametric Study Result Discussion 3
% 3. For Changing the values of D2 we notice that by increasing the value
% we increase the Value of Entrainment ratio and vice versa.
%% 4. For Changing Pp from 6.04 to 5 and 7 and Ent_ratio changes from 0.3571 to 0.4133 and 0.318
%%
% 
% <<C:\Users\shafw\Downloads\Graph\1st Pp 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\1st Pp 2.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd Pp 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd Pp 2.jpg>>
% 
%% Parametric Study Result Discussion 4
% 4. For Changing the values of P_primary we notice that by increasing the value
% we decrease the Value of Entrainment ratio and vice versa.
%% 5. For Changing Ps from 0.4 to 0.3 and 0.5 and Ent_ratio changes from 0.3571 to 0.281 and 0.4246
%%
% 
% <<C:\Users\shafw\Downloads\Graph\1st Ps 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\1st Ps 2.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd Ps 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd Ps 2.jpg>>
% 
%% Parametric Study Result Discussion 5
% 5. For Changing the values of P_secondary we notice that by increasing the value
% we increase the Value of Entrainment ratio and vice versa.
%% 6. For Changing Tp from 95 to 90 and 100 and Ent_ratio changes from 0.3571 to 0.3566 and 0.3594
%%
% 
% <<C:\Users\shafw\Downloads\Graph\1st Tp 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\1st Tp 2.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd Tp 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd Tp 2.jpg>>
% 

%% Parametric Study Result Discussion 6
% 6. For Changing the values of T_primary we notice that by increasing the value
% we increase the Value of Entrainment ratio and vice versa.
%% 7. For Changing Ts from 8 to 0 and 15 and Ent_ratio changes from 0.3571 to 0.3623 and 0.3527
%%
% 
% <<C:\Users\shafw\Downloads\Graph\1st Ts 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\1st Ts 2.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd Ts 1.jpg>>
% 
% 
% <<C:\Users\shafw\Downloads\Graph\2nd Ts 2.jpg>>
% 

%% Parametric Study Result Discussion 7
% 7. For Changing the values of T_secondary we notice that by increasing the value
% we decrease the Value of Entrainment ratio and vice versa.
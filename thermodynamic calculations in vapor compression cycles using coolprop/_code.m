clear all
close all

pyversion
[v,e] = pyversion; system([e,' -m pip uninstall -y CoolProp']);
[v,e] = pyversion; system([e,' -m pip install --user -U CoolProp']);

import py.CoolProp.CoolProp.*
%% all given values
Qe = 50 * 1000;
Teva = -20 + 273.15;                %Given in the question
Peva = PropsSI('P', 'T',Teva,'Q',1,'R134a');

Tcon = 40 + 273.15;                 %Given in the question
Pcon = PropsSI('P','T',Tcon,'Q',0,'R134a');
P9 = Peva;                          %Given in the question
P4 = Pcon;                          %Given in the question
P5 = P4;                            %No change of pressure in Condenser
P1 = P9;                            %No change of pressure in Evaporator
beta = 0.5;
Pi = (1-beta) * Peva + beta * Pcon; %This equation was given in the 
                                    %question and we get the 
                                    %Intermediate pressure from this

                                    
P2 = Pi;
P3 = Pi; 
P7 = Pi;
P6 = Pi;
P8 = Pi;                            %All of these points have intermediate
                                    % pressure the value of which we
                                    % determined previously

etas = 0.85;                        %Given value in the question
etame = 0.9;
Tsc = 35 + 273.15; 

%% Point 1
T1 = -15 + 273.15;
Tsh = T1;                           %Point 1 is at superheated region, 
                                    % so we can take T1 as Tsh

S1 = PropsSI('S','P',P1, 'T', T1, 'R134a');

H1 = PropsSI('H','P',P1,'T',T1,'R134a');
disp(T1);
disp(S1);
disp(P1);
disp(H1);
%% Point 2
S2s = S1;                                    %Isentropic compression in
                                             %compressor

H2s = PropsSI('H','P',P2, 'S', S2s,'R134a');
%Here we had to calculate 2s first, because the compressor is not ideal
%and has some efficiency. If the efficiency was 100% then we could have
%gotten point 2 directly, instead of 2s.

H2 = ((H2s-H1) / etas) + H1;
S2 = PropsSI('S','P', P2, 'H', H2, 'R134a');
T2 = PropsSI('T','P', P2, 'H', H2, 'R134a');
disp(S2s);
disp(H2s);
disp(T2);
disp(S2);
disp(P2);
disp(H2);
%% Point 3
T3 = PropsSI('T', 'P', P3, 'Q', 1, 'R134a'); %Its on the vapor line, so
H3 = PropsSI('H', 'P', P3, 'Q', 1, 'R134a'); % by only knowing the value
S3 = PropsSI('S', 'P', P3, 'Q', 1, 'R134a'); % of P3 we can get other values
disp(T3);
disp(S3);
disp(P3);
disp(H3);
%% Point 5
T5 = Tsc;                                    %Point 5 is at subcooled region, 
                                             % so we can take T5 as Tsc
S5 = PropsSI('S', 'P', P5, 'T', T5, 'R134a');
H5 = PropsSI('H','P',P5,'T',T5,'R134a');
disp(T5);
disp(S5);
disp(P5);
disp(H5);
%% Point 6
H6 = H5;                                     %Isenthalpic expansion in valve

T6 = PropsSI('T', 'P', P6, 'H', H6, 'R134a');
S6 = PropsSI('S', 'P', P6, 'H', H6, 'R134a');
disp(T6);
disp(S6);
disp(P6);
disp(H6);
%% Point 8
H8 = PropsSI('H', 'P', P8, 'Q', 0, 'R134a'); %Its on the liquid line, so
T8 = PropsSI('T', 'P', P8, 'Q', 0, 'R134a'); %by only knowing the value
S8 = PropsSI('S', 'P', P8, 'Q', 0, 'R134a'); %of P8 we can get other values
disp(T8);
disp(S8);
disp(P8);
disp(H8);
%% Point 9
H9 = H8;                                     %Isenthalpic expansion in valve
S9 = PropsSI('S', 'P', P9, 'H', H9, 'R134a');
T9 = PropsSI('T','P', P9, 'H', H9, 'R134a');
disp(T9);
disp(S9);
disp(P9);
disp(H9);

%% mass balance and energy balance to calculate mass flow
me = Qe/ (H1 - H9);                          %energy balance in the 
                                             %evaporator to calculate me

mi = (me*H8 - me*H6)/(H6 - H3);              % energy balance in the flash
                                             % in the flash intercooler

mc = me + mi;                                % mass balance at point 7

disp(me);
disp(mi);
disp(mc);

%% Point 7
H7 = ((H3 * mi) + (H2 * me)) / (mc);         %Energy Balance is used at
                                             %point 7 to calculate H7
T7 = PropsSI('T', 'P', P7, 'H', H7, 'R134a');
S7 = PropsSI('S', 'P', P7, 'H', H7, 'R134a');
disp(T7);
disp(S7);
disp(P7);
disp(H7);
%% Point 4
S4s = S7;                                       %isentropic compression in
                                                % compressor 
                                                
                                                
H4s = PropsSI('H', 'P', P4, 'S', S4s, 'R134a'); % for similar
                                                % reason as 2s,which was
                                                % previously mentioned,
                                                % we calculated 4s
H4 = ((H4s-H7) / etas) + H7;
T4 = PropsSI('T','P', P4, 'H', H4, 'R134a');
S4 = PropsSI('S','P', P4, 'H', H4, 'R134a');
disp(S4s);
disp(H4s);
disp(T4);
disp(S4);
disp(P4);
disp(H4);
%% Calculation of Work 
Qc = H4 - H5;
We1 = H2-H1;
We2 = H4 - H7;
disp(Qc);
disp(We1);
disp(We2);

%% Not used in calculation, only to plot the T-S diagram
Sref = PropsSI('S', 'P', P4, 'Q', 1, 'R134a');  
Tref = PropsSI('T', 'P', P4, 'Q', 1, 'R134a');

Sref1 = PropsSI('S', 'P', P4, 'Q', 0, 'R134a');
Tref1 = PropsSI('T', 'P', P4, 'Q', 0, 'R134a');

Sref2 = PropsSI('S', 'P', P9, 'Q', 1, 'R134a');
Tref2 = PropsSI('T', 'P', P9, 'Q', 1, 'R134a');
% These values (Sref,Tref, Sref1, Tref1, Sref2,Tref2) are not
% used anywhere in the calculation. they were only calculated to plot
% the graph, as the graph passes through these points in the fully
% liquid or fully vapor region

%% Code for plotting the diagrams
Pl=[];
Hl=[];
Tl=[];
Sl=[];

Pg=[];
Hg=[];
Tg=[];
Sg=[];

for P = 5000.0: 100000.0 :4059280.0
    H   = PropsSI( 'H','P', P , 'Q', 0, 'R134a');
    Pl=[Pl P];
    Hl=[Hl H];

end

for P = 5000.0: 100000.0 :4059280.0
    H   = PropsSI( 'H','P', P , 'Q', 1, 'R134a');
    Pg=[Pg P];
    Hg=[Hg H];

end

for T = 170.0: 1.0 : 374.0
    S   = PropsSI( 'S','T', T , 'Q', 0, 'R134a');
    Tl = [Tl T];
    Sl =[Sl S];
end

for T = 170.0: 1.0 : 374.0
    S   = PropsSI( 'S','T', T , 'Q', 1, 'R134a');
    Tg = [Tg T];
    Sg =[Sg S];
end

Hmat = [H1, H2, H7, H4, H5, H6, H3, H2, H8, H9, H1];
Pmat = [P1, P2, P7, P4, P5, P6, P3, P2, P8, P9, P1];
Smat = [S1, S2, S7, S4, Sref, Sref1, S5, S6, S3, S7, S3, S8, S9, Sref2, S1];
Tmat = [T1, T2, T7, T4, Tref, Tref1, T5, T6, T3, T7, T3, T8, T9, Tref2, T1];
figure
label1 = {'1', '2', '7', '4', '5', '6' , '3', '2', '8', '9', '1'};
plot(Hl,Pl)
xlabel('Specific enthalpy (J/kg)')
ylabel('Pressure')
title(' P-h diagram ')
hold all;
plot(Hg,Pg)
hold all;
plot(Hmat,Pmat)
text (Hmat, Pmat, label1)

figure
label2 = {'1', '2', '7', '4', '', '' , '5', '6', '3', '7', '3', '8', '9', '', '1'};
plot(Sl,Tl)
xlabel('Specific entropy (J/kg.K)')
ylabel('Temperature')
title(' T-s diagram ')
hold all;
plot(Sg,Tg)
hold all;
plot(Smat,Tmat)
text (Smat, Tmat, label2)

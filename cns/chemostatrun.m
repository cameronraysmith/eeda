%-----------------------------%
% chemostat model with nutrient switching
% depends on nutrientswitch.m
% nutrientswitch(Tfin, S0, n0, D, um, Ks, gam)
% nutrientswitch(100, 3, 5, 0.5, 1, 2.5, 1)
%
% VARIABLE NAME: DESCRIPTION [UNITS]
%
% Tfin: simulation runtime [time]
% S0: initial substrate concentration [mass/volume]
% n0: initial cell density [mass of organisms / volume]
% D: dilution constant = flow rate / reaction volume [1/time]
% um: maximum specific growth rate [1/time]
% Ks: half-saturation (Michaelis-Menten) constant [mass/volume]
% gam: mass of organisms formed / mass of substrate used [1/mass]
%
% Time measured in minutes
% 
% sucrose molecular weight 342.3 g/mol
% average sucrose concentrated stock 730mM or 250g/L
% average sucrose concentration in supplemented M9 media is 6mM or 2g/L
%
% glucose molecular weight 180.16 g/mol
% glucose concentrated stock at 20 mass% (i.e. 20g/100mL)
% 20 mL of 20% glucose in 1L of media is 22mM or 4g/L or 0.4%
% 
% One unit of OD600 corresponds to a cell dry weight of 0.3 g L-1
%-----------------------------%

Ts = 20; % switching period in minutes
T = 500; % total simulation time in minutes
nTi = T/Ts; % number of intervals

S0= 4; %initial substrate concentration in g/L
n0= 5; % initial cell density in g/L

D = 0.2; % 200 mL/min / 1000mL culture
um = 10;
Ks = 2;
gam = 1;

S=[];
N=[];
St = [];
Nt = [];

CarbonPulse = 3;

for i=Ts:Ts:T
    [St Nt]=nutrientswitch(Ts, S0, n0, D, um, Ks, gam);
    S = [S; St];
    N = [N; Nt];
    S0 = St(end) + CarbonPulse;
    n0 = Nt(end);
end

ydat = [S N];

linecolor = [0 0.5000 0.4000; 0.4921    0.7460    0.4000; 1.0000    1.0000    0.4000];
%%
figure(); hold on;
for i=1:2 plot(ydat(:,i),'Color',linecolor(i,:),'LineWidth',2); end
legend('substrate','cell density');
hold off;
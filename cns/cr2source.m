%-----------------------------%
% chemostat model with nutrient switching
% depends on ns2source.m
% ns2source(Tfin, S0, n0, D, um1, um2, Ks, gam)
% ns2source(100, 3, 5, 0.5, 1, 1, 2.5, 1)
%
% VARIABLE NAME: DESCRIPTION [UNITS]
%
% Tfin: simulation runtime [time]
% S0: initial substrate concentration [mass/volume]
% n0: initial cell density [wet mass of organisms / volume]
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
% One unit of OD600 corresponds to a cell dry weight of 0.3g/L
% a 1 L, overnight shaker-culture of E. coli with a 
% cell density of 3-4 x 10^9/ml corresponds to a pellet wet weight of approximately 3 g/liter.
% http://www.qiagen.com/plasmid/bacterialcultures.aspx
%
% Strain MW162 of E. coli K-12 
% Generation time in M9-0.4% glucose medium - 48 +- 1.6 minutes ~ .0144
% http://bionumbers.hms.harvard.edu/bionumber.aspx?s=y&id=104332&ver=3
%
% E. coli BL21 (DE3)
% Generation time in M9-0.4% glycerol medium - 102 +- 9 minutes ~ .006796
% http://bionumbers.hms.harvard.edu/bionumber.aspx?s=n&id=105395&ver=4
% 
% 0.05 min^-1 specific growth rate gives a doubling time of about 14 minutes
% 0.03 min^-1 specific growth rate gives a doubling time of about 23 minutes
% 0.015 min^-1 specific growth rate gives a doubling time of about 46 minutes
% 0.012 min^-1 specific growth rate gives a doubling time of about 58 minutes
%
%-----------------------------%
clear; close; 
set(0,'defaultaxesfontsize',30);
set(0,'defaulttextfontsize',30);

Ts = 50; % switching period in minutes
Tfin = 2000; % total simulation time in minutes
Tt = [];
T = [];

%% Artificial test case
% S10 = 3; %initial substrate 1 concentration in g/L
% S20 = 0; %initial substrate 1 concentration in g/L
% S1c = 1; %substrate 1 source concentration in g/L
% S2c = 1; %substrate 2 source concentration in g/L
% n0= 3; % initial cell density in g/L

% D = 0.02; % 20 mL/min / 1000mL culture
% um1 = .05; % 0.05 1/min gives a doubling time of about 14 minutes
% um2 = .05; % 0.05 1/min gives a doubling time of about 14 minutes
% Ks = 3;
% gam = 1;

%% Paramters for glucose/glycerol switching in M9
S10 = 3; %initial substrate 1 concentration in g/L
S20 = 0; %initial substrate 1 concentration in g/L
S1c = 0; %substrate 1 continuous source concentration in g/L
S2c = 0; %substrate 2 continuous source concentration in g/L
n0= 3; % initial cell density in g/L

D = 0.01; % 10 mL/min / 1000mL culture
um1 = .0144; % 0.014 1/min gives a doubling time of about 50 minutes
um2 = .006796; % 0.007 1/min gives a doubling time of about 99 minutes
Ks = 2;
gam = 1;

S1=[]; S2=[];
N=[];
S1t = []; S2t = [];
Nt = [];

CarbonPulse = 3;

for i=Ts:Ts:Tfin
    [S1t S2t Nt Tt]=ns2source(Ts, S10, S20, S1c, S2c, n0, D, um1, um2, Ks, gam);
    S1 = [S1; S1t];
    S2 = [S2; S2t];
    N = [N; Nt];
    T = [T; Tt+(i-Ts)];
    if mod(i/Ts,2)==0
        S10 = S1t(end) + CarbonPulse;
        S20 = S2t(end);
    else
        S10 = S1t(end);
        S20 = S2t(end) + CarbonPulse;
    end
    n0 = Nt(end);    
end

ydat = [S1+S1c S2+S2c N];

linecolor = [0 0.5000 0.4000; 0.4921    0.7460    0.4000; 1.0000    1.0000    0.4000];
linspec = {'-','--','-.'}; cls={['S1'],['S2'],['CD']}; 

h=figure(); set(h,'Color','w');
hold on;
for i=1:3 
    plot(T,ydat(:,i),char(linspec(i)),'Color',linecolor(i,:),'LineWidth',4); 
    text(Tfin/2,max(ydat(:,i)),char(cls(i)),'Color',linecolor(i,:));
end
ylabel('concentrations (g/L)'); xlabel('time (minutes)');
axis([0 Tfin 0 max(max(ydat))]);
hold off;

%-----------------------------%
% chemostat nutrient switching
% 
% usage: nutrientswitch(Tfin, S0, n0, D, um, Ks, gam)
% example parameters: nutrientswitch(100, 3, 5, 0.5, 1, 2.5, 1)
% 
% VARIABLE NAME: DESCRIPTION [UNITS]
%
% Tfin: simulation runtime [time]
% S0: initial substrate concentration [mass/volume]
% n0: initial cell density [number of organisms / volume]
% D: dilution constant = flow rate / reaction volume [1/time]
% um: maximum specific growth rate [1/time]
% Ks: half-saturation (Michaelis-Menten) constant [mass/volume]
% gam: mass of organisms formed / mass of substrate used [1/mass]
%
%----------------------------%

function [S1, S2, N, T] = ns2source(Tfin, S10, S20, S1c, S2c, n0, D, um1, um2, Ks, gam)

[T,Y] = ode45(@chswi,[0 Tfin],[S10 S20 n0],[],S1c,S2c,D,um1,um2,Ks,gam);
%fprintf('final substrate concentration: %0.2f \nfinal cell density (g/L): %0.2f \n',Y(end,1),Y(end,2));
%plot(T,Y,'-o');
%legend('substrate concentration','cell number');
S1 = Y(:,1);
S2 = Y(:,2);
N = Y(:,3);

function dy = chswi(t,y,S1c,S2c,D,um1,um2,Ks,gam)
S1 = y(1); S2 = y(2);
CD = y(3);

dS1_dt = S1c*D-D*S1-um1*S1/(Ks+S1)*CD/gam;
dS2_dt = S2c*D-D*S2-um2*S2/(Ks+S2)*CD/gam;
dCD_dt = CD*(um1*S1/(Ks+S1)+um2*S2/(Ks+S2))-D*CD;
dy = [dS1_dt; dS2_dt; dCD_dt];
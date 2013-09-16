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

function [S, N] = nutrientswitch(Tfin, S0, n0, D, um, Ks, gam)

[T,Y] = ode15s(@chswi,[0 Tfin],[S0 n0],[],S0,D,um,Ks,gam);
%fprintf('final substrate concentration: %0.2f \nfinal cell density (g/L): %0.2f \n',Y(end,1),Y(end,2));
%plot(T,Y,'-o');
%legend('substrate concentration','cell number');
S = Y(:,1);
N = Y(:,2);

function dy = chswi(t,y,S0,D,um,Ks,gam)
%y(1) = substrate
%y(2) = cell density
dy = zeros(2,1); 
dy(1) = S0*D-y(1)*D-um*y(1)/(Ks+y(1))*y(2)/gam;
%dy(1) = -y(1)*D-um*y(1)/(Ks+y(1))*y(2)/gam;
dy(2) = um*y(1)*y(2)/(Ks+y(1))-D*y(2);
close all;
clear all;

% plot growth rate from temperature
T = linspace(4,48,100);
gg = growthfromtemp(T);
figure(1);
plot(T,gg);

% plot population density from growth rate
tt = linspace(0,24);
dd = gompertz(tt,0,0,growthfromtemp(35),3*1.0071);
figure(2);
plot(tt,dd)
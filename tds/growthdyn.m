function [tt,dd,T,timehist] = growthdyn(lowtemp,hightemp,...
                                        tempperiod,dilperiod,...
                                        dilfactor,totalperiods,...
                                        initialdens,lagtime,...
                                        dt, tinit,...
                                        savefigFlag, opensavefigFlag)
%GROWTHDYN fluctuating temperature growth dynamics.
%   GROWTHDYN(lowtemp,hightemp,...
%             tempperiod,dilperiod,...
%              dilfactor,totalperiods,dt,...
%              initialdens,lagtime,...
%              savefigFlag, opensavefigFlag) 
%   simulates growth dynamics where the temperature
%   fluctuates between lowtemp and hightemp with period tempperiod.
%   Dilution occurs every dilperiod with factor dilfactor.
%   The simulation is run for a total time of tempperiod*totalperiods
%   
%   For example, 
%           GROWTHDYN(25,48,24,12,10,0.1,0,2,1,1)
%
%   See also

%   Author(s): Cameron Smith, 16 Sep 2013
%   $Revision: 0.0 $  $Date: 2013/09/17 19:16:21 $
%-------------------------------


%% Initialize parameters (uncomment to run as script)
%-------------------------------
% close all;
% clear all;
% 
% initialdens = 0; %initial density
% lagtime = 2; %hours
% lowtemp = 25; % degrees celsius
% hightemp = 48; % degrees celsius
% tempperiod = 24; % temperature switching time (hours)
% dilperiod = 12; % dilution time (hours)
% dilfactor = 10; % dilution factor
% totalperiods = 3; % total number of temperature switching periods
% dt = 0.1; % time step
% tinit = 0.01; % initial time
%-------------------------------

%% Compute time and temperature vectors

% Time vector
%-------------------------------
avgtemp = (hightemp+lowtemp)/2;
tt = linspace(tinit,tempperiod*totalperiods,(tempperiod*totalperiods-tinit)/dt);
%-------------------------------

% Mixed square-sinusoid temperature vector
%-------------------------------
% 0 => sinuosoid, 1 => square
%ssmix=0.9;

%Low then high
%T1 = avgtemp + (hightemp-avgtemp)*(-1).^(ceil(tt/tempperiod));
%T2 = avgtemp - (hightemp-avgtemp)*sin(tt*pi/tempperiod);

%High then low
%T1 = avgtemp + (hightemp-avgtemp)*(-1).^(floor(tt/tempperiod));
%T2 = avgtemp + (hightemp-avgtemp)*sin(tt*pi/tempperiod);

%T = ssmix*T1 + (1-ssmix)*T2;
%-------------------------------

% Trapezoidal temperature vector
%-------------------------------
T = avgtemp + (hightemp-avgtemp)*trapezoid(tt*2*pi/tempperiod,24);
%-------------------------------

%-----Print Summary-----
fprintf('\nSummary:\n------------------------\n')
fprintf('low temp: %0.0f, high temp: %0.0f\n',lowtemp,hightemp);
fprintf('low_t gr: %0.3f, high_t gr: %0.3f\n',growthfromtemp(lowtemp),growthfromtemp(hightemp));
fprintf('low_t d_gr: %0.3f, high_t d_gr: %0.3f\n',dgrowthfromtemp(lowtemp),dgrowthfromtemp(hightemp));
fprintf('temperature period: %0.0f, dilution period: %0.0f\n',tempperiod,dilperiod);
fprintf('dilution factor: %0.0f\n',dilfactor);
fprintf('time step: %0.5f, number of periods: %0.0f\n',dt,totalperiods);
fprintf('simulation time: %0.0f hours\n------------------------\n\n',tempperiod*totalperiods);
%-------------------------------

%% Main loop
% Vectorized
%-------------------------------
%dd = gompertz(tt,zeros(1,length(tt)),0,growthfromtemp(T),3*1.0071);
%-------------------------------

% Looped
%-------------------------------
toff = 0;
dd = zeros(1,length(tt));
dilcount=1;
for i=1:length(tt)
    if tt(i) >= dilcount*dilperiod
        dildens = dd(i-1)/dilfactor;
        fprintf('dilution at t=%0.2f, old dens=%0.2f, new dens=%0.2f, T=%0.0f\n',...
            tt(i),dd(i-1),dildens,T(i));
        toff = tgompertz(dildens,initialdens,lagtime,growthfromtemp(T(i)),3*1.0071);
        if toff >= 0
            toff = -tt(i) + toff;
        else
            fprintf('t_off less than 0: %0.3f\n',toff);
            toff = -tt(i);
        end
        dilcount=dilcount+1;
        dd(i-1) = dildens;
    end
    if i==1
        dd(i) = dt*dgompertz(tt(i)+toff,initialdens,lagtime,growthfromtemp(T(i)),3*1.0071);
    else
        dd(i) = dd(i-1) + dt*dgompertz(tt(i)+toff,initialdens,lagtime,growthfromtemp(T(i)),3*1.0071);
    end
    timehist(i) = tt(i)+toff;
end
%-------------------------------

%% Plotting

% Settings
%-------------------------------
set(0,'defaultaxesfontsize',20);
set(0,'defaulttextfontsize',20);
scrsz = get(0,'ScreenSize');
%-------------------------------

% Plot data
%-------------------------------
if opensavefigFlag
    h=figure('Visible','off','Position',[0 0 scrsz(3)/4.5 scrsz(4)/1.4]); set(h,'Color','w');
else
    h=figure('Position',[0 0 scrsz(3)/4.5 scrsz(4)/1.4]); set(h,'Color','w');
end
subplot(2,1,1);
plot(tt,dd,'-ko','LineWidth',4,...
    'MarkerFaceColor','r',...
    'MarkerEdgeColor','none',...
    'MarkerSize',5);
ylabel('optical density');
title(sprintf('T_{min}=%0.0f, T_{max}=%0.0f, t_{T}=%0.0f, t_{D}=%0.0f, %0.0fX',...
    lowtemp,hightemp,tempperiod,dilperiod,dilfactor))
subplot(2,1,2);
plot(tt,T,'-ko','LineWidth',4,...
    'MarkerFaceColor','r',...
    'MarkerEdgeColor','none',...
    'MarkerSize',5);
xlabel('time (hours)');
%ylabel(sprintf('temperature (%cC)',char(176)));
ylabel('temperature (\circC)');
%-------------------------------

% Save and display PDF version of figure
%-------------------------------
if savefigFlag
    figname = sprintf('Tmin%0.0fTmax%0.0ft_T%0.0ft_D%0.0fdil%0.0fX.pdf',...
        lowtemp,hightemp,tempperiod,dilperiod,dilfactor);
    fprintf('figure name: %s\n\n',figname);
    export_fig(figname,h);
    close(h);
    if opensavefigFlag
        system(sprintf('evince %s &',figname));
    end
end
%-------------------------------
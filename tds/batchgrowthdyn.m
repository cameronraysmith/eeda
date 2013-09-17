% function [tt,dd,T,timehist] = growthdyn(1:lowtemp,2:hightemp,...
%                                         3:tempperiod,4:dilperiod,...
%                                         5:dilfactor,6:totalperiods,...
%                                         7:initialdens,8:lagtime,...
%                                         9:dt, 10:tinit,...
%                                         11:savefigFlag, 12:opensavefigFlag)
%growthdyn(1,  2,  3,  4,  5, 6, 7, 8,   9,   10, 11, 12)

% Simulation notes (Temperature range)
%-------------------------------
% low   high    notes
% 33    45      gives equivalent growth rate 
% 20    47      gives equivalent growth rate
% 25    48      arbitrary test
% 35    45      arbitrary test
%-------------------------------

% Test simulations
%-------------------------------
%growthdyn(1,  2,  3,  4,  5, 6, 7, 8,   9,   10, 11, 12)
growthdyn(33, 45, 24, 12, 10, 3, 0, 2, 0.1, 0.01, 1, 1)
growthdyn(20, 47, 24, 12, 10, 3, 0, 2, 0.1, 0.01, 1, 1)
growthdyn(25, 48, 24, 12, 10, 3, 0, 2, 0.1, 0.01, 1, 1)
growthdyn(35, 45, 24, 12, 10, 3, 0, 2, 0.1, 0.01, 1, 1)

%-------------------------------
% analyze growth curves using the modified Gompertz
% function from:
% 1. Zwietering MH, Jongenburger I, Rombouts FM, ’t Riet K van (1990)
% Modeling of the bacterial growth curve. Applied and environmental
% microbiology 56:1875–81.
% 
% y = initial + carcap exp [-exp[linslope*exp[1]/carcap * (lagtime-time)+1]]
% Solve[y = n + c*e[-e[linslope*exp[1]/c*(l-t)+1]],t]
% t = lagtime + c/(m*exp(1))*(1 - ln(-ln((y-initial)/carcap)))
%
% initial
% lagtime = lambda ~ 300 -
% linslope = um ~ 0.05 - linear region slope
% carcap = A ~ 1.2 - carrying capacity


function d = gompertz(time,initial,lagtime,linslope,carcap)

d = initial + carcap * exp (-exp(linslope*exp(1)/carcap .* (lagtime-time)+1));
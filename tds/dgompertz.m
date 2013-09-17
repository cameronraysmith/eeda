% Derivative of the Gompertz function

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


% derivative computed with wolfram alpha
% http://goo.gl/eAjjOl
% m exp((m (e l-e t))/c-e^((m (e l-e t))/c+1)+2)

% derivative computed with matlab
% syms c m l t
% diff(c*exp(-exp(m*exp(1)/c*(l-t)+1)),t)
% (3060513257434037*m*exp((3060513257434037*m*(l - t))/(1125899906842624*c) + 1))/(1125899906842624*exp(exp((3060513257434037*m*(l - t))/(1125899906842624*c) + 1)))

function d = dgompertz(time,initial,lagtime,linslope,carcap)

t=time; l=lagtime; m=linslope; c=carcap;
e = exp(1);

d = m * exp((m*(e*l - e*t))/c-exp((m*(e*l - e*t))/c + 1) + 2);
%d = (3060513257434037*m*exp((3060513257434037*m*(l - t))/(1125899906842624*c) + 1))/(1125899906842624*exp(exp((3060513257434037*m*(l - t))/(1125899906842624*c) + 1)));
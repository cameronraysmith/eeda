function grsens = growthsenstemp(T)

% derivative of the growth rate vs temperature curve
% generated using matlab symbolic math toolbox
%
% syms T 
% diff(((T-276)*(1-exp(0.323352*(T-322))))^2,T)
%
T = T+273;
grsens = (2*T - 552)*(exp((2912495893419009*T)/9007199254740992 ...
    - 468911838840460449/4503599627370496) - 1)^2 ...
    + (2912495893419009*exp((2912495893419009*T)/9007199254740992 ...
    - 468911838840460449/4503599627370496)*(exp((2912495893419009*T)/9007199254740992 ...
    - 468911838840460449/4503599627370496) - 1)*(T - 276)^2)/4503599627370496;
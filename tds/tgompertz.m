function t = tgompertz(d,initial,lagtime,linslope,carcap)

t = lagtime + carcap/(linslope*exp(1))*(1 - log(-log((d-initial)/carcap)));

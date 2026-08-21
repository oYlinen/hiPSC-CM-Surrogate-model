function costo=fcosto(par, afit, tfit, t0, MDDS, C)

% global afit tfit t0 MDDS C

DDvero= afit;
DDteor= par(1)*exp(-(t0-tfit)/par(2))+ MDDS*tfit +C;

err=DDvero-DDteor;
costo= err'*err;
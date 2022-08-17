function params = getparams()
params.mp = 0.0185; %kg
params.mb = 0.0086; %kg
params.lp = 0.0245; %m
params.lb = 0.027;  %m
params.g = 9.81;   % m/s^2
params.d = 0.0024/1000; %Nmm/s/1000

params.M2 = 0.3778; %kg
params.M4 = 0.6366; %kg

params.mM = (params.M4-params.M2-2*params.mp-2*params.mb)/2;
params.mP = params.M2-2*params.mM-2*params.mp-2*params.mb;

params.k2 = 14.7467;
params.k4 = 24.3431;

params.c2 = 0.0185;
params.c4 = 0.0285;

MF = [params.M2 params.M4];
kF = [params.k2 params.k4];
cF = [params.c2 params.c4];

params.LsF = mean(MF.*params.g./kF);
params.ctF = mean(cF./MF);

params.Ls = 0.35;
params.L = 0.5;
params.W = 0.4;

params.JM33 = params.mM*(0.05^2+0.03^2)/12;
params.JP33 = params.mP*(params.L^2+params.W^2)/12;

params.tau = 0.0451/1000;  % Nmm
params.theta_s = -0.1920;  % rad
params.theta_e = 0.2793;   % rad
params.epsilon = 0.05;

params.time = 0;



end
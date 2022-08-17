function params = getMassParams(params,N)
%% Report params:
params.MT = params.mP + N*(params.mM+params.mb+params.mp);
params.MR = params.mb*params.lb-params.mp*params.lp;
params.MR2 = params.mb*params.lb^2+params.mp*params.lp^2;
params.MMT = params.mM+params.mb+params.mp;
params.JT = params.JP33+N*params.JM33;

params.k = (params.MT*params.g)/params.Ls;
params.c = params.ctF*params.MT;


params.kt = params.k*(params.L^2+params.W^2);
params.ct = params.c*(params.L^2+params.W^2);
end
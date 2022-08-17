prepEnv
load('../data/tang_xyphi_3_HBEREF_Ham_2/sim_1')

sF = s.results.analysis.fourier.Fzkn(1:6,1:3);
somega = s.results.analysis.fourier.omega
sA = 2*abs(sF)
sP = (angle(sF)+0.5*pi)
%% Set some data
if(~exist('temp','dir'))
    mkdir('temp')
end
addpath('temp')
Nt = 3;
h = HBEOptim('three_tang_inphase');
h.setModel(Nt);

R = 0.5*ones(Nt,1);
gamma = (0:(1/Nt):1-1e-6).'*2*pi;
psi = gamma-0.5*pi;
h.setConfiguration(R,gamma,psi);

U_gains = -0.1*ones(Nt,1);
HStar = -0.0016*ones(Nt,1);
h.setEscapement(HStar,U_gains);
% R = s.R;
% gamma = s.gamma;
% psi = s.psi;

% Generate transfer function and derivative function 
h.getTransferFunctions;

%% Optimization search with in phase
h.setOrder(3);

h.setOrderToZero(2);

h.setPlatToZero(1,1);
h.setPlatToZero(2,1);
h.setPlatToZero(1,3);
h.setPlatToZero(2,3);
h.setMetSync(2,1,0)
h.setMetSync(3,1,0)
% h.setMetSync(4,1,0)
% h.setMetSync(5,1,0)
% h.setMetSync(6,1,0)

h.setPlatAmpToInit(3,1,0.01);
h.setPlatAmpToInit(3,3,0.002);
h.setPlatPhaseToInit(3,1,pi);
h.setPlatPhaseToInit(3,3,pi);

h.setMetAmpToInit(1,1,0.72);
h.setMetAmpToInit(1,3,0.01);

h.setMetPhaseToInit(1,3,0.2);
h.setOToInit(10.8);

h.assembleInitial();

[J,Ai,Pi,Oi] = h.evalError(h.rho_i,10^4)
%%
h.runOptimization(10^3);
%%
h.continueOptimization(10^4);
%%
h.continueOptimization(10^5);

rho_in = h.rho_o
h.save('desc_optims/')
%%
prepEnv
load('../data/tang_xyphi_3_HBEREF_Ham_1/sim_1')

sF = s.results.analysis.fourier.Fzkn(1:6,1:3);
somega = s.results.analysis.fourier.omega;
sA = 2*abs(sF)
sP = (angle(sF)+0.5*pi)
%% Set some data
if(~exist('temp','dir'))
    mkdir('temp')
end
addpath('temp')
Nt = 3;
h = HBEOptim('three_tang_distphase');
h.setModel(Nt);

R = 0.01*ones(Nt,1);
gamma = (0:(1/Nt):1-1e-6).'*2*pi;
psi = gamma-0.5*pi;
h.setConfiguration(R,gamma,psi);

U_gains = -0.1*ones(Nt,1);
HStar = -0.0016*ones(Nt,1);
h.setEscapement(HStar,U_gains);
% R = s.R;
% gamma = s.gamma;
% psi = s.psi;

% Generate transfer function and derivative function 
h.getTransferFunctions;

%% Optimization search with dist phase
h.setOrder(3);

h.setOrderToZero(2);

h.setPlatToZero(3,1);

h.setMetSync(2,1,2*pi/3*1)
h.setMetSync(3,1,2*pi/3*2)

h.setPlatAmpToInit(1,1,0.001);
h.setPlatAmpToInit(2,1,0.001);
h.setPlatAmpToInit(1,3,0.00002);
h.setPlatAmpToInit(2,3,0.00002);
h.setPlatAmpToInit(3,3,10^(-4));

h.setMetAmpToInit(1,1,0.72);
h.setMetAmpToInit(1,3,0.01);


h.setPlatPhaseToInit(1,1,0.5*pi);
h.setPlatPhaseToInit(2,1,pi);
h.setPlatPhaseToInit(1,3,pi);
h.setPlatPhaseToInit(2,3,pi);
h.setPlatPhaseToInit(3,3,0.2);

h.setMetPhaseToInit(1,3,0.1);
h.setOToInit(10.8);

h.assembleInitial();

[J,Ai,Pi,Oi] = h.evalError(h.rho_i,10^4)
%%
h.runOptimization(10^3);
%%
rho_dist = h.rho_o

h.save('desc_optims/')

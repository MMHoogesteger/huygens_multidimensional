%clear all
N = 3;
R = 0.1;
Aref = 0.72;
nreps = 1;
pert = 1e-3*(1:N).';
t0in = repmat(Aref,N,1)+pert;
%% tangential settings generation
simSettings.name = ['tang_xyphi_' num2str(N) '_HBEREF'];

simSettings.N = N;
simSettings.x = true;
simSettings.y = true;
simSettings.phi = true;

simSettings.R = repmat(R,N,1);
simSettings.gamma = (0:(1/N):1-1e-6).'*2*pi;
simSettings.psi = simSettings.gamma-0.5*pi;
simSettings.overridePsi = -0.5*pi; %tangent = -0.5pi / perpendicular = 0

simSettings.fs = 1e3;
simSettings.Tend = 600;
simSettings.Tss = 550;
simSettings.escapement = "Original";
simSettings.nReps = nreps;

simSettings.xyphi0 = zeros(3,nreps);

simSettings.Ntheta0 = [t0in];
                   
simSettings.hasVars = true;
simSettings.varNameCell = {'R'};
Rvec = [0.01,0.5];
simSettings.varVecCell{1} = mat2cell(repmat(Rvec,N,1),N,ones(1,numel(Rvec)));

save(['jobs/' simSettings.name])
%%
load('../data/tang_xyphi_3_HBEREF/sim_1')
t = s.results.solution.t;
Hmtheta1 = s.results.solution.Hm(:,4);
omega = s.results.analysis.fourier.omega;
Hmmean1 = fourierDC(t(t>550),Hmtheta1(t>550),omega);

HstarR001 = Hmmean1;
omegaR001 = omega;

load('../data/tang_xyphi_3_HBEREF/sim_2')
t = s.results.solution.t;
Hmtheta1 = s.results.solution.Hm(:,4);
omega = s.results.analysis.fourier.omega;
Hmmean1 = fourierDC(t(t>550),Hmtheta1(t>550),omega);

HstarR050 = Hmmean1;
omegaR050 = omega;

%% tangential settings generation
clear simSettings
simSettings.name = ['tang_xyphi_' num2str(N) '_HBEREF_Ham_1'];

simSettings.N = N;
simSettings.x = true;
simSettings.y = true;
simSettings.phi = true;
R = 0.01;
simSettings.R = repmat(R,N,1);
simSettings.gamma = (0:(1/N):1-1e-6).'*2*pi;
simSettings.psi = simSettings.gamma-0.5*pi;
simSettings.overridePsi = -0.5*pi; %tangent = -0.5pi / perpendicular = 0

simSettings.fs = 1e3;
simSettings.Tend = 1500;
simSettings.Tss = 1450;
simSettings.escapement = "Hamiltonian";
simSettings.Hstars = repmat(HstarR001,N,1)+0.03e-3;
simSettings.Ugains = repmat(-0.1,N,1);
simSettings.nReps = nreps;

simSettings.xyphi0 = zeros(3,nreps);

simSettings.Ntheta0 = [0.6;0.4;-0.5];
                   
simSettings.hasVars = false;

save(['jobs/' simSettings.name])

%% tangential settings generation
simSettings.name = ['tang_xyphi_' num2str(N) '_HBEREF_Ham_2'];

simSettings.N = N;
simSettings.x = true;
simSettings.y = true;
simSettings.phi = true;
R = 0.5;
simSettings.R = repmat(R,N,1);
simSettings.gamma = (0:(1/N):1-1e-6).'*2*pi;
simSettings.psi = simSettings.gamma-0.5*pi;
simSettings.overridePsi = -0.5*pi; %tangent = -0.5pi / perpendicular = 0

simSettings.fs = 1e3;
simSettings.Tend = 1500;
simSettings.Tss = 1450;
simSettings.escapement = "Hamiltonian";
simSettings.Hstars = repmat(HstarR050,N,1)+0.03e-3;
simSettings.Ugains = repmat(-0.1,N,1);
simSettings.nReps = nreps;

simSettings.xyphi0 = zeros(3,nreps);

simSettings.Ntheta0 = [0.6;0.4;-0.5];
                   
simSettings.hasVars = false;

save(['jobs/' simSettings.name])
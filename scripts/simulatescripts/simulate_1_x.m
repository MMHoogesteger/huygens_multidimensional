%simulate_1_x - Simulation to find reference and explore initial describing
% function analysis
%
% Author: M.M. Hoogesteger - Msc student of Mechanical Engineering
% Eindhoven University of Technology, Mechanical Engineering, Dynamics and
% control
% email address: m.m.hoogesteger@student.tue.nl
% May 2017; Last revision: 14-May-2018


clear all

prepEnv;

saveFolder = 'data\sim_1_xyphi_ref\';



%% Generate reference data by simulating one metronome on the platform
s = Simulation('1xyphi_original');
s.setModel(1,1,1,1);

% Prepare configuration
R = 0;
gamma = 0;
psi = 0;

s.setConfiguration(R,gamma,psi);
s.generateModel;

% Initial conditions and other settings
s.setSimParams(1000,1000)
s.setEscapement("Original")
s.setInitialConditionsZeroVelocity(1,[0;0;0],0.25*pi);

s.simulate(true);

s.analyze(900,true);

s.save(saveFolder);





%% Generate reference data for Hamiltonian escapement
t = s.results.solution.t;
Hmtheta1 = s.results.solution.Hm(:,2);
omega = s.results.analysis.fourier.omega;
Hmmean = fourierDC(t(t>900),Hmtheta1(t>900),omega);
%%
s = Simulation('1x_hamiltonian');
s.setModel(1,1,0,0);

% Prepare configuration
R = 0;
gamma = 0;
psi = 0;

s.setConfiguration(R,gamma,psi);
s.generateModel;

% Initial conditions and other settings
s.setSimParams(1000,1000)
s.setEscapement("Hamiltonian",Hmmean,-0.5);
s.setInitialConditionsZeroVelocity(1,0,0.25*pi);

s.simulate();
%%
s.analyze(900);

s.save(saveFolder);







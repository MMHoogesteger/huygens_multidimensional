%simulate_1_x - Simulation to find reference and explore initial describing
% function analysis
%
% Author: M.M. Hoogesteger - Msc student of Mechanical Engineering
% Eindhoven University of Technology, Mechanical Engineering, Dynamics and
% control
% email address: m.m.hoogesteger@student.tue.nl
% May 2017; Last revision: 14-May-2018

close all

prepEnv;

saveFolder = 'data\sim_2_xyphi_ref\';


%% Generate reference data by simulating one metronome on the platform
s = Simulation('2xyphi_original');
s.setModel(2,1,1,1);

% Prepare configuration
R = [-0.1;0.1];
gamma = [0;0];
psi = [0;0];

s.setConfiguration(R,gamma,psi);
s.generateModel;

% Initial conditions and other settings
s.setSimParams(1000,1000)
s.setEscapement("Original")
s.setInitialConditionsZeroVelocity(1,[0;0;0],[0.25*pi;0.23*pi]);

s.simulate(true);

s.analyze(900,true);

s.save(saveFolder);


%% Generate reference data for Hamiltonian escapement
t = s.results.solution.t;
Hmtheta1 = s.results.solution.Hm(:,2);
Hmtheta2 = s.results.solution.Hm(:,3);
omega = s.results.analysis.fourier.omega;
Hmmean1 = fourierDC(t(t>900),Hmtheta1(t>900),omega);
Hmmean2 = fourierDC(t(t>900),Hmtheta2(t>900),omega);
%%
s = Simulation('2x_hamiltonian');
s.setModel(2,1,0,0);

% Prepare configuration
R = [-0.1;0.1];
gamma = [0;0];
psi = [0;0];

s.setConfiguration(R,gamma,psi);
s.generateModel;

% Initial conditions and other settings
s.setSimParams(1000,1000)
s.setEscapement("Hamiltonian",[Hmmean1;Hmmean2]+0.03e-3,[-0.1;-0.1]);
s.setInitialConditionsZeroVelocity(1,0,[0.25*pi;0.23*pi]);

s.simulate();
%%
s.analyze(900);

s.save(saveFolder);


%% Plot some stuff

s.plotEnergy(1,0,150,100,savePlots,'1x_ham')
s.plotEnergyDeviations(1,0,1000,savePlots,'1x_ham')
s.plotCircleFit(1,800,savePlots,'1x_ham')

s.plotFourierAnalysis(1,savePlots,'1x')




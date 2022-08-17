%simulate_1_x - Simulation to find reference and explore initial describing
% function analysis
%
% Author: M.M. Hoogesteger - Msc student of Mechanical Engineering
% Eindhoven University of Technology, Mechanical Engineering, Dynamics and
% control
% email address: m.m.hoogesteger@student.tue.nl
% May 2017; Last revision: 16-July-2018

close all

prepEnv;

saveFolder = '../data/sim_linper_2_long/';

savePlots = false;



%% Plot settings and initialization




% To get correct figures in Matlab, not important
set(groot, 'defaulttextinterpreter','latex');  
set(groot, 'defaultAxesTickLabelInterpreter','latex');  
set(groot, 'defaultLegendInterpreter','latex');

% PGFplots line width same as axis
set(groot, 'DefaultLineLineWidth', 0.4);
% Color for background lines etc.
cGray = [0.4 0.4 0.4];

Rvec = 0.05:0.001:0.09;
Svec = cell(size(Rvec));
%% Generate reference data by simulating one metronome on the platform
for Rid = 1:length(Rvec)
    Rvar = Rvec(Rid);
    s = Simulation(['R_' num2str(Rid)]);
    s.setModel(2,1,1,1);

    % Prepare configuration
    R = [Rvar;Rvar];
    gamma = [pi;0];
    psi = [0.5*pi;-0.5*pi];

    s.setConfiguration(R,gamma,psi);
    s.generateModel;
    %%
    % Initial conditions and other settings
    s.setSimParams(1000,500)
    s.setEscapement("Original")
    s.setInitialConditionsZeroVelocity(11,zeros(3,11),[0.25*pi*ones(1,11);(-0.249999:0.0499:0.25)*pi]);

    s.simulate();
    %%
    s.analyze(450);

    s.save(saveFolder);

    Svec{Rid} = s;
end
%% Plot some stuff

Phi = zeros(11,length(Rvec));
for Rid = 1:length(Rvec)
    load([saveFolder 'R_' num2str(Rid) '.mat'])
    s.analyze(450);
    for idx = 1:11
        
        Phi(idx,Rid) = angle(s.results(idx).analysis.fourier.Fzkn(5,1));
    end
    
end



%%
Tend = 200;

s.plotEnergy(11,0,Tend,100,savePlots,'1x')
s.plotEnergyDeviations(11,0,Tend,savePlots,'1x')
s.plotCircleFit(11,280,savePlots,'1x')
%%
s.plotFourierAnalysis(11,savePlots,'1x')
    



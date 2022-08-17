clear all
%% linper
simSettings.name = 'linper_xyphi_2_R';

simSettings.N = 2;
simSettings.x = true;
simSettings.y = true;
simSettings.phi = true;

simSettings.R = [0;0];
simSettings.gamma = [pi;0];
simSettings.psi = [0.5*pi;-0.5*pi];

simSettings.fs = 1e3;
simSettings.Tend = 800;
simSettings.Tss = 750;
simSettings.escapement = "Original";
simSettings.nReps = 12;
simSettings.xyphi0 = zeros(3,12);
Aref = 0.72; %rad
simSettings.Ntheta0 = [Aref*ones(1,12);...
                       linspace(-Aref+1e-6,Aref-1e-6,12)];
                   
simSettings.hasVars = true;
simSettings.varNameCell = {'R'};
Rvec = 0.02:0.02:0.5;
simSettings.varVecCell{1} = mat2cell([Rvec;Rvec],2,ones(1,numel(Rvec)));

save(['jobs/' simSettings.name])
%
simSettings.name = 'linper_yphi_2_R';
simSettings.x = false;
simSettings.y = true;
simSettings.phi = true;
simSettings.xyphi0 = zeros(2,12);
save(['jobs/' simSettings.name])
%
simSettings.name = 'linper_xphi_2_R';
simSettings.x = true;
simSettings.y = false;
simSettings.phi = true;
save(['jobs/' simSettings.name])
%
simSettings.name = 'linper_xy_2_R';
simSettings.x = true;
simSettings.y = true;
simSettings.phi = false;
save(['jobs/' simSettings.name])
%
simSettings.name = 'linper_x_2_R';
simSettings.x = true;
simSettings.y = false;
simSettings.phi = false;
simSettings.xyphi0 = zeros(1,12);
save(['jobs/' simSettings.name])
%
simSettings.name = 'linper_y_2_R';
simSettings.x = false;
simSettings.y = true;
simSettings.phi = false;
save(['jobs/' simSettings.name])
%
simSettings.name = 'linper_phi_2_R';
simSettings.x = false;
simSettings.y = false;
simSettings.phi = true;
save(['jobs/' simSettings.name])

% %%
% simSettings.name = 'circtan_3_R';
% 
% simSettings.N = 3;
% simSettings.x = true;
% simSettings.y = true;
% simSettings.phi = true;
% 
% simSettings.R = [0;0;0];
% simSettings.gamma = [pi;pi/3;-pi/3];
% simSettings.psi = simSettings.gamma-0.5*pi;
% 
% simSettings.fs = 1e3;
% simSettings.Tend = 1000;
% simSettings.Tss = 950;
% simSettings.escapement = "Original";
% simSettings.nReps = 12;
% simSettings.xyphi0 = zeros(3,12);
% Aref = 0.72; %rad
% simSettings.Ntheta0 = [Aref*ones(1,12);...
%                        linspace(-Aref+1e-6,Aref-1e-6,12)];
%                    
% simSettings.hasVars = true;
% simSettings.varNameCell = {'R'};
% Rvec = 0.01:0.001:0.5;
% simSettings.varVecCell = mat2cell([Rvec;Rvec],2,ones(1,numel(Rvec)));
% 
% save(['jobs/' simSettings.name])
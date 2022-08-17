%clear all
N = 10;
R = 0.1;
Aref = 0.72;
nreps = 12;
pert = 1e-3*(1:N).';
t0in = repmat(Aref,N,1)+pert;
t0anti = Aref*((-1).^(1:N)).'+pert;
t0dist = linspace(-Aref,Aref,N).'+pert;
t0randinanti = Aref*(round(rand(N,2))*2-1)+pert;
t0rand = Aref*(rand(N,7)*2-1)+pert;
%% tangential settings generation
simSettings.name = ['per_xyphi_' num2str(N) '_R'];

simSettings.N = N;
simSettings.x = true;
simSettings.y = true;
simSettings.phi = true;

simSettings.R = repmat(R,N,1);
simSettings.gamma = (0:(1/N):1-1e-6).'*2*pi;
simSettings.psi = simSettings.gamma-0.5*pi;
simSettings.overridePsi = 0; %tangent = -0.5pi / perpendicular = 0 (psi = gamma + overridepsi)

simSettings.fs = 1e3;
simSettings.Tend = 600;
simSettings.Tss = 550;
simSettings.escapement = "Original";
simSettings.nReps = nreps;

simSettings.xyphi0 = zeros(3,nreps);

simSettings.Ntheta0 = [t0in,t0anti,t0dist,t0randinanti,t0rand];
                   
simSettings.hasVars = true;
simSettings.varNameCell = {'R'};
Rvec = [0.01,0.05:0.1:0.5];
simSettings.varVecCell{1} = mat2cell(repmat(Rvec,N,1),N,ones(1,numel(Rvec)));

save(['jobs/' simSettings.name])
%%
simSettings.name = 'linper_yphi_2_R';
simSettings.x = false;
simSettings.y = true;
simSettings.phi = true;
simSettings.xyphi0 = zeros(2,12);
save(['jobs/' simSettings.name])
%%
simSettings.name = 'linper_xphi_2_R';
simSettings.x = true;
simSettings.y = false;
simSettings.phi = true;
save(['jobs/' simSettings.name])
%%
simSettings.name = 'linper_xy_2_R';
simSettings.x = true;
simSettings.y = true;
simSettings.phi = false;
save(['jobs/' simSettings.name])
%%
simSettings.name = 'linper_x_2_R';
simSettings.x = true;
simSettings.y = false;
simSettings.phi = false;
simSettings.xyphi0 = zeros(1,12);
save(['jobs/' simSettings.name])
%%
simSettings.name = 'linper_y_2_R';
simSettings.x = false;
simSettings.y = true;
simSettings.phi = false;
save(['jobs/' simSettings.name])
%%
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
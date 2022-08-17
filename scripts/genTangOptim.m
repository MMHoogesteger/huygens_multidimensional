function [h] = genTangOptim(N,radius,name)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Set some data
if(~exist('temp','dir'))
    mkdir('temp')
end
addpath('temp')

h = HBEOptim(name);
h.setModel(N);

R = radius*ones(N,1);
gamma = (0:(1/N):1-1e-6).'*2*pi;
psi = gamma-0.5*pi;
h.setConfiguration(R,gamma,psi);

U_gains = -0.1*ones(N,1);
HStar = -0.0016*ones(N,1);
h.setEscapement(HStar,U_gains);
% R = s.R;
% gamma = s.gamma;
% psi = s.psi;


end


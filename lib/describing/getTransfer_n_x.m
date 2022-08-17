function [G0,G] = getTransfer_n_x(Nt,omega,r,params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Calculate transfer function
% Get mass parameters
Nq = Nt+1;
M = zeros(Nq,Nq);
M(1,1) = params.MT;

for tid = 1:Nt
    M(1,tid+1) = params.MR;
    M(tid+1,1) = params.MR;
    M(tid+1,tid+1) = params.MR2;
end

C = diag([params.c params.d*ones(1,Nt)]);
%C = zeros(3,3);
K = diag([params.k -params.MR*params.g*ones(1,Nt)]);

G = zeros(3,3,r)*omega;
for k = 1:r
    G(:,:,k) = inv(-omega^2*k^2*M + 1i*omega*k*C + K);
end
G0 = inv(K);
end


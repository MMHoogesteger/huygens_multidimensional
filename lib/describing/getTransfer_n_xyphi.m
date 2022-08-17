function [G0,G] = getTransfer_n_xyphi(Nt,omega,r,params,R,gamma,psi)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%% Calculate transfer function
% Get mass parameters
% Nq = Nt+1;
% M = zeros(Nq,Nq);
% M(1,1) = params.MT;
% 
% for tid = 1:Nt
%     M(1,tid+1) = params.MR;
%     M(tid+1,1) = params.MR;
%     M(tid+1,tid+1) = params.MR2;
% end
MT = params.MT;
MR = params.MR;
MMT = params.MMT;
MR2 = params.MR2;
JT = params.JT;

Nq = Nt+3;
M = zeros(Nq,Nq);
M(1,1) = MT;
M(2,2) = MT;

M(1,3) = -sum(MMT.*R.*sin(gamma));
M(3,1) = M(1,3);

M(2,3) = sum(MMT.*R.*cos(gamma));
M(3,2) = M(2,3);

M(3,3) = sum(MMT.*R.^2) + JT;

for tid = 1:Nt
    
    M(1,3+tid) = MR*cos(psi(tid));
    M(3+tid,1) = M(1,3+tid);

    M(2,3+tid) = MR*sin(psi(tid));
    M(3+tid,2) = M(2,3+tid);

    M(3,3+tid) = MR.*R(tid).*sin(psi(tid)-gamma(tid));
    M(3+tid,3) = M(3,3+tid);
    M(3+tid,3+tid) = MR2;
end


C = diag([params.c params.c params.ct params.d*ones(1,Nt)]);

K = diag([params.k params.k params.kt -params.MR*params.g*ones(1,Nt)]);

G = zeros(Nq,Nq,r)*omega;
for k = 1:r
    G(:,:,k) = inv(-omega^2*k^2*M + 1i*omega*k*C + K);
end
G0 = inv(K);
end


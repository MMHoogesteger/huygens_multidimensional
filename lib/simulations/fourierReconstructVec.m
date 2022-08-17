function [f,df,ddf] = fourierReconstructVec(t,a0,ak,omega)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Reconstruct time signals
Nq = size(ak,1);
f = zeros(Nq,numel(t));
df = zeros(Nq,numel(t));
ddf = zeros(Nq,numel(t));

for n = 1:Nq
    a = ak(n,:);
    f(n,:) = fourierReconstruct(t,a0(n),a,omega);
    df(n,:) = fourierReconstructDer(t,a,omega);
    ddf(n,:) = fourierReconstructDDer(t,a,omega);
end
end


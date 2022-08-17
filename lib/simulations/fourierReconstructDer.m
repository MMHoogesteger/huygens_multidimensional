function [df] = fourierReconstructDer(t,ck,omega)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dfk = zeros(numel(ck),numel(t));

for k = 1:numel(ck)
    dfk(k,:) = 1i*k*omega*(ck(k)*exp(1i*omega*k*t) - conj(ck(k))*exp(-1i*omega*k*t));

end
df = sum(dfk,1);

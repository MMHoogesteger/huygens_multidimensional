function [ddf] = fourierReconstructDDer(t,ck,omega)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
ddfk = zeros(numel(ck),numel(t));

for k = 1:numel(ck)
    ddfk(k,:) = -k^2*omega^2*(ck(k)*exp(1i*omega*k*t) + conj(ck(k))*exp(-1i*omega*k*t));

end
ddf = sum(ddfk,1);

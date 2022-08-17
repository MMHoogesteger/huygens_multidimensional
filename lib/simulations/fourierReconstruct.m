function [f] = fourierReconstruct(t,c0,ck,omega)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
fk = zeros(numel(ck),numel(t));

for k = 1:numel(ck)
    fk(k,:) = ck(k)*exp(1i*omega*k*t) + conj(ck(k))*exp(-1i*omega*k*t);

end
f = c0+sum(fk,1);

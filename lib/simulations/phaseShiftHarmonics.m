function [ak_shifted] = phaseShiftHarmonics(ak,phase_shift)
% a_ke^{jk\omega t} + a_k^*e^{-jk\omega t} = Amplitude*sin(k\omega t + phase)
c = exp(1i*phase_shift);
[m,n] = size(ak);
v = 1:n;
cv = c.^v;

ak_shifted = ak.*cv;
end


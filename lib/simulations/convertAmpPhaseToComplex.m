function [ak] = convertAmpPhaseToComplex(Amplitude,Phase)
% Converts the Amplitude / phase representation to a fourier coefficient so
% that:
% a_ke^{jk\omega t} + a_k^*e^{-jk\omega t} = Amplitude*sin(k\omega t + phase)
abs_ak =  abs(0.5*Amplitude);
angle = Phase - 0.5*pi;

ak = abs_ak.*(cos(angle)+1i*sin(angle));
end


function [J] = errDescnxyphi_r_six_tang(Nt,rho_f,Nr,G_func,dGdomega_func,params,U_gains,HStar,R,gamma,psi,W,rho_ids,rho_o,res)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Nq = 3+Nt;
% Put optimization variables in place
rho_f(rho_ids) = rho_o;

% Descale platform movements
omega = rho_f(end);
rho_s = reshape(rho_f(1:end-1),Nq*2,Nr);

A = rho_s(1:Nq,:);
P = rho_s(Nq+(1:Nq),:);
if(Nr>2)
A(:,3) = A(:,3)*1e-2;
end



ak = convertAmpPhaseToComplex(A,P);

%Hack: All metronomes same amplitude and shifted phase
% ak(5,:) = phaseShiftHarmonics(ak(4,:),2*pi/6);
% ak(6,:) = phaseShiftHarmonics(ak(4,:),2*pi/6*2);
% ak(7,:) = phaseShiftHarmonics(ak(4,:),2*pi/6*3);
% ak(8,:) = phaseShiftHarmonics(ak(4,:),2*pi/6*4);
% ak(9,:) = phaseShiftHarmonics(ak(4,:),2*pi/6*5);

% ak(5,:) = ak(4,:);
% ak(6,:) = ak(4,:);
% ak(7,:) = ak(4,:);
% ak(8,:) = ak(4,:);
% ak(9,:) = ak(4,:);

ak(5,:) = phaseShiftHarmonics(ak(4,:),2*pi/2);
ak(6,:) = ak(4,:);
ak(7,:) = phaseShiftHarmonics(ak(4,:),2*pi/2);
ak(8,:) = ak(4,:);
ak(9,:) = phaseShiftHarmonics(ak(4,:),2*pi/2);


res = 2^nextpow2(res)-1;

[Jk] = harmonic_n_xyphi_eval_r_fft(Nt,ak,omega,G_func,dGdomega_func,params,U_gains,HStar,R,gamma,psi,res);

Jwk = W*Jk;
J = sqrt(sum(abs(Jwk).^2));



end


function [J,trnablaJ] = errDescnxyphi(Nt,rho_f,G_func,dGdomega_func,params,U_gains,HStar,R,gamma,psi,W,rho_ids,rho_o,res)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Nq = 3+Nt;
% Put optimization variables in place
rho_f(rho_ids) = rho_o;

% Descale platform movements
rho_f(1) = rho_f(1)*1e-3;
rho_f(2) = rho_f(2)*1e-3;
rho_f(3) = rho_f(3)*1e-2;

A = rho_f(1:Nq);
P = rho_f(Nq+(1:Nq))*1e-2;

omega = rho_f(end);

a1 = convertAmpPhaseToComplex(A,P);

[J1,nablaJ1] = harmonic_n_xyphi_eval(Nt,a1,omega,G_func,dGdomega_func,params,U_gains,HStar,R,gamma,psi,res);

Jw1 = W*J1(:,1);
nablaJw1 = W*nablaJ1;

J = sqrt(sum(abs(Jw1).^2));
trnablaJ = (1/J*real(Jw1'*nablaJw1)).';
% Rescale nabla to correct for scaling
trnablaJ(1) = trnablaJ(1)*1e-3;
trnablaJ(2) = trnablaJ(2)*1e-3;
trnablaJ(3) = trnablaJ(3)*1e-2;
trnablaJ(Nq+(1:Nq)) = trnablaJ(Nq+(1:Nq))*1e-2;
trnablaJ = trnablaJ(rho_ids);

end


function [J,trnablaJ] = errDescnxyphi_r_lsq(Nt,rho_f,Nr,G_func,dGdomega_func,params,U_gains,HStar,R,gamma,psi,W,rho_ids,rho_o,res)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
Nq = 3+Nt;
% Put optimization variables in place
rho_f(rho_ids) = rho_o;

% Descale platform movements
omega = 1e-2*rho_f(end);
rho_s = reshape(rho_f(1:end-1),Nq*2,Nr);

A = rho_s(1:Nq,:);
P = rho_s(Nq+(1:Nq),:);
if(Nr>2)
    A(:,3) = A(:,3)*1e-2;
end


ak = convertAmpPhaseToComplex(A,P);

[Jk,nablaJk] = harmonic_n_xyphi_eval_r(Nt,ak,omega,G_func,dGdomega_func,params,U_gains,HStar,R,gamma,psi,res);

Jwk = W*Jk;
%nablaJk = nablaJk(:,[1:Nq*Nr,Nq*Nr+[1:3 5:Nq*Nr]]);
nablaJwk = W*nablaJk;

J = [real(Jwk);imag(Jwk)];

trnablaJ = nablaJwk;

trnablaJA = nablaJwk(:,1:Nq*Nr);
trnablaJP = nablaJwk(:,Nq*Nr+1:2*Nq*Nr);
trnablaJomega = trnablaJ(:,end);
if(Nr>2)
trnablaJA(:,2*Nq+1:3*Nq) = trnablaJA(:,2*Nq+1:3*Nq)*1e-2;
end
for nr = 1:Nr
    trnablaJ(:,((nr-1)*Nq*2+1):((nr-1)*Nq*2+Nq)) = trnablaJA(:,((nr-1)*Nq+1):((nr-1)*Nq+Nq));
    trnablaJ(:,((nr-1)*Nq*2+Nq+1):((nr)*Nq*2)) = trnablaJP(:,((nr-1)*Nq+1):((nr-1)*Nq+Nq));
end
trnablaJ(:,end+1) = trnablaJomega;
% Rescale nabla to correct for scaling
trnablaJ(:,end) = trnablaJ(:,end)*1e-2;
trnablaJ = trnablaJ(:,rho_ids);
trnablaJ = [real(trnablaJ);imag(trnablaJ)];




end


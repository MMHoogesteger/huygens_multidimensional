function [J,trnablaJ] = errDescnxyphi_r(Nt,rho_f,Nr,G_func,dGdomega_func,params,U_gains,HStar,R,gamma,psi,W,rho_ids,rho_o,res)
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


J = sqrt(sum(abs(Jwk).^2));
trnablaJ = (1/J*real(Jwk'*nablaJwk)).';

trnablaJA = trnablaJ(1:Nq*Nr);
trnablaJP = trnablaJ(Nq*Nr+1:2*Nq*Nr);
trnablaJomega = trnablaJ(end);
trnablaJAr = reshape(trnablaJA,Nq,[]);
if(Nr>2)
trnablaJAr(:,3) = trnablaJAr(:,3)*1e-2;
end
trnablaJPr = reshape(trnablaJP,Nq,[]);
trnablaJr = [trnablaJAr ;trnablaJPr];

trnablaJ = [reshape(trnablaJr,[],1);trnablaJomega];
% Rescale nabla to correct for scaling
trnablaJ(end) = trnablaJ(end)*1e-2;
trnablaJ = trnablaJ(rho_ids);




end


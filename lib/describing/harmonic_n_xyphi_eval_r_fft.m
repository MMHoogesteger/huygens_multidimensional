function [Jk] = harmonic_n_xyphi_eval_r_fft(Nt,akm,omega,G_func,dGdomega_func,params,U_gains,HStar,R,gamma,psi,res)
bDebug = false;
Nq = Nt+3;
tstep = (2*pi/omega)/res;
Nr = size(akm,2);
NA = Nr*Nq;
if(~(Nq==size(akm,1)))
    error('wrong number and fourier coefficients')
end
ak = reshape(akm,[],1);


t = 0:tstep:2*pi/omega;
% Reconstruct time signals
[q,dq,ddq] = fourierReconstructVec(t,zeros(Nq,1),akm,omega);
% Lift dimensions for matrix multiplication etc:
t = reshape(t,1,1,[]);
q = reshape(q,Nq,1,[]);
dq = reshape(dq,Nq,1,[]);
ddq = reshape(ddq,Nq,1,[]);
%% Describing functions
x = q(1,:,:);
dx = dq(1,:,:);
ddx = ddq(1,:,:);

y = q(2,:,:);
dy = dq(2,:,:);
ddy = ddq(2,:,:);

phi = q(3,:,:);
dphi = dq(3,:,:);
ddphi = ddq(3,:,:);

qtid = (1:Nt) + 3;
theta = q(qtid,:,:);
dtheta = dq(qtid,:,:);
ddtheta = ddq(qtid,:,:);

MT = params.MT;
MR = params.MR;
MMT = params.MMT;
MR2 = params.MR2;
JT = params.JT;
g = params.g;

%% Mass matrix parts
%x
Mxphi = -sum(MR.*sin(theta).*sin(psi+phi) + MMT.*R.*(sin(gamma+phi)-sin(gamma))).*ddphi;
Mxtheta = sum(MR.*(cos(theta).*cos(psi+phi)-cos(psi)).*ddtheta);
Mx  = Mxphi+Mxtheta;

%y
Myphi = sum(MR.*sin(theta).*cos(psi+phi) + MMT.*R.*(cos(gamma+phi)-cos(gamma))).*ddphi;
Mytheta = sum(MR.*(cos(theta).*sin(psi+phi)-sin(psi)).*ddtheta);
My = Myphi+Mytheta;

%phi
Mphix = -sum(MR.*sin(theta).*sin(psi+phi) + MMT.*R.*(sin(gamma+phi)-sin(gamma))).*ddx;
Mphiy = sum(MR.*sin(theta).*cos(psi+phi) + MMT.*R.*(cos(gamma+phi)-cos(gamma))).*ddy;
Mphiphi = sum(2.*MR.*R.*cos(psi-gamma).*sin(theta) + MR2.*sin(theta).^2).*ddphi;
Mphitheta = sum(MR.*R.*sin(psi-gamma).*(cos(theta)-1).*ddtheta);
Mphi = Mphix+Mphiy+Mphiphi+Mphitheta;

%theta
Mt = zeros(Nt,1,numel(t));
for tid = 1:Nt
    Mtx = MR.*(cos(theta(tid,:,:)).*cos(psi(tid)+phi)-cos(psi(tid))).*ddx;
    Mty = MR.*(cos(theta(tid,:,:)).*sin(psi(tid)+phi)-sin(psi(tid))).*ddy;
    Mtphi = MR.*R(tid).*sin(psi(tid)-gamma(tid)).*(cos(theta(tid,:,:))-1).*ddphi;
    Mt(tid,1,:) = Mtx+Mty+Mtphi;
end

%% Nonlinear parts
Fx = -MMT*sum(dphi.^2.*R.*cos(gamma+phi))...
    - MR*sum(dtheta.^2.*sin(theta).*cos(psi+phi)...
    + 2.*dtheta.*dphi.*cos(theta).*sin(psi+phi)...
    + dphi.^2.*sin(theta).*cos(psi+phi));

Fy = -MMT*sum(dphi.^2.*R.*sin(gamma+phi))...
    + MR*sum(-dtheta.^2.*sin(theta).*sin(psi+phi)...
    + 2.*dtheta.*dphi.*cos(theta).*cos(psi+phi)...
    - dphi.^2.*sin(theta).*sin(psi+phi));

Fphi = MR*sum(-dtheta.^2.*R.*sin(psi-gamma).*sin(theta)...
    + 2.*dphi.*dtheta.*R.*cos(psi-gamma).*cos(theta))...
    + 2*MR2.*sum(dphi.*dtheta.*sin(theta).*cos(theta));


Ft = zeros(Nt,1,numel(t));
for tid = 1:Nt
    Ft(tid,1,:) = -MR.*dphi.^2.*R(tid).*cos(psi(tid)-gamma(tid)).*cos(theta(tid,:,:))...
        - MR2.*dphi.^2.*sin(theta(tid,:,:)).*cos(theta(tid,:,:))...
        - MR.*g.*(sin(theta(tid,:,:))-theta(tid,:,:));
end



Ux  = 0*t;
Uy  = 0*t;
Uphi  = 0*t;
Ut = U_gains.*dtheta.*(0.5*MR2.*dtheta.^2 + MR*g*cos(theta) - HStar);

M = -[Mx;My;Mphi;Mt];
F = -[Fx;Fy;Fphi;Ft];
U = [Ux;Uy;Uphi;Ut];
Ftot = M+F+U;

Ftot_fft = fft(Ftot,[],3)./length(t);
for nr = 1:Nr
    ekc = exp(-1i*nr*omega*t);
    %Dk(:,nr) = omega/(2*pi)*trapz(Ftot.*ekc,3)*tstep;
    Dk(:,nr) = Ftot_fft(:,:,nr+1);
    if(bDebug)
        DkM(:,nr) = omega/(2*pi)*trapz(M.*ekc,3)*tstep;
        DkF(:,nr) = omega/(2*pi)*trapz(F.*ekc,3)*tstep;
        DkU(:,nr) = omega/(2*pi)*trapz(U.*ekc,3)*tstep;
    end
end




G = zeros(Nq,Nq,Nr);

for nr = 1:Nr
    G(:,:,nr) = G_func(omega*nr);
end
%% Calculate residual
for nr = 1:Nr
    
    nids = (nr-1)*Nq+(1:Nq).';
    Jk(nids,1) = akm(:,nr)-G(:,:,nr)*Dk(:,nr);
    
end
end


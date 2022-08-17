function [Jk,nablaJk] = harmonic_n_xyphi_eval_r(Nt,akm,omega,G_func,dGdomega_func,params,U_gains,HStar,R,gamma,psi,res)
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


for nr = 1:Nr
    ekc = exp(-1i*nr*omega*t);
    Dk(:,nr) = omega/(2*pi)*trapz(Ftot.*ekc,3)*tstep;
    if(bDebug)
        DkM(:,nr) = omega/(2*pi)*trapz(M.*ekc,3)*tstep;
        DkF(:,nr) = omega/(2*pi)*trapz(F.*ekc,3)*tstep;
        DkU(:,nr) = omega/(2*pi)*trapz(U.*ekc,3)*tstep;
    end
end



%%
A = 2*abs(ak);
P = angle(ak)+0.5*pi;
t = reshape(t,1,1,[]);
%decomega = -1i*t.*ekc;

% derivative of a_k ii denotes only diagonal entries, rest is zero
iidadA = ak./A;
iidadP = 1i*ak;
iidadP(4) = NaN;
dadomega = 0*ak;
nablaa = [diag(iidadA),diag(iidadP),dadomega];

% Derivatives for amplitude and phase (works only for 1 harmonic)
iidqdA = zeros(Nq,Nr,numel(t));
iiddqdA = zeros(Nq,Nr,numel(t));
iidddqdA = zeros(Nq,Nr,numel(t));

iidqdP = zeros(Nq,Nr,numel(t));
iiddqdP = zeros(Nq,Nr,numel(t));
iidddqdP = zeros(Nq,Nr,numel(t));

dqdomega = zeros(Nq,1,numel(t));
ddqdomega = zeros(Nq,1,numel(t));
dddqdomega = zeros(Nq,1,numel(t));

for nq = 1:Nq
    for nr = 1:Nr
        nqa = (nr-1)*Nq+nq;
        
        Sqa = sin(nr.*omega.*t+P(nqa));
        Cqa = cos(nr.*omega.*t+P(nqa));
        
        ASqa = A(nqa).*Sqa;
        ACqa = A(nqa).*Cqa;
        
        iidqdA(nq,nr,:) = Sqa;
        iiddqdA(nq,nr,:) = nr.*omega.*Cqa;
        iidddqdA(nq,nr,:) = -nr.^2.*omega.^2.*Sqa;

        iidqdP(nq,nr,:) = ACqa;
        iiddqdP(nq,nr,:) = -nr.*omega.*ASqa;
        iidddqdP(nq,nr,:) = -nr.^2.*omega.^2.*ACqa;
        
        dqdomega(nq,1,:) = dqdomega(nq,1,:) + ACqa.*nr.*t;
        ddqdomega(nq,1,:) = ddqdomega(nq,1,:) + ACqa.*nr - ASqa.*nr.^2.*omega.*t;
        dddqdomega(nq,1,:) = dddqdomega(nq,1,:) - nr.^2.*omega.*(2.*ASqa + nr.*omega.*t.*Cqa);
    end
end



% Mass x
dMxdx = zeros(1,1,numel(t));
dMxddx = zeros(1,1,numel(t));
dMxdddx = zeros(1,1,numel(t));

dMxdy = zeros(1,1,numel(t));
dMxddy = zeros(1,1,numel(t));
dMxdddy = zeros(1,1,numel(t));

dMxdphi = zeros(1,1,numel(t));
dMxddphi = zeros(1,1,numel(t));
dMxdddphi = zeros(1,1,numel(t));

dMxdt = zeros(1,Nt,numel(t));
dMxddt = zeros(1,Nt,numel(t));
dMxdddt = zeros(1,Nt,numel(t));

% Mass y
dMydx = zeros(1,1,numel(t));
dMyddx = zeros(1,1,numel(t));
dMydddx = zeros(1,1,numel(t));

dMydy = zeros(1,1,numel(t));
dMyddy = zeros(1,1,numel(t));
dMydddy = zeros(1,1,numel(t));

dMydphi = zeros(1,1,numel(t));
dMyddphi = zeros(1,1,numel(t));
dMydddphi = zeros(1,1,numel(t));

dMydt = zeros(1,Nt,numel(t));
dMyddt = zeros(1,Nt,numel(t));
dMydddt = zeros(1,Nt,numel(t));

% Mass phi
dMphidx = zeros(1,1,numel(t));
dMphiddx = zeros(1,1,numel(t));
dMphidddx = zeros(1,1,numel(t));

dMphidy = zeros(1,1,numel(t));
dMphiddy = zeros(1,1,numel(t));
dMphidddy = zeros(1,1,numel(t));

dMphidphi = zeros(1,1,numel(t));
dMphiddphi = zeros(1,1,numel(t));
dMphidddphi = zeros(1,1,numel(t));

dMphidt = zeros(1,Nt,numel(t));
dMphiddt = zeros(1,Nt,numel(t));
dMphidddt = zeros(1,Nt,numel(t));

% Mass theta

dMtdx = zeros(Nt,1,numel(t));
dMtddx = zeros(Nt,1,numel(t));
dMtdddx = zeros(Nt,1,numel(t));

dMtdy = zeros(Nt,1,numel(t));
dMtddy = zeros(Nt,1,numel(t));
dMtdddy = zeros(Nt,1,numel(t));

dMtdphi = zeros(Nt,1,numel(t));
dMtddphi = zeros(Nt,1,numel(t));
dMtdddphi = zeros(Nt,1,numel(t));

dMtdt = zeros(Nt,Nt,numel(t));
dMtddt = zeros(Nt,Nt,numel(t));
dMtdddt = zeros(Nt,Nt,numel(t));

% F x
dFxdx = zeros(1,1,numel(t));
dFxddx = zeros(1,1,numel(t));
dFxdddx = zeros(1,1,numel(t));

dFxdy = zeros(1,1,numel(t));
dFxddy = zeros(1,1,numel(t));
dFxdddy = zeros(1,1,numel(t));

dFxdphi = zeros(1,1,numel(t));
dFxddphi = zeros(1,1,numel(t));
dFxdddphi = zeros(1,1,numel(t));

dFxdt = zeros(1,Nt,numel(t));
dFxddt = zeros(1,Nt,numel(t));
dFxdddt = zeros(1,Nt,numel(t));

% F y
dFydx = zeros(1,1,numel(t));
dFyddx = zeros(1,1,numel(t));
dFydddx = zeros(1,1,numel(t));

dFydy = zeros(1,1,numel(t));
dFyddy = zeros(1,1,numel(t));
dFydddy = zeros(1,1,numel(t));

dFydphi = zeros(1,1,numel(t));
dFyddphi = zeros(1,1,numel(t));
dFydddphi = zeros(1,1,numel(t));

dFydt = zeros(1,Nt,numel(t));
dFyddt = zeros(1,Nt,numel(t));
dFydddt = zeros(1,Nt,numel(t));

% F phi
dFphidx = zeros(1,1,numel(t));
dFphiddx = zeros(1,1,numel(t));
dFphidddx = zeros(1,1,numel(t));

dFphidy = zeros(1,1,numel(t));
dFphiddy = zeros(1,1,numel(t));
dFphidddy = zeros(1,1,numel(t));

dFphidphi = zeros(1,1,numel(t));
dFphiddphi = zeros(1,1,numel(t));
dFphidddphi = zeros(1,1,numel(t));

dFphidt = zeros(1,Nt,numel(t));
dFphiddt = zeros(1,Nt,numel(t));
dFphidddt = zeros(1,Nt,numel(t));

% F theta

dFtdx = zeros(Nt,1,numel(t));
dFtddx = zeros(Nt,1,numel(t));
dFtdddx = zeros(Nt,1,numel(t));

dFtdy = zeros(Nt,1,numel(t));
dFtddy = zeros(Nt,1,numel(t));
dFtdddy = zeros(Nt,1,numel(t));

dFtdphi = zeros(Nt,1,numel(t));
dFtddphi = zeros(Nt,1,numel(t));
dFtdddphi = zeros(Nt,1,numel(t));

dFtdt = zeros(Nt,Nt,numel(t));
dFtddt = zeros(Nt,Nt,numel(t));
dFtdddt = zeros(Nt,Nt,numel(t));

% U x
dUxdx = zeros(1,1,numel(t));
dUxddx = zeros(1,1,numel(t));
dUxdddx = zeros(1,1,numel(t));

dUxdy = zeros(1,1,numel(t));
dUxddy = zeros(1,1,numel(t));
dUxdddy = zeros(1,1,numel(t));

dUxdphi = zeros(1,1,numel(t));
dUxddphi = zeros(1,1,numel(t));
dUxdddphi = zeros(1,1,numel(t));

dUxdt = zeros(1,Nt,numel(t));
dUxddt = zeros(1,Nt,numel(t));
dUxdddt = zeros(1,Nt,numel(t));

% U y
dUydx = zeros(1,1,numel(t));
dUyddx = zeros(1,1,numel(t));
dUydddx = zeros(1,1,numel(t));

dUydy = zeros(1,1,numel(t));
dUyddy = zeros(1,1,numel(t));
dUydddy = zeros(1,1,numel(t));

dUydphi = zeros(1,1,numel(t));
dUyddphi = zeros(1,1,numel(t));
dUydddphi = zeros(1,1,numel(t));

dUydt = zeros(1,Nt,numel(t));
dUyddt = zeros(1,Nt,numel(t));
dUydddt = zeros(1,Nt,numel(t));

% U phi
dUphidx = zeros(1,1,numel(t));
dUphiddx = zeros(1,1,numel(t));
dUphidddx = zeros(1,1,numel(t));

dUphidy = zeros(1,1,numel(t));
dUphiddy = zeros(1,1,numel(t));
dUphidddy = zeros(1,1,numel(t));

dUphidphi = zeros(1,1,numel(t));
dUphiddphi = zeros(1,1,numel(t));
dUphidddphi = zeros(1,1,numel(t));

dUphidt = zeros(1,Nt,numel(t));
dUphiddt = zeros(1,Nt,numel(t));
dUphidddt = zeros(1,Nt,numel(t));

% U theta

dUtdx = zeros(Nt,1,numel(t));
dUtddx = zeros(Nt,1,numel(t));
dUtdddx = zeros(Nt,1,numel(t));

dUtdy = zeros(Nt,1,numel(t));
dUtddy = zeros(Nt,1,numel(t));
dUtdddy = zeros(Nt,1,numel(t));

dUtdphi = zeros(Nt,1,numel(t));
dUtddphi = zeros(Nt,1,numel(t));
dUtdddphi = zeros(Nt,1,numel(t));

dUtdt = zeros(Nt,Nt,numel(t));
dUtddt = zeros(Nt,Nt,numel(t));
dUtdddt = zeros(Nt,Nt,numel(t));

%%
% M
dMxdphi = sum(MR.*sin(theta).*cos(psi+phi)...
    + MMT.*R.*cos(gamma+phi)).*ddphi...
    +sum(MR.*cos(theta).*sin(psi+phi).*ddtheta);
dMxdddphi = sum(MR.*sin(theta).*sin(psi+phi) + MMT.*R.*(sin(gamma+phi)-sin(gamma)));


dMydphi = sum(MR.*sin(theta).*sin(psi+phi)...
    + MMT.*R.*sin(gamma+phi)).*ddphi...
    -sum(MR.*cos(theta).*cos(psi+phi).*ddtheta);
dMydddphi = -sum(MR.*sin(theta).*cos(psi+phi) + MMT.*R.*(cos(gamma+phi)-cos(gamma)));


dMphidphi = sum(MR.*sin(theta).*cos(psi+phi) + MMT.*R.*cos(gamma+phi)).*ddx...
    + sum(MR.*sin(theta).*sin(psi+phi) + MMT.*R.*sin(gamma+phi)).*ddy;
dMphidddx = sum(MR.*sin(theta).*sin(psi+phi) + MMT.*R.*(sin(gamma+phi)-sin(gamma)));
dMphidddy = -sum(MR.*sin(theta).*cos(psi+phi) + MMT.*R.*(cos(gamma+phi)-cos(gamma)));
dMphidddphi = -sum(2.*MR.*R.*cos(psi-gamma).*sin(theta) + MR2.*sin(theta).^2);

% F

dFxdphi = -MMT*sum(dphi.^2.*R.*sin(gamma+phi))...
    + MR*sum(-(dtheta.^2+dphi.^2).*sin(theta).*sin(psi+phi)...
            + 2.*dtheta.*dphi.*cos(theta).*cos(psi+phi));
dFxddphi = 2*MMT*sum(dphi.*R.*cos(gamma+phi))...
    + 2*MR*sum(dphi.*sin(theta).*cos(psi+phi)...
            + dtheta.*cos(theta).*sin(psi+phi));

dFydphi = MMT*sum(dphi.^2.*R.*cos(gamma+phi))...
    + MR*sum((dtheta.^2+ dphi.^2).*sin(theta).*cos(psi+phi)...
            + 2.*dtheta.*dphi.*cos(theta).*sin(psi+phi));
dFyddphi = 2*MMT*sum(dphi.*R.*sin(gamma+phi))...
    + 2*MR*sum(dphi.*sin(theta).*sin(psi+phi)...
            - dtheta.*cos(theta).*cos(psi+phi));
        
dFphiddphi = -MR*sum(2.*dtheta.*R.*cos(psi-gamma).*cos(theta))...
    - 2*MR2.*sum(dtheta.*sin(theta).*cos(theta));


for tid = 1:Nt
    theta_t = theta(tid,:,:);
    dtheta_t = dtheta(tid,:,:);
    ddtheta_t = ddtheta(tid,:,:);
    


    % M
    dMxdt(1,tid,:) = MR.*cos(theta_t).*sin(psi(tid)+phi).*ddphi...
        + MR.*sin(theta_t).*cos(psi(tid)+phi).*ddtheta_t;
    dMxdddt(1,tid,:) = -MR.*(cos(theta_t).*cos(psi(tid)+phi)-cos(psi(tid)));
    
    
    dMydt(1,tid,:) = -MR.*cos(theta_t).*cos(psi(tid)+phi).*ddphi...
        + MR.*sin(theta_t).*sin(psi(tid)+phi).*ddtheta_t;
    dMydddt(1,tid,:) = -MR.*(cos(theta_t).*sin(psi(tid)+phi)-sin(psi(tid)));
    
    
    dMphidt(1,tid,:) = MR.*cos(theta_t).*sin(psi(tid)+phi).*ddx...
        - MR.*cos(theta_t).*cos(psi(tid)+phi).*ddy...
        -(2.*MR.*R(tid).*cos(psi(tid)-gamma(tid)).*cos(theta_t)...
        + 2.*MR2.*sin(theta_t).*cos(theta_t)).*ddphi...
        +MR.*R(tid).*sin(psi(tid)-gamma(tid)).*sin(theta_t).*ddtheta_t;
    
    dMphidddt(1,tid,:) = - MR.*R(tid).*sin(psi(tid)-gamma(tid)).*(cos(theta_t)-1);
    
    
    dMtdphi(tid,1,:) = MR.*cos(theta_t).*sin(psi(tid)+phi).*ddx...
                       -MR.*cos(theta_t).*cos(psi(tid)+phi).*ddy;
                  
    dMtdt(tid,tid,:) = MR.*sin(theta_t).*cos(psi(tid)+phi).*ddx...
                +MR.*sin(theta_t).*sin(psi(tid)+phi).*ddy...
                +MR.*R(tid).*sin(psi(tid)-gamma(tid)).*sin(theta_t).*ddphi;

    dMtdddx(tid,1,:) = - MR.*(cos(theta_t).*cos(psi(tid)+phi)-cos(psi(tid)));
    dMtdddy(tid,1,:) = - MR.*(cos(theta_t).*sin(psi(tid)+phi)-sin(psi(tid)));
    dMtdddphi(tid,1,:) = - MR.*R(tid).*sin(psi(tid)-gamma(tid)).*(cos(theta_t)-1);
    
    
    % F
    dFxdt(1,tid,:) = MR*((dtheta_t.^2+dphi.^2).*cos(theta_t).*cos(psi(tid)+phi)...
    - 2.*dtheta_t.*dphi.*sin(theta_t).*sin(psi(tid)+phi));
    dFxddt(1,tid,:) = MR*(2.*dtheta_t.*sin(theta_t).*cos(psi(tid)+phi)...
    + 2.*dphi.*cos(theta_t).*sin(psi(tid)+phi));

    dFydt(1,tid,:) = MR*((dtheta_t.^2+ dphi.^2).*cos(theta_t).*sin(psi(tid)+phi)...
    + 2.*dtheta_t.*dphi.*sin(theta_t).*cos(psi(tid)+phi));
    dFyddt(1,tid,:) = MR*(2.*dtheta_t.*sin(theta_t).*sin(psi(tid)+phi)...
    - 2.*dphi.*cos(theta_t).*cos(psi(tid)+phi));

    dFphidt(1,tid,:) = MR*(dtheta_t.^2.*R(tid).*sin(psi(tid)-gamma(tid)).*cos(theta_t)...
    + 2.*dphi.*dtheta_t.*R(tid).*cos(psi(tid)-gamma(tid)).*sin(theta_t))...
    + 2*MR2.*(dphi.*dtheta_t.*(sin(theta_t).^2 - cos(theta_t).^2));
    dFphiddt(1,tid,:) = 2*MR*sum(dtheta_t.*R(tid).*sin(psi(tid)-gamma(tid)).*sin(theta_t)...
    - dphi.*R(tid).*cos(psi(tid)-gamma(tid)).*cos(theta_t))...
    - 2*MR2.*sum(dphi.*sin(theta_t).*cos(theta_t));
    
    dFtdt(tid,tid,:) = -MR.*dphi.^2.*R(tid).*cos(psi(tid)-gamma(tid)).*sin(theta_t)...
        - MR2.*dphi.^2.*(sin(theta_t).^2 - cos(theta_t).^2)...
        + MR.*g.*(cos(theta_t)-1);
    dFtddphi(tid,1,:) = 2.*MR.*dphi.*R(tid).*cos(psi(tid)-gamma(tid)).*cos(theta_t)...
        + 2.*MR2.*dphi.*sin(theta_t).*cos(theta_t);
    
    dUtdt(tid,tid,:)  = -U_gains(tid).*dtheta_t.*MR.*g.*sin(theta_t);
    dUtddt(tid,tid,:) = U_gains(tid).*((3/2)*MR2*dtheta_t.^2 + MR*g*cos(theta_t) - HStar(tid));
    
end

%%
% M
dMxdq = [dMxdx,dMxdy,dMxdphi,dMxdt];
dMydq = [dMydx,dMydy,dMydphi,dMydt];
dMphidq = [dMphidx,dMphidy,dMphidphi,dMphidt];
dMtdq = [dMtdx,dMtdy,dMtdphi,dMtdt];

dMqdq = [dMxdq;dMydq;dMphidq;dMtdq];

dMxddq = [dMxddx,dMxddy,dMxddphi,dMxddt];
dMyddq = [dMyddx,dMyddy,dMyddphi,dMyddt];
dMphiddq = [dMphiddx,dMphiddy,dMphiddphi,dMphiddt];
dMtddq = [dMtddx,dMtddy,dMtddphi,dMtddt];

dMqddq = [dMxddq;dMyddq;dMphiddq;dMtddq];

dMxdddq = [dMxdddx,dMxdddy,dMxdddphi,dMxdddt];
dMydddq = [dMydddx,dMydddy,dMydddphi,dMydddt];
dMphidddq = [dMphidddx,dMphidddy,dMphidddphi,dMphidddt];
dMtdddq = [dMtdddx,dMtdddy,dMtdddphi,dMtdddt];

dMqdddq = [dMxdddq;dMydddq;dMphidddq;dMtdddq];

% F
dFxdq = [dFxdx,dFxdy,dFxdphi,dFxdt];
dFydq = [dFydx,dFydy,dFydphi,dFydt];
dFphidq = [dFphidx,dFphidy,dFphidphi,dFphidt];
dFtdq = [dFtdx,dFtdy,dFtdphi,dFtdt];

dFqdq = [dFxdq;dFydq;dFphidq;dFtdq];

dFxddq = [dFxddx,dFxddy,dFxddphi,dFxddt];
dFyddq = [dFyddx,dFyddy,dFyddphi,dFyddt];
dFphiddq = [dFphiddx,dFphiddy,dFphiddphi,dFphiddt];
dFtddq = [dFtddx,dFtddy,dFtddphi,dFtddt];

dFqddq = [dFxddq;dFyddq;dFphiddq;dFtddq];

dFxdddq = [dFxdddx,dFxdddy,dFxdddphi,dFxdddt];
dFydddq = [dFydddx,dFydddy,dFydddphi,dFydddt];
dFphidddq = [dFphidddx,dFphidddy,dFphidddphi,dFphidddt];
dFtdddq = [dFtdddx,dFtdddy,dFtdddphi,dFtdddt];

dFqdddq = [dFxdddq;dFydddq;dFphidddq;dFtdddq];

% U
dUxdq = [dUxdx,dUxdy,dUxdphi,dUxdt];
dUydq = [dUydx,dUydy,dUydphi,dUydt];
dUphidq = [dUphidx,dUphidy,dUphidphi,dUphidt];
dUtdq = [dUtdx,dUtdy,dUtdphi,dUtdt];

dUqdq = [dUxdq;dUydq;dUphidq;dUtdq];

dUxddq = [dUxddx,dUxddy,dUxddphi,dUxddt];
dUyddq = [dUyddx,dUyddy,dUyddphi,dUyddt];
dUphiddq = [dUphiddx,dUphiddy,dUphiddphi,dUphiddt];
dUtddq = [dUtddx,dUtddy,dUtddphi,dUtddt];

dUqddq = [dUxddq;dUyddq;dUphiddq;dUtddq];

dUxdddq = [dUxdddx,dUxdddy,dUxdddphi,dUxdddt];
dUydddq = [dUydddx,dUydddy,dUydddphi,dUydddt];
dUphidddq = [dUphidddx,dUphidddy,dUphidddphi,dUphidddt];
dUtdddq = [dUtdddx,dUtdddy,dUtdddphi,dUtdddt];

dUqdddq = [dUxdddq;dUydddq;dUphidddq;dUtdddq];


dFtotdq = dMqdq + dFqdq + dUqdq;
dFtotddq = dMqddq + dFqddq + dUqddq;
dFtotdddq = dMqdddq + dFqdddq + dUqdddq;

dFtotdA = zeros(Nq,NA,numel(t));
dFtotdP = zeros(Nq,NA,numel(t));

for nq = 1:Nq
    for nr = 1:Nr
        nqa = (nr-1)*Nq+nq;
        qid = nq;
        dFtotdA(:,nqa,:) = dFtotdq(:,qid,:).*iidqdA(qid,nr,:) + dFtotddq(:,qid,:).*iiddqdA(qid,nr,:) + dFtotdddq(:,qid,:).*iidddqdA(qid,nr,:);
        dFtotdP(:,nqa,:) = dFtotdq(:,qid,:).*iidqdP(qid,nr,:) + dFtotddq(:,qid,:).*iiddqdP(qid,nr,:) + dFtotdddq(:,qid,:).*iidddqdP(qid,nr,:);
        if(bDebug)
            dMqdA(:,nqa,:) = dMqdq(:,qid,:).*iidqdA(qid,nr,:) + dMqddq(:,qid,:).*iiddqdA(qid,nr,:) + dMqdddq(:,qid,:).*iidddqdA(qid,nr,:);
            dFqdA(:,nqa,:) = dFqdq(:,qid,:).*iidqdA(qid,nr,:) + dFqddq(:,qid,:).*iiddqdA(qid,nr,:) + dFqdddq(:,qid,:).*iidddqdA(qid,nr,:);
            dUqdA(:,nqa,:) = dUqdq(:,qid,:).*iidqdA(qid,nr,:) + dUqddq(:,qid,:).*iiddqdA(qid,nr,:) + dUqdddq(:,qid,:).*iidddqdA(qid,nr,:);

            dMqdP(:,nqa,:) = dMqdq(:,qid,:).*iidqdP(qid,nr,:) + dMqddq(:,qid,:).*iiddqdP(qid,nr,:) + dMqdddq(:,qid,:).*iidddqdP(qid,nr,:);
            dFqdP(:,nqa,:) = dFqdq(:,qid,:).*iidqdP(qid,nr,:) + dFqddq(:,qid,:).*iiddqdP(qid,nr,:) + dFqdddq(:,qid,:).*iidddqdP(qid,nr,:);
            dUqdP(:,nqa,:) = dUqdq(:,qid,:).*iidqdP(qid,nr,:) + dUqddq(:,qid,:).*iiddqdP(qid,nr,:) + dUqdddq(:,qid,:).*iidddqdP(qid,nr,:);
        end
    end
end

dFtotdomega = mtimesx(dFtotdq,dqdomega) + mtimesx(dFtotddq,ddqdomega) + mtimesx(dFtotdddq,dddqdomega);

% derivative of Gc
for nr = 1:Nr
    ekc = exp(-1i*nr*omega*t);
    nids = (nr-1)*Nq+(1:Nq).';
    dcdA(nids,:) = omega/(2*pi)*trapz(dFtotdA(:,:,:).*ekc,3)*tstep;
    dcdP(nids,:) = omega/(2*pi)*trapz(dFtotdP(:,:,:).*ekc,3)*tstep;
    
    dcdomega(nids,1) = 1/omega*(Dk(:,nr)-Ftot(:,:,end).*ekc(end))...
            + omega/(2*pi)*trapz(dFtotdomega.*ekc - Ftot.*1i.*nr.*t.*ekc,3)*tstep;
        
    G(:,:,nr) = G_func(omega*nr);
    dGdomega(:,:,nr) = dGdomega_func(omega*nr);
    
    dGcdA(nids,:) = G(:,:,nr)*dcdA(nids,:);
    dGcdP(nids,:) = G(:,:,nr)*dcdP(nids,:);
    dGcdomega(nids,1) = G(:,:,nr)*dcdomega(nids,1) + dGdomega(:,:,nr)*Dk(:,nr);
            
    if(bDebug)
        dDMdA(nids,:) = omega/(2*pi)*trapz(dMqdA(:,:,:).*ekc,3)*tstep;
        dDFdA(nids,:) = omega/(2*pi)*trapz(dFqdA(:,:,:).*ekc,3)*tstep;
        dDUdA(nids,:) = omega/(2*pi)*trapz(dUqdA(:,:,:).*ekc,3)*tstep;

        dDMdP(nids,:) = omega/(2*pi)*trapz(dMqdP(:,:,:).*ekc,3)*tstep;
        dDFdP(nids,:) = omega/(2*pi)*trapz(dFqdP(:,:,:).*ekc,3)*tstep;
        dDUdP(nids,:) = omega/(2*pi)*trapz(dUqdP(:,:,:).*ekc,3)*tstep;
    end
end
for nr = 1:Nr
    %% Calculate residual
    nids = (nr-1)*Nq+(1:Nq).';
    Jk(nids,1) = akm(:,nr)-G(:,:,nr)*Dk(:,nr);
    
end
nablaGc_nr = [dGcdA, dGcdP, dGcdomega];
nablaJk = [diag(iidadA),diag(iidadP),dadomega]-nablaGc_nr;
end


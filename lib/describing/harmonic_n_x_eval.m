function [J1,nablaJ1] = harmonic_n_x_eval(Nt,a1,omega,G_func,dGdomega_func,params,U_gains,HStar,res)
bDebug = false;
Nq = Nt+1;
tstep = (2*pi/omega)/res;

t = 0:tstep:2*pi/omega;
% Reconstruct time signals
[q,dq,ddq] = fourierReconstructVec(t,zeros(Nq,1),a1,omega);
% Lift dimensions for matrix multiplication etc:
t = reshape(t,1,1,[]);
q = reshape(q,Nq,1,[]);
dq = reshape(dq,Nq,1,[]);
ddq = reshape(ddq,Nq,1,[]);
%% Describing functions
x = q(1,:,:);
dx = dq(1,:,:);
ddx = ddq(1,:,:);

qtid = (1:Nt) + 1;
theta = q(qtid,:,:);
dtheta = dq(qtid,:,:);
ddtheta = ddq(qtid,:,:);

MR = params.MR;
MR2 = params.MR2;
g = params.g;

Mx  = sum(-MR*ddtheta.*(cos(theta)-1));
Mt = zeros(Nt,1,numel(t));
for tid = 1:Nt
    Mt(tid,1,:) = -MR*ddx.*(cos(theta(tid,:,:))-1);
end

Fx  = sum(MR*dtheta.^2.*sin(theta));
Ft = MR*g*(sin(theta)-theta);

Ux  = 0*t;
Ut = U_gains.*dtheta.*(0.5*MR2.*dtheta.^2 + MR*g*cos(theta) - HStar);

M = [Mx;Mt];
F = [Fx;Ft];
U = [Ux;Ut];
Ftot = M+F+U;

ec = exp(-1i*omega*t);
D1 = omega/(2*pi)*trapz(Ftot.*ec,3)*tstep;
if(bDebug)
D1M = omega/(2*pi)*trapz(M.*ec,3)*tstep;
D1F = omega/(2*pi)*trapz(F.*ec,3)*tstep;
D1U = omega/(2*pi)*trapz(U.*ec,3)*tstep;
end

%%
A = 2*abs(a1);
t = reshape(t,1,1,[]);
decdomega = -1i*t.*ec;

% derivative of a_k ii denotes only diagonal entries, rest is zero
iidadA = a1./A;
iidadP = 1i*a1;
iidadP(2) = NaN;
dadomega = 0*a1;
nablaa = [diag(iidadA),diag(iidadP),dadomega];

% Derivatives for amplitude and phase
iidqdA = q./A;
iiddqdA = dq./A;
iidddqdA = ddq./A;

iidqdP = dq./omega;
iiddqdP = -omega.*q;
iidddqdP = -omega.*dq;

dqdomega = dq.*t./omega;
ddqdomega = dq./omega -omega.*t.*q;
dddqdomega = -omega.*(2.*q + t.*dq);

% Mass
dMxdx = zeros(1,1,numel(t));
dMxddx = zeros(1,1,numel(t));
dMxdddx = zeros(1,1,numel(t));

dMxdt = zeros(1,Nt,numel(t));
dMxddt = zeros(1,Nt,numel(t));
dMxdddt = zeros(1,Nt,numel(t));


dMtdx = zeros(Nt,1,numel(t));
dMtddx = zeros(Nt,1,numel(t));
dMtdddx = zeros(Nt,1,numel(t));

dMtdt = zeros(Nt,Nt,numel(t));
dMtddt = zeros(Nt,Nt,numel(t));
dMtdddt = zeros(Nt,Nt,numel(t));

% F
dFxdx = zeros(1,1,numel(t));
dFxddx = zeros(1,1,numel(t));
dFxdddx = zeros(1,1,numel(t));

dFxdt = zeros(1,Nt,numel(t));
dFxddt = zeros(1,Nt,numel(t));
dFxdddt = zeros(1,Nt,numel(t));


dFtdx = zeros(Nt,1,numel(t));
dFtddx = zeros(Nt,1,numel(t));
dFtdddx = zeros(Nt,1,numel(t));

dFtdt = zeros(Nt,Nt,numel(t));
dFtddt = zeros(Nt,Nt,numel(t));
dFtdddt = zeros(Nt,Nt,numel(t));

% U
dUxdx = zeros(1,1,numel(t));
dUxddx = zeros(1,1,numel(t));
dUxdddx = zeros(1,1,numel(t));

dUxdt = zeros(1,Nt,numel(t));
dUxddt = zeros(1,Nt,numel(t));
dUxdddt = zeros(1,Nt,numel(t));


dUtdx = zeros(Nt,1,numel(t));
dUtddx = zeros(Nt,1,numel(t));
dUtdddx = zeros(Nt,1,numel(t));

dUtdt = zeros(Nt,Nt,numel(t));
dUtddt = zeros(Nt,Nt,numel(t));
dUtdddt = zeros(Nt,Nt,numel(t));


for tid = 1:Nt
    dMxdt(1,tid,:)  = MR*ddtheta(tid,:,:).*sin(theta(tid,:,:));
    dMxdddt(1,tid,:)= -MR.*(cos(theta(tid,:,:))-1);
    
    dMtdt(tid,tid,:)= MR*ddx.*sin(theta(tid,:,:));
    dMtdddx(tid,1,:)= -MR.*(cos(theta(tid,:,:))-1);
    
    dFxdt(1,tid,:)  = MR*dtheta(tid,:,:).^2.*cos(theta(tid,:,:));
    dFxddt(1,tid,:) = 2*MR*dtheta(tid,:,:).*sin(theta(tid,:,:));
    
    dFtdt(tid,tid,:)  = MR*g*(cos(theta(tid,:,:))-1);
    
    dUtdt(tid,tid,:)  = -U_gains(tid).*dtheta(tid,:,:).*MR.*g.*sin(theta(tid,:,:));
    dUtddt(tid,tid,:) = U_gains(tid).*((3/2)*MR2*dtheta(tid,:,:).^2 + MR*g*cos(theta(tid,:,:)) - HStar(tid));
    
end

dMqdq = [dMxdx,dMxdt;dMtdx,dMtdt];
dMqddq = [dMxddx,dMxddt;dMtddx,dMtddt];
dMqdddq = [dMxdddx,dMxdddt;dMtdddx,dMtdddt];

dFqdq = [dFxdx,dFxdt;dFtdx,dFtdt];
dFqddq = [dFxddx,dFxddt;dFtddx,dFtddt];
dFqdddq = [dFxdddx,dFxdddt;dFtdddx,dFtdddt];

dUqdq = [dUxdx,dUxdt;dUtdx,dUtdt];
dUqddq = [dUxddx,dUxddt;dUtddx,dUtddt];
dUqdddq = [dUxdddx,dUxdddt;dUtdddx,dUtdddt];

dFtotdq = dMqdq + dFqdq + dUqdq;
dFtotddq = dMqddq + dFqddq + dUqddq;
dFtotdddq = dMqdddq + dFqdddq + dUqdddq;

dFtotdA = zeros(Nq,Nq,numel(t));
dFtotdP = zeros(Nq,Nq,numel(t));
for qid = 1:Nq
    dFtotdA(:,qid,:) = dFtotdq(:,qid,:).*iidqdA(qid,:,:) + dFtotddq(:,qid,:).*iiddqdA(qid,:,:) + dFtotdddq(:,qid,:).*iidddqdA(qid,:,:);
    dFtotdP(:,qid,:) = dFtotdq(:,qid,:).*iidqdP(qid,:,:) + dFtotddq(:,qid,:).*iiddqdP(qid,:,:) + dFtotdddq(:,qid,:).*iidddqdP(qid,:,:);
    if(bDebug)
    dMqdA(:,qid,:) = dMqdq(:,qid,:).*iidqdA(qid,:,:) + dMqddq(:,qid,:).*iiddqdA(qid,:,:) + dMqdddq(:,qid,:).*iidddqdA(qid,:,:);
    dFqdA(:,qid,:) = dFqdq(:,qid,:).*iidqdA(qid,:,:) + dFqddq(:,qid,:).*iiddqdA(qid,:,:) + dFqdddq(:,qid,:).*iidddqdA(qid,:,:);
    dUqdA(:,qid,:) = dUqdq(:,qid,:).*iidqdA(qid,:,:) + dUqddq(:,qid,:).*iiddqdA(qid,:,:) + dUqdddq(:,qid,:).*iidddqdA(qid,:,:);
    
    dMqdP(:,qid,:) = dMqdq(:,qid,:).*iidqdP(qid,:,:) + dMqddq(:,qid,:).*iiddqdP(qid,:,:) + dMqdddq(:,qid,:).*iidddqdP(qid,:,:);
    dFqdP(:,qid,:) = dFqdq(:,qid,:).*iidqdP(qid,:,:) + dFqddq(:,qid,:).*iiddqdP(qid,:,:) + dFqdddq(:,qid,:).*iidddqdP(qid,:,:);
    dUqdP(:,qid,:) = dUqdq(:,qid,:).*iidqdP(qid,:,:) + dUqddq(:,qid,:).*iiddqdP(qid,:,:) + dUqdddq(:,qid,:).*iidddqdP(qid,:,:);
    end
end
dFtotdomega = mtimesx(dFtotdq,dqdomega) + mtimesx(dFtotddq,ddqdomega) + mtimesx(dFtotdddq,dddqdomega);

% derivative of Gc

dcdA = omega/(2*pi)*trapz(dFtotdA.*ec,3)*tstep;
if(bDebug)
dDMdA = omega/(2*pi)*trapz(dMqdA.*ec,3)*tstep;
dDFdA = omega/(2*pi)*trapz(dFqdA.*ec,3)*tstep;
dDUdA = omega/(2*pi)*trapz(dUqdA.*ec,3)*tstep;

dDMdP = omega/(2*pi)*trapz(dMqdP.*ec,3)*tstep;
dDFdP = omega/(2*pi)*trapz(dFqdP.*ec,3)*tstep;
dDUdP = omega/(2*pi)*trapz(dUqdP.*ec,3)*tstep;
end
dcdP = omega/(2*pi)*trapz(dFtotdP.*ec,3)*tstep;
dcdomega = 1/omega*(D1-Ftot(:,:,end).*ec(end))...
           + omega/(2*pi)*trapz(dFtotdomega.*ec + Ftot.*decdomega,3)*tstep;
c = D1;
G = G_func(omega);
dGdomega = dGdomega_func(omega);

dGcdA = G*dcdA;
dGcdP = G*dcdP;
dGcdomega = G*dcdomega + dGdomega*c;

nablaGc = [dGcdA, dGcdP, dGcdomega];


%% Calculate residual

J1 = a1-G*D1;

nablaJ1 = nablaa-nablaGc;

end


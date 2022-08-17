%% Check eqs of motion
%% Circle mass matrix
clear all
n =1;

theta = sym('theta',[1,n]);
psi = sym('psi',[1,n]);
gamma = sym('gamma',[1,n]);

syms phi
syms MT MMT MR MR2 R JT

M = sym(zeros(3+n,3+n));
M(1,1) = MT;
M(2,2) = MT;
M(1,3) = -MR*sum(sin(theta).*sin(psi+phi));
M(3,1) = M(1,3);

M(2,3) = MR*sum(sin(theta).*cos(psi+phi));
M(3,2) = M(2,3);
M(3,3) = n*MMT*R^2 + MR2*sum(sin(theta).^2) + JT;

for met_k = 1:n
    M(1,3+met_k) = MR*cos(theta(met_k))*cos(psi(met_k)+phi);
    M(3+met_k,1) = M(1,3+met_k);
    
    M(2,3+met_k) = MR*cos(theta(met_k))*sin(psi(met_k)+phi);
    M(3+met_k,2) = M(2,3+met_k);
    
    M(3,3+met_k) = MR*R*cos(theta(met_k));
    M(3+met_k,3) = M(3,3+met_k);
    M(3+met_k,3+met_k) = MR2;
end

%%  General Mass matrix
clear all
n =2;

theta = sym('theta',[n,1]);
dtheta = sym('dtheta',[n,1]);
psi = sym('psi',[n,1]);
gamma = sym('gamma',[n,1]);
R = sym('R',[n,1]);
H = sym(zeros(n,1));
conf = [R;gamma;psi;H];

syms phi dphi x dx y dy
syms ux uy uphi
Utheta = sym('Utheta',[n,1]);

%% model params
mdl.parameters = {'mp','mb','lp','lb','g','d','mP','mM','k','kt','c','ct',...
                    'epsilon','tau','theta_e','theta_s','JP33','JM33'};

%Reserve par index 11 for period time in AUTO
if length(mdl.parameters) >= 11
    mdl.parameters = {mdl.parameters{1:10} 'time' mdl.parameters{11:end}};
end

par = [];
for i=1:length(mdl.parameters)
    eval(strcat(mdl.parameters{i}, '=sym(''', mdl.parameters{i}, "');"));
    eval(sprintf(['par = [par;' mdl.parameters{i} '];']));
end

%% Conversion to paper masses
MT = mP + n*(mM+mb+mp);
MMT = mM+mb+mp;
MR = mb*lb-mp*lp;
MR2 = mb*lb^2+mp*lp^2;
JT = JP33+n*JM33;

%% Mass matrix
M = sym(zeros(3+n,3+n));
M(1,1) = MT;
M(2,2) = MT;
M(1,3) = -sum(MR.*sin(theta).*sin(psi+phi) + MMT.*R.*sin(gamma+phi));
M(3,1) = M(1,3);

M(2,3) = sum(MR*sin(theta).*cos(psi+phi) + MMT.*R.*cos(gamma+phi));
M(3,2) = M(2,3);

M(3,3) = sum(MMT.*R.^2 + 2.*MR.*R.*cos(psi-gamma).*sin(theta) + MR2.*sin(theta).^2) + JT;

for met_k = 1:n
    M(1,3+met_k) = MR*cos(theta(met_k))*cos(psi(met_k)+phi);
    M(3+met_k,1) = M(1,3+met_k);
    
    M(2,3+met_k) = MR*cos(theta(met_k))*sin(psi(met_k)+phi);
    M(3+met_k,2) = M(2,3+met_k);
    
    M(3,3+met_k) = MR.*R(met_k).*sin(psi(met_k)-gamma(met_k)).*cos(theta(met_k));
    M(3+met_k,3) = M(3,3+met_k);
    M(3+met_k,3+met_k) = MR2;
end

%% RHS
RHS = sym(zeros(3+n,1));
RHS(1) = -MMT*sum(dphi.^2.*R.*cos(gamma+phi))...
    - MR*sum(dtheta.^2.*sin(theta).*cos(psi+phi)...
    + 2.*dtheta.*dphi.*cos(theta).*sin(psi+phi)...
    + dphi.^2.*sin(theta).*cos(psi+phi))...
    + c*dx + k*x - ux;

RHS(2) = -MMT*sum(dphi.^2.*R.*sin(gamma+phi))...
    + MR*sum(-dtheta.^2.*sin(theta).*sin(psi+phi)...
    + 2.*dtheta.*dphi.*cos(theta).*cos(psi+phi)...
    - dphi.^2.*sin(theta).*sin(psi+phi))...
    + c*dy + k*y - uy;

RHS(3) = MR*sum(-dtheta.^2.*R.*sin(psi-gamma).*sin(theta)...
    + 2.*dphi.*dtheta.*R.*cos(psi-gamma).*cos(theta))...
    + 2*MR2.*sum(dphi.*dtheta.*sin(theta).*cos(theta))...
    + ct*dphi + kt*phi - uphi;

for met_k = 1:n
    theta_k = theta(met_k);
    gamma_k = gamma(met_k);
    psi_k = psi(met_k);
    R_k = R(met_k);
    RHS(3+met_k) = -MR*dphi^2*R_k*cos(psi_k-gamma_k)*cos(theta_k)...
        - MR2*dphi^2*sin(theta_k)*cos(theta_k)...
        + d*dtheta(met_k) - MR*g*sin(theta_k) - Utheta(met_k);
end

RHS = -RHS;



%% Check
q = [x;y;phi;theta;dx;dy;dphi;dtheta];
M_matlab = m3d_2_mass(0,q,par,conf);
M_matlab = M_matlab(3+n+1:end,3+n+1:end)
simplify(M_matlab-M)

RHS_matlab = m3d_2_row(0,q,par,conf,[ux;uy;uphi;Utheta]);
RHS_matlab = RHS_matlab(3+n+1:end)
simplify(RHS_matlab-RHS)
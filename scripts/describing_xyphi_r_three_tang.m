prepEnv
load('../data/tang_xyphi_3_HBEREF_Ham_1/sim_1')
%% Set some data
if(~exist('temp','dir'))
    mkdir('temp')
end
addpath('temp')

Nt = 3;
params = getMassParams(getparams(),Nt);
R = 0.01*ones(Nt,1);
R = 0.5*ones(Nt,1);
gamma = (0:(1/Nt):1-1e-6).'*2*pi;
psi = gamma-0.5*pi;

R = s.R;
gamma = s.gamma;
psi = s.psi;

% Generate transfer function and derivative function 
syms symomega

[G0s,Gs] = getTransfer_n_xyphi(Nt,symomega,1,params,R,gamma,psi);
Gds = diff(Gs,symomega);
G_func = matlabFunction(Gs,'file',['temp/G_func'],'vars',symomega,'outputs',{'G'});
dGdomega_func = matlabFunction(Gds,'file',['temp/dGdomega_func'],'vars',symomega,'outputs',{'dGdomega'});
%%
% Gradient provided optimization
fmuoptions = optimoptions(@fminunc,...
                          'TolFun',1e-4,'TolX',1e-4,...
                          'Display','iter-detailed','Algorithm','quasi-newton',...
                          'StepTolerance',1e-4,...
                          'MaxFunctionEvaluations',1500,...
                          'SpecifyObjectiveGradient',true,...
                          'DerivativeCheck','off');
                      
% Gradient provided optimization
fmuoptionsp = optimoptions(@fminunc,...
                          'TolFun',1e-12,'TolX',1e-12,...
                          'Display','iter-detailed','Algorithm','quasi-newton',...
                          'StepTolerance',1e-12,...
                          'MaxFunctionEvaluations',1500,...
                          'SpecifyObjectiveGradient',true,...
                          'DerivativeCheck','off');
                      
fmuoptionslsqnl = optimoptions('lsqnonlin',...
                          'TolFun',1e-12,'TolX',1e-12,...
                          'Display','iter-detailed','Algorithm','levenberg-marquardt',...
                          'StepTolerance',1e-15,...
                          'MaxFunctionEvaluations',1500,...
                          'SpecifyObjectiveGradient',true,...
                          'DerivativeCheck','off');
% Gradient provided optimization
fmuoptionst = optimoptions(@fminunc,...
                          'TolFun',1e-10,'TolX',1e-10,...
                          'Display','iter-detailed','Algorithm','trust-region',...
                          'StepTolerance',1e-10,...
                          'MaxFunctionEvaluations',1500,...
                          'SpecifyObjectiveGradient',true,...
                          'DerivativeCheck','off');       
                      
% Gradient provided optimization - quiet
fmuoptionsq = optimoptions(@fminunc,...
                          'TolFun',1e-10,'TolX',1e-10,...
                          'Display','none','Algorithm','quasi-newton',...
                          'StepTolerance',1e-10,...
                          'MaxFunctionEvaluations',1500,...
                          'SpecifyObjectiveGradient',true,...
                          'DerivativeCheck','off');

% Add numerical gradient check
fmuoptionsc = optimoptions(@fminunc,...
                          'TolFun',1e-10,'TolX',1e-10,...
                          'Display','iter-detailed','Algorithm','quasi-newton',...
                          'StepTolerance',1e-10,...
                          'MaxFunctionEvaluations',1500,...
                          'SpecifyObjectiveGradient',true,...
                          'FiniteDifferenceType','central',...
                          'FiniteDifferenceStepSize',1e-10,...
                          'DerivativeCheck','on');  

% Find numerical gradient
fmuoptionscn = optimoptions(@fminunc,...
                          'TolFun',1e-10,'TolX',1e-10,...
                          'Display','iter-detailed','Algorithm','quasi-newton',...
                          'StepTolerance',1e-10,...
                          'MaxFunctionEvaluations',5,...
                          'SpecifyObjectiveGradient',false,...
                          'DerivativeCheck','off');  
                      
% Numerical optimization
fmuoptionsn = optimoptions(@fminunc,...
                          'TolFun',1e-10,'TolX',1e-10,...
                          'Display','iter-detailed','Algorithm','quasi-newton',...
                          'StepTolerance',1e-10,...
                          'MaxFunctionEvaluations',1500,...
                          'SpecifyObjectiveGradient',false,...
                          'DerivativeCheck','off');

%% Optimization search with x phase
Nr = 3;
sF = s.results.analysis.fourier.Fzkn(1:6,1:Nr);
somega = s.results.analysis.fourier.omega;
sA = 2*abs(sF);
sP = (angle(sF)+0.5*pi);


% Perturbations
 sA(1:2,:) = 0;
 sA(3,3) = 0.003;

% sP(1,1) = sP(1,1)+0.1*pi;
% sP(5,1) = sP(5,1)+0.1*pi;
% sA(:,1) = sA(:,1)+0.0002;

if(Nr>2)
    sA(:,3) = sA(:,3)*1e2;
end
% Simplifications
sA(:,2) = 0;
sP(:,2) = 0;

% sA(3,1) = 0;
% sP(3,1) = 0;


rho_f = [reshape([sA;sP],[],1);somega];



A_ids = ones(6,Nr);
P_ids = ones(6,Nr);
A_ids(5:6,:) = 0;
P_ids(5:6,:) = 0;
A_ids(:,2) = 0;
P_ids(:,2) = 0;

P_ids(4,1) = 0;

O_ids = [1];

rho_ids = logical([reshape([A_ids;P_ids],[],1);O_ids]);

rho_i = rho_f(rho_ids);


fJ = @(rho) errDescnxyphi_r_three_tang(Nt,rho_f,Nr,G_func,dGdomega_func,params,s.U_gains,s.HStar,R,gamma,psi,eye(6*Nr),rho_ids,rho,1e4);
fJf = @(rho) errDescnxyphi_r_three_tang(Nt,rho_f,Nr,G_func,dGdomega_func,params,s.U_gains,s.HStar,R,gamma,psi,eye(6*Nr),rho_ids,rho,1e3);
fJlsq = @(rho) errDescnxyphi_r_lsq(Nt,rho_f,Nr,G_func,dGdomega_func,params,s.U_gains,s.HStar,R,gamma,psi,eye(6*Nr),rho_ids,rho,1e2);
fJ(rho_f(rho_ids))
%%
Amu = fminunc(fJ,rho_i,fmuoptionsp)
[J,gJ] = fJ(Amu)
%%
Amu = lsqnonlin(fJlsq,Amu,[],[],fmuoptionslsqnl)
[J,gJ] = fJ(Amu)
%%
Amu = fminsearch(fJ,Amu,optimset('Display','iter','TolFun',1e-15,'TolX',1e-15,'MaxFunEvals',5e6,'MaxIter',2e6,'PlotFcns',@optimplotfval))
[J,gJ] = fJ(Amu)
Amu = fminsearch(fJf,Amu,optimset('Display','iter','TolFun',1e-15,'TolX',1e-15,'MaxFunEvals',5e6,'MaxIter',2e6,'PlotFcns',@optimplotfval))
[J,gJ] = fJ(Amu)
%%
Amu = fminunc(fJ,Amu,fmuoptions)
[J,gJ] = fJ(Amu)
%%
Amu = fminunc(fJ,rho_i,fmuoptionsn)
[J,gJ] = fJ(Amu)

%%
Amu = fminunc(fJ,Amu,fmuoptionsn)
[J,gJ] = fJ(Amu)
%%
Amu = fminsearch(fJ,Amu,optimset('Display','iter','TolFun',1e-14,'TolX',1e-10,'MaxFunEvals',5e4,'MaxIter',2e4,'PlotFcns',@optimplotfval))
[J,gJ] = fJ(Amu)
%%
fJf = @(rho) errDescnxyphi_r(Nt,rho_f,Nr,G_func,dGdomega_func,params,s.U_gains,s.HStar,R,gamma,psi,eye(6*Nr),rho_ids,rho,1e2);


Amu = fminunc(fJf,Amu,fmuoptionsc)
[J,gJ] = fJ(Amu)
%%
Ax = 0.58;
Ay = 0.58;
Aphi = 0.03;
At1 = 0.7;
At2 = 0.7;
At3 = 0.7;
Px = 0.5*pi*1e2;
Py = pi*1e2;
Pphi = 1;
Pt1 = 0;
Pt2 = 2*1e2;
Pt3 = 4*1e2;

A = [Ax;Ay;Aphi;At1;At2;At3];
P = [Px;Py;Pphi;Pt1;Pt2;Pt3];
O = 10.9;

rho_f = [A;P;1e-3*A;P;O];



A_ids = [1 1 1 1 1 1].';
P_ids = [1 1 1 0 1 1].';
O_ids = [1];
rho_ids = logical([A_ids;P_ids;A_ids;A_ids;O_ids]);

rho_i = rho_f(rho_ids);

fJ = @(rho) errDescnxyphi_r(Nt,rho_f,Nr,G_func,dGdomega_func,params,s.U_gains,s.HStar,R,gamma,psi,eye(6),rho_ids,rho,1e3);
fJf = @(rho) errDescnxyphi_r(Nt,rho_f,Nr,G_func,dGdomega_func,params,s.U_gains,s.HStar,R,gamma,psi,eye(6),rho_ids,rho,1e5);
Amu = fminunc(fJ,rho_i,fmuoptions)
[J,gJ] = fJ(Amu)
Amu = fminunc(fJf,Amu,fmuoptionsc)
[J,gJ] = fJ(Amu)

prepEnv
load('../data/tang_xyphi_3_HBEREF_Ham_1/sim_1')
%% Set some data
if(~exist('temp','dir'))
    mkdir('temp')
end
addpath('temp')

Nt = 2;
params = getMassParams(getparams(),Nt);
R = 0.2*ones(Nt,1);
gamma = zeros(Nt,1);
psi = zeros(Nt,1);

%R = s.R;
%gamma = s.gamma;
%psi = s.psi;

% Generate transfer function and derivative function 
syms symomega

[G0s,Gs] = getTransfer_n_xyphi(Nt,symomega,1,params,R,gamma,psi);
Gds = diff(Gs,symomega);
G_func = matlabFunction(Gs,'file',['temp/G_func'],'vars',symomega,'outputs',{'G'});
dGdomega_func = matlabFunction(Gds,'file',['temp/dGdomega_func'],'vars',symomega,'outputs',{'dGdomega'});
%%
% Gradient provided optimization
fmuoptions = optimoptions(@fminunc,...
                          'TolFun',1e-10,'TolX',1e-10,...
                          'Display','iter-detailed','Algorithm','quasi-newton',...
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
sF = s.results.analysis.fourier.Fzkn(1:6,1);
somega = s.results.analysis.fourier.omega;
sA = 2*abs(sF(1:5));
sP = angle(sF(1:5))+0.5*pi;
sA(1:2) = sA(1:2)*1e3;
sA(3) = sA(3)*1e2;
sA(2) =0;
sA(3) = 0;


Ax = 0.8;
Ay = 0;
Aphi =0;
At1 = 0.7;
At2 = 0.7;
Px = 1;
Py = 0;
Pphi = 0;
Pt1 = 0;
Pt2 = 1;

rho_f = [Ax; Ay; Aphi; At1; At2;Px; Py; Pphi; Pt1; Pt2;somega];



A_ids = [1 0 0 1 1];
P_ids = [1 0 0 0 1];
O_ids = [1];

rho_ids = logical([A_ids,P_ids,O_ids]);

rho_i = rho_f(rho_ids);

fJ = @(rho) errDescnxyphi(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains(1:2),s.HStar(1:2),R,gamma,psi,eye(5),rho_ids,rho,1e3);
fJf = @(rho) errDescnxyphi(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains(1:2),s.HStar(1:2),R,gamma,psi,eye(5),rho_ids,rho,1e5);
Amu = fminunc(fJ,rho_i,fmuoptions)
[J,gJ] = fJ(Amu)
Amu = fminunc(fJf,Amu,fmuoptionsc)
[J,gJ] = fJ(Amu)
%%
Ax = 0.8;
Ay = 0.8;
Aphi = 0.0;
At1 = 0.7;
At2 = 0.7;
At3 = 0.7;
Px = 1;
Py = 0.5;
Pphi = 1;
Pt1 = 0;
Pt2 = 2;
Pt3 = 4;

A = [Ax;Ay;Aphi;At1;At2;At3];
P = [Px;Py;Pphi;Pt1;Pt2;Pt3];
O = 10.9;

rho_f = [A;P;O];



A_ids = [1 1 0 1 1 1];
P_ids = [1 1 0 0 1 1];
O_ids = [1];
rho_ids = logical([A_ids,P_ids,O_ids]);

rho_i = rho_f(rho_ids);

fJ = @(rho) errDescnxyphi(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,R,gamma,psi,eye(6),rho_ids,rho,1e3);
fJf = @(rho) errDescnxyphi(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,R,gamma,psi,eye(6),rho_ids,rho,1e5);
Amu = fminunc(fJ,rho_i,fmuoptionsn)
[J,gJ] = fJ(Amu)
Amu = fminunc(fJf,Amu,fmuoptionsc)
[J,gJ] = fJ(Amu)
%% Optimization search with theta2 phase
Ax = 0.8;
At1 = 0.7;
At2 = 0.7;
Px = sP1(1);
Pt1 = 0;
Pt2 = 1;

A = [Ax;At1;At2];
P = [Px;Pt1;Pt2];
O = somega;

rho_f = [A;P;O];


A_ids = [1 1 1];
P_ids = [0 0 1];
O_ids = [0];
rho_ids = logical([A_ids,P_ids,O_ids]);

rho_i = rho_f(rho_ids);

fJ = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e3);
fJf = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e5);
Amu = fminunc(fJ,rho_i,fmuoptions)
[J,gJ] = fJ(Amu)
Amu = fminunc(fJf,Amu,fmuoptions)
[J,gJ] = fJ(Amu)
%% Optimization search with all phase

Ax = 0.8;
At1 = 0.7;
At2 = 0.7;
Px = 1;
Pt1 = 0;
Pt2 = 1;

A = [Ax;At1;At2];
P = [Px;Pt1;Pt2];
O = somega;

rho_f = [A;P;O];


A_ids = [1 1 1];
P_ids = [1 0 1];
O_ids = [0];
rho_ids = logical([A_ids,P_ids,O_ids]);

rho_i = rho_f(rho_ids);

fJ = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e3);
fJf = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e5);
Amu = fminunc(fJ,rho_i,fmuoptions)
[J,gJ] = fJ(Amu)
Amu = fminunc(fJf,Amu,fmuoptions)
[J,gJ] = fJ(Amu)

%% Optimization search with x phase - anti

Ax = 0.001;
At1 = 0.7;
At2 = 0.7;
Px = 1;
Pt1 = 0;
Pt2 = pi*1e2;

A = [Ax;At1;At2];
P = [Px;Pt1;Pt2];
O = somega;

rho_f = [A;P;O];


A_ids = [1 1 1];
P_ids = [0 0 0];
O_ids = [0];
rho_ids = logical([A_ids,P_ids,O_ids]);

rho_i = rho_f(rho_ids);

fJ = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e3);
fJf = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e5);
Amu = fminunc(fJ,rho_i,fmuoptionsc)
[J,gJ] = fJ(Amu)
Amu = fminunc(fJf,Amu,fmuoptions)
[J,gJ] = fJ(Amu)
%% Optimization search with theta2 phase - anti
Ax = 0;
At1 = 0.7;
At2 = 0.7;
Px = sP1(1);
Pt1 = 0;
Pt2 = 3.2e2;

A = [Ax;At1;At2];
P = [Px;Pt1;Pt2];
O = somega;

rho_f = [A;P;O];


A_ids = [0 1 1];
P_ids = [0 0 1];
O_ids = [0];
rho_ids = logical([A_ids,P_ids,O_ids]);

rho_i = rho_f(rho_ids);

fJ = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e3);
fJf = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e5);
Amu = fminunc(fJ,rho_i,fmuoptions)
[J,gJ] = fJ(Amu)
Amu = fminunc(fJf,Amu,fmuoptions)
[J,gJ] = fJ(Amu)
%% Optimization search with all phase - anti

Ax = 0.001;
At1 = 0.7;
At2 = 0.7;
Px = 1;
Pt1 = 0;
Pt2 = 3.2e2;

A = [Ax;At1;At2];
P = [Px;Pt1;Pt2];
O = somega;

rho_f = [A;P;O];


A_ids = [1 1 1];
P_ids = [1 0 1];
O_ids = [0];
rho_ids = logical([A_ids,P_ids,O_ids]);

rho_i = rho_f(rho_ids);

fJ = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e3);
fJf = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e5);
Amu = fminunc(fJ,rho_i,fmuoptions)
[J,gJ] = fJ(Amu)
Amu = fminunc(fJf,Amu,fmuoptions)
[J,gJ] = fJ(Amu)


%% Optimization search with all phase and omega

Ax = 0.8;
At1 = 0.7;
At2 = 0.7;
Px = 1;
Pt1 = 0;
Pt2 = 1;

A = [Ax;At1;At2];
P = [Px;Pt1;Pt2];
O = 10.95;

rho_f = [A;P;O];


A_ids = [1 1 1];
P_ids = [1 0 1];
O_ids = [1];
rho_ids = logical([A_ids,P_ids,O_ids]);

rho_i = rho_f(rho_ids);

fJ = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e3);
fJf = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e5);
Amu = fminunc(fJ,rho_i,fmuoptions)
[J,gJ] = fJ(Amu)
Amu = fminunc(fJf,Amu,fmuoptions)
[J,gJ] = fJ(Amu)

%% Optimization search with all phase and omega - anti

Ax = 0.0;
At1 = 0.7;
At2 = 0.7;
Px = 1;
Pt1 = 0;
Pt2 = 3.2e2;

A = [Ax;At1;At2];
P = [Px;Pt1;Pt2];
O = 10.75;

rho_f = [A;P;O];


A_ids = [0 1 1];
P_ids = [0 0 1];
O_ids = [1];
rho_ids = logical([A_ids,P_ids,O_ids]);

rho_i = rho_f(rho_ids);

fJ = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e3);
fJf = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e5);
Amu = fminunc(fJ,rho_i,fmuoptions)
[J,gJ] = fJ(Amu)
Amu = fminunc(fJf,Amu,fmuoptions)
[J,gJ] = fJ(Amu)
%% Optimization search with omega

Ax = 0.8;
At1 = 0.7;
At2 = 0.7;
Px = 1e2;
Pt1 = 0;
Pt2 = 1e2;

omv = (-0.2:0.005:0.2) + somega;
Jm = zeros(size(omv));

oAt1ma = Jm;
oAt2ma = Jm;
oPt2ma = Jm;
oJma = Jm;

oAxmi = Jm;
oAt1mi = Jm;
oAt2mi = Jm;
oPxmi = Jm;
oPt2mi = Jm;
oJmi = Jm;

oAxmia = Jm;
oAt1mia = Jm;
oAt2mia = Jm;
oPxmia = Jm;
oPt2mia = Jm;
oJmia = Jm;

N = numel(Jm)
%%

n= 0;
for omidx = 1:N
    
    omidx
    omega = omv(omidx);
    A = [0;At1;At2];
    P = [Px;Pt1;3.2e2];

    A_ids = [0 1 1];
    P_ids = [0 0 1];
    O_ids = [0];
    rho_ids = logical([A_ids,P_ids,O_ids]);
    
    rho_f = [A;P;omega];
    rho_i = rho_f(rho_ids);
    
    fJ = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e3);
    fJf = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e5);
    Amu = fminunc(fJ,rho_i,fmuoptionsq);
    [J,gJ] = fJ(Amu);
    Amu = fminunc(fJf,Amu,fmuoptionsq);
    [J,gJ] = fJ(Amu);
    
    oAt1ma(omidx) = Amu(1);
    oAt2ma(omidx) = Amu(2);
    oPt2ma(omidx) = Amu(3)*1e-2;
    oJma(omidx) = fJ(Amu);
    
    
    A = [Ax;At1;At2];
    P = [Px;Pt1;Pt2];

    A_ids = [1 1 1];
    P_ids = [1 0 1];
    O_ids = [0];
    rho_ids = logical([A_ids,P_ids,O_ids]);
    
    rho_f = [A;P;omega];
    rho_i = rho_f(rho_ids);
    
    fJ = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e3);
    fJf = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e5);
    Amu = fminunc(fJ,rho_i,fmuoptionsq);
    [J,gJ] = fJ(Amu);
    Amu = fminunc(fJf,Amu,fmuoptionsq);
    [J,gJ] = fJ(Amu);
    
    oAxmi(omidx) = Amu(1);
    oAt1mi(omidx) = Amu(2);
    oAt2mi(omidx) = Amu(3);
    oPxmi(omidx) = Amu(4)*1e-2;
    oPt2mi(omidx) = Amu(5)*1e-2;
    oJmi(omidx) = fJ(Amu);
    
%     A = [Ax*1e-2;At1;At2];
%     P = [Px;Pt1;3.2e2];
% 
%     A_ids = [1 1 1];
%     P_ids = [1 0 1];
%     O_ids = [0];
%     rho_ids = logical([A_ids,P_ids,O_ids]);
%     
%     rho_f = [A;P;omega];
%     rho_i = rho_f(rho_ids);
%     
%     fJ = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e3);
%     fJf = @(rho) errDescnx(Nt,rho_f,G_func,dGdomega_func,params,s.U_gains,s.HStar,diag([1e3 1 1]),rho_ids,rho,1e5);
%     Amu = fminunc(fJ,rho_i,fmuoptions)
%     [J,gJ] = fJ(Amu)
%     Amu = fminunc(fJf,Amu,fmuoptions)
%     [J,gJ] = fJ(Amu)
%     
%     oAxmia(omidx) = Amu(1);
%     oAt1mia(omidx) = Amu(2);
%     oAt2mia(omidx) = Amu(3);
%     oPxmia(omidx) = Amu(4);
%     oPt2mia(omidx) = Amu(5);
%     oJmia(omidx) = fJ(Amu);
    
end


%% Anti-phase
figure;

subplot(4,1,1)
plot(omv,oAt1ma,'kx')
ylabel('$A_{1\theta_1}$','Interpreter','Latex')
legend('hide')

subplot(4,1,2)
plot(omv,oAt2ma,'kx')
ylabel('$A_{1\theta_2}$','Interpreter','Latex')
legend('hide')


subplot(4,1,3)
plot(omv,abs(oPt2ma),'kx')
ylabel('$\psi_{1\theta_2}$','Interpreter','Latex')
legend('hide')

subplot(4,1,4)
plot(omv,log10(oJma),'kx')
xlabel('$\omega$')
ylabel('$\log_{10}(J_o)$','Interpreter','Latex')
legend('hide')
matlab2tikz('2x_hbe_optim_omega_anti.tex','parseStrings',false,...
                                    'height','\figureheight',...
                                    'width','\figurewidth',...
                                    'showInfo', false);
%% In-phase
figure;
subplot(6,1,1)
plot(omv,oAxmi,'kx')
ylabel('$A_{1x}$')
legend('hide')

subplot(6,1,2)
plot(omv,oAt1mi,'kx')
ylabel('$A_{1\theta_1}$','Interpreter','Latex')
legend('hide')

subplot(6,1,3)
plot(omv,oAt2mi,'kx')
ylabel('$A_{1\theta_2}$','Interpreter','Latex')
legend('hide')

subplot(6,1,4)
plot(omv,oPxmi,'kx')
ylabel('$\psi_{1x}$','Interpreter','Latex')
legend('hide')
subplot(6,1,5)
plot(omv,zeros(size(omv)),'r')
hold on
plot(omv,pi*ones(size(omv)),'r')
plot(omv,oPt2mi,'kx')
ylabel('$\psi_{1\theta_2}$','Interpreter','Latex')
legend('hide')

subplot(6,1,6)
plot(omv,log10(oJmi),'kx')
xlabel('$\omega$')
ylabel('$\log_{10}(J_o)$')
legend('hide')
matlab2tikz('2x_hbe_optim_omega_in.tex','parseStrings',false,...
                                    'height','\figureheight',...
                                    'width','\figurewidth',...
                                    'showInfo', false);
%% Anti-phase
figure;
subplot(6,1,1)
plot(omv,oAxmia,'kx')
ylabel('$A_{1x}$')
legend('hide')

subplot(6,1,2)
plot(omv,oAt1mia,'kx')
ylabel('$A_{1\theta_1}$')
legend('hide')

subplot(6,1,3)
plot(omv,oAt2mia,'kx')
ylabel('$A_{1\theta_2}$')
legend('hide')

subplot(6,1,4)
plot(omv,oPxmia,'kx')
ylabel('$\psi_{1x}$')
legend('hide')
subplot(6,1,5)
plot(omv,oPt2mia,'kx')
ylabel('$\psi_{1\theta_2}$')
legend('hide')

subplot(6,1,6)
plot(omv,log10(oJmia),'kx')
xlabel('$\omega$')
ylabel('$\log_{10}(J_o)$')
legend('hide')
matlab2tikz('2x_hbe_optim_omega_anti_full.tex','parseStrings',false,...
                                    'height','\figureheight',...
                                    'width','\figurewidth',...
                                    'showInfo', false);

%
%%
fmuoptions = optimoptions(@fminunc,'TolFun',1e-8,'TolX',1e-8,'Display','iter','Algorithm','quasi-newton','StepTolerance',1e-8,'MaxFunctionEvaluations',2500);
rD = 1;

a0 = [0;0;0];

Axi = 1e-3;
At1i = 0.72;
At2i = 0.72;
Pxi = 0.01;
Pt1 = 0;


fak = @(A,P) convertAmpPhaseToComplex(A,P);


params = getMassParams(getparams(),2);


fJ = @(A) errDescWithTransfer_2_x(a0,fak([A(1:3)],[A(4); 0; A(5)]),fak([A(1:3)],[A(4); 0; A(5)]),rD,A(6),params,s.U_gains,s.HStar,diag([1e3 1 1]));
Amu = fminunc(fJ,[Axi;At1i;At2i;0.001;0.001;10.95],fmuoptions)
fJ(Amu)



    

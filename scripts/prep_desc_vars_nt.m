prepEnv

%% Setup some data
if(~exist('temp','dir'))
    mkdir('temp')
end
addpath('temp')

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
%%


params = getMassParams(getparams(),Nt);
R = ring_diam*ones(Nt,1);
gamma = (0:(1/Nt):1-1e-6).'*2*pi;
if(ring_orientation=="perp")
    psi = gamma;
elseif(ring_orientation=="tang")
    psi = gamma-0.5*pi;
else
    error('orientation?');
end

U_gains = met_U_gains*ones(6,1);
HStar = met_HStar*ones(6,1);

% Generate transfer function and derivative function 
syms symomega

[G0s,Gs] = getTransfer_n_xyphi(Nt,symomega,1,params,R,gamma,psi);
Gds = diff(Gs,symomega);
G_func = matlabFunction(Gs,'file',['temp/G_func'],'vars',symomega,'outputs',{'G'});
dGdomega_func = matlabFunction(Gds,'file',['temp/dGdomega_func'],'vars',symomega,'outputs',{'dGdomega'});
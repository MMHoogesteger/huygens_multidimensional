function mdl = derive_model(mdl,derive_matlab)

%derive_model - Derives model for a certain amount of metronomes
%Calculates all energies and the equations of motion for a certain model
%
% Syntax:  derive_model
%
% Preconditions
%   mdl - struct with at least the following fields set:
%
%       mdl.N - The number of metronomes
%
%       Logicals: Does the model contain ??? as degree of freedom?
%       mdl.DOF_x
%       mdl.DOF_y
%       mdl.DOF_phi
%       mdl.PSI_zero    Are all metronomes swinging in the x-direction,
%                       then PSI = zero
%       mdl.GAMMA_zero  Are all metronomes standing on the x-axis, then 
%                       GAMMA = zero
%
%       mdl.modelName   String with the name of the model, used for file
%                       outputs and so on. Must be unique.
%
%   derive_matlab - determines wether manually entered expressions are
%                   used, or matlab derives the Euler-Lagrange equations of
%                   motion itself. 0 = manual expressions, 1 = automatic
%                   derivation, >1 = both and checks.
%
% Other m-files required: f_derive_model, escapement_ham, escapement_or
% Subfunctions: none
% MAT-files required: none
%
% See also: derive_models, f_derive_model

% Author: M.M. Hoogesteger - Msc student of Mechanical Engineering
% Eindhoven University of Technology, Mechanical Engineering, Dynamics and
% control
% email address: m.m.hoogesteger@student.tue.nl
% December 2017; Last revision: 3-May-2018



%% Setup model - Prepare symbolic variables
% Time
syms t
psi = 0;
gamma = 0;
% Parameters - strings
mdl.configuration = {'R','gamma','psi'};
mdl.parameters = {'mp','mb','lp','lb','g','d','mP','mM','k','kt','c','ct',...
                    'epsilon','tau','theta_e','theta_s','JP33','JM33'};
if length(mdl.parameters) >= 11
    mdl.parameters = {mdl.parameters{1:10} 'time' mdl.parameters{11:end}};
end

% Nonzero DOFs
mdl.DOFs = find([mdl.DOF_x mdl.DOF_y mdl.DOF_phi]);
mdl.NqP = length(mdl.DOFs);
mdl.Nq = mdl.NqP+mdl.N;

% Primary DOFs
q_arrays = {'x';'y';'phi'};
q_array=q_arrays(mdl.DOFs);

% Metronome DOFs
for i=1:mdl.N
    q_array{mdl.NqP+i,1} = ['theta' num2str(i)];
end


mdl.q_array = q_array;

%% Syms
% Physical configuration - Generate one symbolic variable of every 
%                          configuration parameter for every metronome

conf = [];
for i=1:length(mdl.configuration)
    eval(strcat(mdl.configuration{i}, '=sym(''', mdl.configuration{i}, "',[mdl.N 1]);"));
    eval(sprintf(['conf = [conf;' mdl.configuration{i} '];']));
end
assume(conf,'real');

% Parameters - Generate one symbolic variable of every parameter
par = [];
for i=1:length(mdl.parameters)
    eval(strcat(mdl.parameters{i}, '=sym(''', mdl.parameters{i}, "');"));
    eval(sprintf(['par = [par;' mdl.parameters{i} '];']));
end
assume(par,'real');
assume(par,'positive');

% Generate two rotational inertia matrices and set the only parameters used
% in the 3,3 place
JP = sym('JP',3);
JP(3,3) = JP33;
JM = sym('JM',3);
JM(3,3) = JM33;


% Genertate input symbolics
syms Ux Uy Uphi;
Utheta = [];
for i=1:mdl.N
    eval(strcat(['Utheta' num2str(i)], "=sym('",['Utheta' num2str(i)], "');"));
    eval(sprintf(['Utheta = [Utheta;Utheta' num2str(i) '];']));
end

% Generalized coordinates and inputs assembly
q = [];
dq = [];
ddq = [];
U = [];
for i=1:length(mdl.q_array)
    eval(strcat(mdl.q_array{i}, '=sym(''', mdl.q_array{i}, "');"));
    eval(sprintf(['q = [q;' mdl.q_array{i} '];']));
    eval(strcat('d', mdl.q_array{i}, '=sym(''d', mdl.q_array{i}, "');"));
    eval(sprintf(['dq = [dq;d' mdl.q_array{i} '];']));
    eval(strcat('dd', mdl.q_array{i}, '=sym(''dd', mdl.q_array{i}, "');"));
    eval(sprintf(['ddq = [ddq;dd' mdl.q_array{i} '];']));
    eval(sprintf(['U = [U;U' mdl.q_array{i} '];']));
end

z = [q;dq];

%% Other inputs
U_gains = sym('U_gains',[mdl.N,1]);
HStar = sym('HStar',[mdl.N,1]);

theta = q(mdl.NqP+1:end);
dtheta = dq(mdl.NqP+1:end);

s.q = q;
s.dq = dq;
s.ddq = ddq;
s.U = U;
s.z = z;
s.t = t;
s.U_gains = U_gains;
s.HStar = HStar;
mdl.symbolics = s;

%% Constraints
% These constraints don't change q and so on, but do change the expressions
% that do not reference q(1) but x for example. Therefore q always has the
% general extensive form, but the expressions and functions become simpler.

if(~mdl.DOF_x)
    x=0;
    dx=0;
end

if(~mdl.DOF_y)
    y=0;
    dy=0;
end

if(~mdl.DOF_phi)
    phi=0;
    dphi=0;
end

if(mdl.PSI_zero)
    psi=psi*0;
end

if(mdl.GAMMA_zero)
    gamma=gamma*0;
end

%% Find symbolic expressions
B = cos(gamma).*R;
A = sin(gamma).*R;

% Derive using Matlab (slow)
if(derive_matlab>0)
    
    [MM, RHS, T, V, H, HM] = f_derive_model(mdl, U, mP, mM, mb, mp, lb, lp, JP, JM, A, B, gamma, q, dq, ddq, x, dx, y, dy, phi, dphi, theta, dtheta, psi, c, ct, k, kt, d, g);
end

% Use handmade expressions (faster)
if(derive_matlab~=1)
    %% Paper masses
    MT = mP + mdl.N*(mM+mb+mp);
    MMT = mM+mb+mp;
    MR = mb*lb-mp*lp;
    MR2 = mb*lb^2+mp*lp^2;
    JT = JP33+mdl.N*JM33;

    %% Functions
    % Mass matrix
    MMh = mass_matrix_hand(mdl, MT, MMT, MR, MR2, JT, R, gamma, phi, theta, psi);

    % RHS
    RHSh = right_side_hand(mdl, Ux, Uy, Uphi, Utheta, MMT, MR, MR2, R, gamma, x, dx, y, dy, phi, dphi, theta, dtheta, psi, c, ct, k, kt, d, g);

    
    % Potential energy
    Vh = potential_hand(x, y, phi, theta, k, kt, g, MR);
    % Kinetic energy
    Th = kinetic_hand(MT, MMT, MR, MR2, JT, A, B, dx, dy, phi, dphi, theta, dtheta, psi);
    
    % Hamiltonian and Energy of the various parts
    Hh = Th+Vh;
    HMh = hamiltonian_parts_hand(mdl,x, dx, y, dy, phi, dphi, theta, dtheta, mP, k, JP33, kt, MR2, MR, g);
    
    
    
end

%% If not derive matlab and both, set values to manual and check
if(derive_matlab==0)
    MM = MMh;
    RHS = RHSh;
    V = Vh;
    T = Th;
    H = Hh;
    HM = HMh;
end


if(derive_matlab>1)
    c1 = simplify(MM-MMh);
    c2 = simplify(RHS-RHSh);
    c3 = simplify(T-Th);
    c4 = simplify(V-Vh);
    c5 = simplify(H-Hh);
    c6 = simplify(HM-HMh);
    c = [c1(:);c2(:);c3;c4;c5;c6(:)];
    
    if(~isequal(c,sym(zeros(size(c)))))
        error('Hand and Matlab not equal')
    end
end

%% Escapements
U_OR = escapement_or(mdl,theta, dtheta, epsilon, theta_s, theta_e, tau);
Htheta = HM(mdl.NqP+1:end);
U_HAM = escapement_ham(mdl,dtheta,Htheta,HStar,U_gains);




%% Conversion to 1-order ODE:


MM_11 = eye(mdl.Nq);
MM_12 = zeros(mdl.Nq);
MM_21 = MM_12;
MM_22 = MM;

MM = [MM_11 MM_12; MM_21 MM_22];

RHS = [dq;RHS];

%% Writing
% Write equations to matlab file
mdl.fMM = matlabFunction(MM,'file',['models\' mdl.modelName '_mass'],'vars',{t,z,par,conf},'outputs',{'M'});
mdl.fRHS = matlabFunction(RHS,'file',['models\' mdl.modelName '_row'],'vars',{t,z,par,conf,U},'outputs',{'dxdt'});
mdl.fU_OR = matlabFunction(U_OR,'file',['models\' mdl.modelName '_esc_or'],'vars',{t,z,par,conf},'outputs',{'u'});
mdl.fU_HAM = matlabFunction(U_HAM,'file',['models\' mdl.modelName '_esc_ham'],'vars',{t,z,par,conf,U_gains,HStar},'outputs',{'u'});
mdl.fH = matlabFunction(H,'file',['models\' mdl.modelName '_ham'],'vars',{t,z,par,conf},'outputs',{'H'});
mdl.fHM = matlabFunction(HM,'file',['models\' mdl.modelName '_ham_m'],'vars',{t,z,par,conf},'outputs',{'HM'});

save(['models\' mdl.modelName], 'mdl')
end

function [MM, RHS, T, V, H, HM] = f_derive_model(mdl, U, mP, mM, mb, mp, lb, lp, JP, JM, A, B, gamma, q, dq, ddq, x, dx, y, dy, phi, dphi, theta, dtheta, psi, c, ct, k, kt, d, g)

%% Position vectors
e0 = eye(3);
AP0 = get_rotation_matrix(3,phi);
eP = AP0*e0;
rP = [x y 0];

for i=1:mdl.N
    % Rotation matrices
    AMiP{i} = get_rotation_matrix(3,psi(i));
    AMi0{i} = simplify(AMiP{i}*AP0);
    AiMi{i} = get_rotation_matrix(2,theta(i));
    Ai0{i} = simplify(AiMi{i}*AMi0{i});

    % Metronome coordinate system
    ei{i} = Ai0{i}*e0;

    % Metronome body position vector
    rMiP_P = [B(i) A(i) 0];
    rMiP{i} = rMiP_P*eP;
    rMi{i} = rMiP{i}+rP;

    % Bob position vector
    rbiMi_i{i} = [0 0 lb];
    rbiMi{i} = rbiMi_i{i}*ei{i};
    rbi{i} = rbiMi{i} + rMi{i};

    % Pendulum position vector
    rpiMi_i{i} = [0 0 -lp];
    rpiMi{i} = rpiMi_i{i}*ei{i};
    rpi{i} = rpiMi{i} + rMi{i};
end

%% Velocities

drP = (jacobian(rP,q)*dq).';
for i=1:mdl.N
    drMi{i} = (jacobian(rMi{i} ,q)*dq).';
    drbi{i} = (jacobian(rbi{i} ,q)*dq).';
    drpi{i} = (jacobian(rpi{i} ,q)*dq).';
end

omegaP0 = [0 0 dphi]*e0;


%% Kinetic energy

% Platform rotation
rTP = 0.5*omegaP0*AP0.'*JP*AP0*omegaP0.';

% Rotation of one metronome (simplified, all the same)
rTM = 0.5*omegaP0*AP0.'*JM*AP0*omegaP0.';

% Platform translation
tTP = simplify(0.5*mP*(drP*drP.'));

tTPx = 0.5.*mP.*dx.^2;
tTPy = 0.5.*mP.*dy.^2;

if(~isequal(simplify(tTPy+tTPx-tTP),sym(zeros(1))))
    error('Some error with kinetic energy of platform translation');
end

for i=1:mdl.N
    % Metronome translation
    tTMi{i} = simplify(0.5*mM*drMi{i}*drMi{i}.');

    % Bob translation
    tTbi{i} = simplify(0.5*mb*drbi{i}*drbi{i}.');

    % Pendulum translation
    tTpi{i} = simplify(0.5*mp*drpi{i}*drpi{i}.');
end

% Total kinetic
TT = rTP+tTP;
for i=1:mdl.N
    TT = TT + rTM + tTMi{i} + tTbi{i} + tTpi{i};
end
T = simplify(expand(TT));

%% Potential energy
g_vec = [0 0 -g];
for i=1:mdl.N
    Vbi{i} = -mb*g_vec*rbi{i}.';
    Vpi{i} = -mp*g_vec*rpi{i}.';
end

Vkx = 0.5*k*x^2;
Vky = 0.5*k*y^2;
Vkphi = 0.5*kt*phi^2;

VT = Vkx+Vky+Vkphi;

for i=1:mdl.N
    VT = VT + Vbi{i} + Vpi{i};
end
V = simplify(expand(VT));

%% Lagrangian / Hamiltonian

H = simplify(T+V);
L = simplify(T-V);


for i=1:mdl.N
    Htheta(i,1) = tTbi{i} + tTpi{i} + Vbi{i} + Vpi{i};
end
% Set all dependency on platform to zero fro metronomes.
Htheta = subs(Htheta,[x,dx,y,dy,phi,dphi],[0 0 0 0 0 0]);

HM = [Vkx+tTPx;Vky+tTPy;Vkphi+rTP;Htheta];
HM = HM([mdl.DOFs 3+(1:mdl.N)]);


%% Non-conservative generalized forces
Qncs = [-c*dx;-c*dy;-ct*dphi];
Qnc = [Qncs(mdl.DOFs);-d*dtheta];


%% Equations of motion
dLdq = simplify(jacobian(L,q).');

dLd_dq = simplify(jacobian(L,dq).');
MM = simplify(jacobian(dLd_dq,dq));
CMatrix = simplify(jacobian(dLd_dq,q)*dq);
ddtdLd_dq = simplify(MM*ddq+CMatrix);


% dd_dqdLd_dq*ddq + ddqdLd_dq*dq -  dLdq- Qnc = 0;
Eq = simplify(ddtdLd_dq - dLdq-Qnc-U);



% dd_dqdLd_dq*ddq = -ddqdLd_dq*dq + dLdq + Qnc
RHS = simplify(-CMatrix+dLdq+Qnc+U);
end
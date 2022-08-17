function u = mphi_2_esc_or(t,in2,in3,in4)
%MPHI_2_ESC_OR
%    U = MPHI_2_ESC_OR(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    09-Aug-2018 13:59:47

dtheta1 = in2(5,:);
dtheta2 = in2(6,:);
epsilon = in3(14,:);
tau = in3(15,:);
theta1 = in2(2,:);
theta2 = in2(3,:);
theta_e = in3(16,:);
theta_s = in3(17,:);
t2 = 1.0./epsilon;
t3 = dtheta1.*t2;
t4 = tanh(t3);
t5 = dtheta2.*t2;
t6 = tanh(t5);
u = [0.0;t4.*tau.*(tanh(t2.*(theta_e-t4.*theta1))-tanh(t2.*(theta_s-t4.*theta1))).*(1.0./2.0);t6.*tau.*(tanh(t2.*(theta_e-t6.*theta2))-tanh(t2.*(theta_s-t6.*theta2))).*(1.0./2.0)];

function u = mxyphi_1_esc_or(t,in2,in3,in4)
%MXYPHI_1_ESC_OR
%    U = MXYPHI_1_ESC_OR(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Sep-2018 09:48:36

dtheta1 = in2(8,:);
epsilon = in3(14,:);
tau = in3(15,:);
theta1 = in2(4,:);
theta_e = in3(16,:);
theta_s = in3(17,:);
t2 = 1.0./epsilon;
t3 = dtheta1.*t2;
t4 = tanh(t3);
u = [0.0;0.0;0.0;t4.*tau.*(tanh(t2.*(theta_e-t4.*theta1))-tanh(t2.*(theta_s-t4.*theta1))).*(1.0./2.0)];

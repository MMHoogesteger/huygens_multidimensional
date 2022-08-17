function u = mxyphi_1_esc_ham(t,in2,in3,in4,U_gains1,HStar1)
%MXYPHI_1_ESC_HAM
%    U = MXYPHI_1_ESC_HAM(T,IN2,IN3,IN4,U_GAINS1,HSTAR1)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Sep-2018 09:48:36

dtheta1 = in2(8,:);
g = in3(5,:);
lb = in3(4,:);
lp = in3(3,:);
mb = in3(2,:);
mp = in3(1,:);
theta1 = in2(4,:);
u = [0.0;0.0;0.0;U_gains1.*dtheta1.*(-HStar1+dtheta1.^2.*(lb.^2.*mb.*(1.0./2.0)+lp.^2.*mp.*(1.0./2.0))+g.*cos(theta1).*(lb.*mb-lp.*mp))];
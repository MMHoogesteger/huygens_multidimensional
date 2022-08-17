function u = mx_3_esc_ham(t,in2,in3,in4,in5,in6)
%MX_3_ESC_HAM
%    U = MX_3_ESC_HAM(T,IN2,IN3,IN4,IN5,IN6)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    09-Aug-2018 14:00:25

HStar1 = in6(1,:);
HStar2 = in6(2,:);
HStar3 = in6(3,:);
U_gains1 = in5(1,:);
U_gains2 = in5(2,:);
U_gains3 = in5(3,:);
dtheta1 = in2(6,:);
dtheta2 = in2(7,:);
dtheta3 = in2(8,:);
g = in3(5,:);
lb = in3(4,:);
lp = in3(3,:);
mb = in3(2,:);
mp = in3(1,:);
theta1 = in2(2,:);
theta2 = in2(3,:);
theta3 = in2(4,:);
t2 = lb.^2;
t3 = mb.*t2.*(1.0./2.0);
t4 = lp.^2;
t5 = mp.*t4.*(1.0./2.0);
t6 = t3+t5;
t7 = lb.*mb;
t9 = lp.*mp;
t8 = t7-t9;
u = [0.0;U_gains1.*dtheta1.*(-HStar1+dtheta1.^2.*t6+g.*t8.*cos(theta1));U_gains2.*dtheta2.*(-HStar2+dtheta2.^2.*t6+g.*t8.*cos(theta2));U_gains3.*dtheta3.*(-HStar3+dtheta3.^2.*t6+g.*t8.*cos(theta3))];

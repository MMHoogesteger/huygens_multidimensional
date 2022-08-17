function HM = mxyphi_3_ham_m(t,in2,in3,in4)
%MXYPHI_3_HAM_M
%    HM = MXYPHI_3_HAM_M(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Sep-2018 09:49:45

JP33 = in3(18,:);
dphi = in2(9,:);
dtheta1 = in2(10,:);
dtheta2 = in2(11,:);
dtheta3 = in2(12,:);
dx = in2(7,:);
dy = in2(8,:);
g = in3(5,:);
k = in3(9,:);
kt = in3(10,:);
lb = in3(4,:);
lp = in3(3,:);
mP = in3(7,:);
mb = in3(2,:);
mp = in3(1,:);
phi = in2(3,:);
theta1 = in2(4,:);
theta2 = in2(5,:);
theta3 = in2(6,:);
x = in2(1,:);
y = in2(2,:);
t2 = lb.^2;
t3 = mb.*t2.*(1.0./2.0);
t4 = lp.^2;
t5 = mp.*t4.*(1.0./2.0);
t6 = t3+t5;
t7 = lb.*mb;
t9 = lp.*mp;
t8 = t7-t9;
HM = [dx.^2.*mP.*(1.0./2.0)+k.*x.^2.*(1.0./2.0);dy.^2.*mP.*(1.0./2.0)+k.*y.^2.*(1.0./2.0);JP33.*dphi.^2.*(1.0./2.0)+kt.*phi.^2.*(1.0./2.0);dtheta1.^2.*t6+g.*t8.*cos(theta1);dtheta2.^2.*t6+g.*t8.*cos(theta2);dtheta3.^2.*t6+g.*t8.*cos(theta3)];

function H = mphi_1_ham(t,in2,in3,in4)
%MPHI_1_HAM
%    H = MPHI_1_HAM(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    09-Aug-2018 13:59:28

JM33 = in3(19,:);
JP33 = in3(18,:);
R1 = in4(1,:);
dphi = in2(3,:);
dtheta1 = in2(4,:);
g = in3(5,:);
gamma1 = in4(2,:);
kt = in3(10,:);
lb = in3(4,:);
lp = in3(3,:);
mM = in3(8,:);
mb = in3(2,:);
mp = in3(1,:);
phi = in2(1,:);
psi1 = in4(3,:);
theta1 = in2(2,:);
t2 = dphi.^2;
t3 = sin(theta1);
t4 = cos(gamma1);
t5 = sin(psi1);
t6 = cos(psi1);
t7 = cos(theta1);
t8 = sin(gamma1);
t9 = lb.*mb;
t10 = t9-lp.*mp;
t11 = R1.^2;
H = kt.*phi.^2.*(1.0./2.0)+(lb.^2.*mb+lp.^2.*mp).*(t2.*t3.^2+dtheta1.^2).*(1.0./2.0)+t2.*(JM33+JP33).*(1.0./2.0)+t10.*(R1.*t2.*t3.*t4.*t6.*2.0+R1.*t2.*t3.*t5.*t8.*2.0+R1.*dphi.*dtheta1.*t4.*t5.*t7.*2.0-R1.*dphi.*dtheta1.*t6.*t7.*t8.*2.0).*(1.0./2.0)+t2.*(t4.^2.*t11+t8.^2.*t11).*(mM+mb+mp).*(1.0./2.0)+g.*t7.*t10;

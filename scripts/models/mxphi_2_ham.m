function H = mxphi_2_ham(t,in2,in3,in4)
%MXPHI_2_HAM
%    H = MXPHI_2_HAM(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    09-Aug-2018 13:59:54

JM33 = in3(19,:);
JP33 = in3(18,:);
R1 = in4(1,:);
R2 = in4(2,:);
dphi = in2(6,:);
dtheta1 = in2(7,:);
dtheta2 = in2(8,:);
dx = in2(5,:);
g = in3(5,:);
gamma1 = in4(3,:);
gamma2 = in4(4,:);
k = in3(9,:);
kt = in3(10,:);
lb = in3(4,:);
lp = in3(3,:);
mM = in3(8,:);
mP = in3(7,:);
mb = in3(2,:);
mp = in3(1,:);
phi = in2(2,:);
psi1 = in4(5,:);
psi2 = in4(6,:);
theta1 = in2(3,:);
theta2 = in2(4,:);
x = in2(1,:);
t2 = cos(gamma1);
t3 = R1.^2;
t4 = sin(gamma1);
t5 = dphi.^2;
t6 = cos(gamma2);
t7 = R2.^2;
t8 = sin(gamma2);
t9 = sin(phi);
t10 = cos(phi);
t11 = phi+psi1;
t12 = phi+psi2;
t13 = sin(theta1);
t14 = sin(theta2);
t15 = cos(theta1);
t16 = sin(psi1);
t17 = cos(psi1);
t18 = cos(theta2);
t19 = sin(psi2);
t20 = cos(psi2);
t21 = lb.*mb;
t22 = t21-lp.*mp;
H = (lb.^2.*mb+lp.^2.*mp).*(t5.*t13.^2+t5.*t14.^2+dtheta1.^2+dtheta2.^2).*(1.0./2.0)+t5.*(JM33.*2.0+JP33).*(1.0./2.0)+dx.^2.*(mM.*2.0+mP+mb.*2.0+mp.*2.0).*(1.0./2.0)+t22.*(dtheta1.*dx.*t15.*cos(t11).*2.0+dtheta2.*dx.*t18.*cos(t12).*2.0-dphi.*dx.*t13.*sin(t11).*2.0-dphi.*dx.*t14.*sin(t12).*2.0+R1.*t2.*t5.*t13.*t17.*2.0+R1.*t4.*t5.*t13.*t16.*2.0+R2.*t5.*t6.*t14.*t20.*2.0+R2.*t5.*t8.*t14.*t19.*2.0+R1.*dphi.*dtheta1.*t2.*t15.*t16.*2.0-R1.*dphi.*dtheta1.*t4.*t15.*t17.*2.0+R2.*dphi.*dtheta2.*t6.*t18.*t19.*2.0-R2.*dphi.*dtheta2.*t8.*t18.*t20.*2.0).*(1.0./2.0)+kt.*phi.^2.*(1.0./2.0)+k.*x.^2.*(1.0./2.0)+(mM+mb+mp).*(t5.*(t2.^2.*t3+t3.*t4.^2)+t5.*(t6.^2.*t7+t7.*t8.^2)-dphi.*dx.*(R1.*t2.*t9+R1.*t4.*t10).*2.0-dphi.*dx.*(R2.*t6.*t9+R2.*t8.*t10).*2.0).*(1.0./2.0)+g.*t22.*(t15+t18);

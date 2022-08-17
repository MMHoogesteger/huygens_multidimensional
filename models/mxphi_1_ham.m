function H = mxphi_1_ham(t,in2,in3,in4)
%MXPHI_1_HAM
%    H = MXPHI_1_HAM(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    09-Aug-2018 13:59:33

JM33 = in3(19,:);
JP33 = in3(18,:);
R1 = in4(1,:);
dphi = in2(5,:);
dtheta1 = in2(6,:);
dx = in2(4,:);
g = in3(5,:);
gamma1 = in4(2,:);
k = in3(9,:);
kt = in3(10,:);
lb = in3(4,:);
lp = in3(3,:);
mM = in3(8,:);
mP = in3(7,:);
mb = in3(2,:);
mp = in3(1,:);
phi = in2(2,:);
psi1 = in4(3,:);
theta1 = in2(3,:);
x = in2(1,:);
t2 = phi+psi1;
t3 = sin(theta1);
t4 = dphi.^2;
t5 = cos(gamma1);
t6 = cos(theta1);
t7 = sin(psi1);
t8 = cos(psi1);
t9 = sin(gamma1);
t10 = R1.^2;
t11 = lb.*mb;
t12 = t11-lp.*mp;
H = dx.^2.*(mM+mP+mb+mp).*(1.0./2.0)+(t4.*(t5.^2.*t10+t9.^2.*t10)-dphi.*dx.*(R1.*t5.*sin(phi)+R1.*t9.*cos(phi)).*2.0).*(mM+mb+mp).*(1.0./2.0)+kt.*phi.^2.*(1.0./2.0)+k.*x.^2.*(1.0./2.0)+(lb.^2.*mb+lp.^2.*mp).*(t3.^2.*t4+dtheta1.^2).*(1.0./2.0)+t4.*(JM33+JP33).*(1.0./2.0)+t12.*(dtheta1.*dx.*t6.*cos(t2).*2.0-dphi.*dx.*t3.*sin(t2).*2.0+R1.*t3.*t4.*t5.*t8.*2.0+R1.*t3.*t4.*t7.*t9.*2.0+R1.*dphi.*dtheta1.*t5.*t6.*t7.*2.0-R1.*dphi.*dtheta1.*t6.*t8.*t9.*2.0).*(1.0./2.0)+g.*t6.*t12;

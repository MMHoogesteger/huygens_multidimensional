function M = mxphi_1_mass(t,in2,in3,in4)
%MXPHI_1_MASS
%    M = MXPHI_1_MASS(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    09-Aug-2018 13:59:32

JM33 = in3(19,:);
JP33 = in3(18,:);
R1 = in4(1,:);
gamma1 = in4(2,:);
lb = in3(4,:);
lp = in3(3,:);
mM = in3(8,:);
mP = in3(7,:);
mb = in3(2,:);
mp = in3(1,:);
phi = in2(2,:);
psi1 = in4(3,:);
theta1 = in2(3,:);
t2 = phi+psi1;
t3 = lb.*mb;
t10 = lp.*mp;
t4 = t3-t10;
t5 = gamma1+phi;
t6 = sin(t5);
t7 = mM+mb+mp;
t8 = sin(t2);
t9 = sin(theta1);
t11 = -R1.*t6.*t7-t4.*t8.*t9;
t12 = gamma1-psi1;
t13 = cos(theta1);
t14 = cos(t2);
t15 = t4.*t13.*t14;
t16 = sin(t12);
t17 = lb.^2;
t18 = mb.*t17;
t19 = lp.^2;
t20 = mp.*t19;
t21 = t18+t20;
M = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,mM+mP+mb+mp,t11,t15,0.0,0.0,0.0,t11,JM33+JP33+R1.^2.*t7+t9.^2.*t21+R1.*t9.*cos(t12).*(lb.*mb.*2.0-lp.*mp.*2.0),-R1.*t4.*t13.*t16,0.0,0.0,0.0,t15,-R1.*t4.*t13.*t16,t21],[6,6]);

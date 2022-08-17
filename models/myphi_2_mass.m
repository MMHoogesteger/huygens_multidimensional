function M = myphi_2_mass(t,in2,in3,in4)
%MYPHI_2_MASS
%    M = MYPHI_2_MASS(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    09-Aug-2018 13:59:56

JM33 = in3(19,:);
JP33 = in3(18,:);
R1 = in4(1,:);
R2 = in4(2,:);
gamma1 = in4(3,:);
gamma2 = in4(4,:);
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
t2 = lb.*mb;
t6 = lp.*mp;
t3 = t2-t6;
t4 = mM+mb+mp;
t5 = phi+psi1;
t7 = phi+psi2;
t8 = cos(t5);
t9 = sin(theta1);
t10 = t3.*t8.*t9;
t11 = cos(t7);
t12 = sin(theta2);
t13 = t3.*t11.*t12;
t14 = gamma1+phi;
t15 = cos(t14);
t16 = R1.*t4.*t15;
t17 = gamma2+phi;
t18 = cos(t17);
t19 = R2.*t4.*t18;
t20 = t10+t13+t16+t19;
t21 = lb.^2;
t22 = mb.*t21;
t23 = lp.^2;
t24 = mp.*t23;
t25 = t22+t24;
t26 = lb.*mb.*2.0;
t27 = t26-lp.*mp.*2.0;
t28 = gamma1-psi1;
t29 = cos(theta1);
t30 = gamma2-psi2;
t31 = cos(theta2);
t32 = sin(t5);
t33 = t3.*t29.*t32;
t34 = sin(t28);
t35 = sin(t7);
t36 = t3.*t31.*t35;
t37 = sin(t30);
M = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,mM.*2.0+mP+mb.*2.0+mp.*2.0,t20,t33,t36,0.0,0.0,0.0,0.0,t20,JM33.*2.0+JP33+R1.^2.*t4+R2.^2.*t4+t9.^2.*t25+t12.^2.*t25+R1.*t9.*t27.*cos(t28)+R2.*t12.*t27.*cos(t30),-R1.*t3.*t29.*t34,-R2.*t3.*t31.*t37,0.0,0.0,0.0,0.0,t33,-R1.*t3.*t29.*t34,t25,0.0,0.0,0.0,0.0,0.0,t36,-R2.*t3.*t31.*t37,0.0,t25],[8,8]);

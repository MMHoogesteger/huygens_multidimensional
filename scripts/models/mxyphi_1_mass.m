function M = mxyphi_1_mass(t,in2,in3,in4)
%MXYPHI_1_MASS
%    M = MXYPHI_1_MASS(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Sep-2018 09:48:35

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
phi = in2(3,:);
psi1 = in4(3,:);
theta1 = in2(4,:);
t2 = phi+psi1;
t3 = lb.*mb;
t8 = lp.*mp;
t4 = t3-t8;
t5 = mM+mP+mb+mp;
t6 = cos(t2);
t7 = sin(theta1);
t9 = gamma1+phi;
t10 = mM+mb+mp;
t11 = sin(t2);
t12 = cos(theta1);
t13 = sin(t9);
t14 = -R1.*t10.*t13-t4.*t7.*t11;
t15 = t4.*t6.*t7;
t16 = cos(t9);
t17 = R1.*t10.*t16;
t18 = t15+t17;
t19 = gamma1-psi1;
t20 = t4.*t6.*t12;
t21 = t4.*t11.*t12;
t22 = sin(t19);
t23 = lb.^2;
t24 = mb.*t23;
t25 = lp.^2;
t26 = mp.*t25;
t27 = t24+t26;
M = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t5,0.0,t14,t20,0.0,0.0,0.0,0.0,0.0,t5,t18,t21,0.0,0.0,0.0,0.0,t14,t18,JM33+JP33+R1.^2.*t10+t7.^2.*t27+R1.*t7.*cos(t19).*(lb.*mb.*2.0-lp.*mp.*2.0),-R1.*t4.*t12.*t22,0.0,0.0,0.0,0.0,t20,t21,-R1.*t4.*t12.*t22,t27],[8,8]);

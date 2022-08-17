function M = mxyphi_3_mass(t,in2,in3,in4)
%MXYPHI_3_MASS
%    M = MXYPHI_3_MASS(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Sep-2018 09:49:42

JM33 = in3(19,:);
JP33 = in3(18,:);
R1 = in4(1,:);
R2 = in4(2,:);
R3 = in4(3,:);
gamma1 = in4(4,:);
gamma2 = in4(5,:);
gamma3 = in4(6,:);
lb = in3(4,:);
lp = in3(3,:);
mM = in3(8,:);
mP = in3(7,:);
mb = in3(2,:);
mp = in3(1,:);
phi = in2(3,:);
psi1 = in4(7,:);
psi2 = in4(8,:);
psi3 = in4(9,:);
theta1 = in2(4,:);
theta2 = in2(5,:);
theta3 = in2(6,:);
t2 = mM+mb+mp;
t3 = lb.*mb;
t5 = lp.*mp;
t4 = t3-t5;
t6 = phi+psi1;
t7 = phi+psi2;
t8 = phi+psi3;
t9 = mM.*3.0;
t10 = mb.*3.0;
t11 = mp.*3.0;
t12 = mP+t9+t10+t11;
t13 = cos(t6);
t14 = sin(theta1);
t15 = cos(t7);
t16 = sin(theta2);
t17 = cos(t8);
t18 = sin(theta3);
t19 = gamma1+phi;
t20 = gamma2+phi;
t21 = gamma3+phi;
t22 = sin(t6);
t23 = cos(theta1);
t24 = sin(t7);
t25 = cos(theta2);
t26 = sin(t8);
t27 = cos(theta3);
t28 = sin(t19);
t29 = sin(t20);
t30 = sin(t21);
t31 = -R1.*t2.*t28-R2.*t2.*t29-R3.*t2.*t30-t4.*t14.*t22-t4.*t16.*t24-t4.*t18.*t26;
t32 = t4.*t13.*t14;
t33 = t4.*t15.*t16;
t34 = t4.*t17.*t18;
t35 = cos(t19);
t36 = R1.*t2.*t35;
t37 = cos(t20);
t38 = R2.*t2.*t37;
t39 = cos(t21);
t40 = R3.*t2.*t39;
t41 = t32+t33+t34+t36+t38+t40;
t42 = lb.^2;
t43 = mb.*t42;
t44 = lp.^2;
t45 = mp.*t44;
t46 = t43+t45;
t47 = lb.*mb.*2.0;
t49 = lp.*mp.*2.0;
t48 = t47-t49;
t50 = gamma1-psi1;
t51 = gamma2-psi2;
t52 = gamma3-psi3;
t53 = t4.*t13.*t23;
t54 = t4.*t22.*t23;
t55 = sin(t50);
t56 = t4.*t15.*t25;
t57 = t4.*t24.*t25;
t58 = sin(t51);
t59 = t4.*t17.*t27;
t60 = t4.*t26.*t27;
t61 = sin(t52);
M = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t12,0.0,t31,t53,t56,t59,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t12,t41,t54,t57,t60,0.0,0.0,0.0,0.0,0.0,0.0,t31,t41,JM33.*3.0+JP33+R1.^2.*t2+R2.^2.*t2+R3.^2.*t2+t14.^2.*t46+t16.^2.*t46+t18.^2.*t46+R1.*t14.*t48.*cos(t50)+R2.*t16.*t48.*cos(t51)+R3.*t18.*t48.*cos(t52),-R1.*t4.*t23.*t55,-R2.*t4.*t25.*t58,-R3.*t4.*t27.*t61,0.0,0.0,0.0,0.0,0.0,0.0,t53,t54,-R1.*t4.*t23.*t55,t46,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t56,t57,-R2.*t4.*t25.*t58,0.0,t46,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t59,t60,-R3.*t4.*t27.*t61,0.0,0.0,t46],[12,12]);

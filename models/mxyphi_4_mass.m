function M = mxyphi_4_mass(t,in2,in3,in4)
%MXYPHI_4_MASS
%    M = MXYPHI_4_MASS(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Sep-2018 09:50:19

JM33 = in3(19,:);
JP33 = in3(18,:);
R1 = in4(1,:);
R2 = in4(2,:);
R3 = in4(3,:);
R4 = in4(4,:);
gamma1 = in4(5,:);
gamma2 = in4(6,:);
gamma3 = in4(7,:);
gamma4 = in4(8,:);
lb = in3(4,:);
lp = in3(3,:);
mM = in3(8,:);
mP = in3(7,:);
mb = in3(2,:);
mp = in3(1,:);
phi = in2(3,:);
psi1 = in4(9,:);
psi2 = in4(10,:);
psi3 = in4(11,:);
psi4 = in4(12,:);
theta1 = in2(4,:);
theta2 = in2(5,:);
theta3 = in2(6,:);
theta4 = in2(7,:);
t2 = mM+mb+mp;
t3 = lb.*mb;
t5 = lp.*mp;
t4 = t3-t5;
t6 = phi+psi1;
t7 = phi+psi2;
t8 = phi+psi3;
t9 = phi+psi4;
t10 = mM.*4.0;
t11 = mb.*4.0;
t12 = mp.*4.0;
t13 = mP+t10+t11+t12;
t14 = cos(t6);
t15 = sin(theta1);
t16 = cos(t7);
t17 = sin(theta2);
t18 = cos(t8);
t19 = sin(theta3);
t20 = cos(t9);
t21 = sin(theta4);
t22 = gamma1+phi;
t23 = gamma2+phi;
t24 = gamma3+phi;
t25 = gamma4+phi;
t26 = sin(t6);
t27 = cos(theta1);
t28 = sin(t7);
t29 = cos(theta2);
t30 = sin(t8);
t31 = cos(theta3);
t32 = sin(t9);
t33 = cos(theta4);
t34 = sin(t22);
t35 = sin(t23);
t36 = sin(t24);
t37 = sin(t25);
t38 = -R1.*t2.*t34-R2.*t2.*t35-R3.*t2.*t36-R4.*t2.*t37-t4.*t15.*t26-t4.*t17.*t28-t4.*t19.*t30-t4.*t21.*t32;
t39 = t4.*t14.*t15;
t40 = t4.*t16.*t17;
t41 = t4.*t18.*t19;
t42 = t4.*t20.*t21;
t43 = cos(t22);
t44 = R1.*t2.*t43;
t45 = cos(t23);
t46 = R2.*t2.*t45;
t47 = cos(t24);
t48 = R3.*t2.*t47;
t49 = cos(t25);
t50 = R4.*t2.*t49;
t51 = t39+t40+t41+t42+t44+t46+t48+t50;
t52 = lb.^2;
t53 = mb.*t52;
t54 = lp.^2;
t55 = mp.*t54;
t56 = t53+t55;
t57 = lb.*mb.*2.0;
t59 = lp.*mp.*2.0;
t58 = t57-t59;
t60 = gamma1-psi1;
t61 = gamma2-psi2;
t62 = gamma3-psi3;
t63 = gamma4-psi4;
t64 = t4.*t14.*t27;
t65 = t4.*t26.*t27;
t66 = sin(t60);
t67 = t4.*t16.*t29;
t68 = t4.*t28.*t29;
t69 = sin(t61);
t70 = t4.*t18.*t31;
t71 = t4.*t30.*t31;
t72 = sin(t62);
t73 = t4.*t20.*t33;
t74 = t4.*t32.*t33;
t75 = sin(t63);
M = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t13,0.0,t38,t64,t67,t70,t73,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t13,t51,t65,t68,t71,t74,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t38,t51,JM33.*4.0+JP33+R1.^2.*t2+R2.^2.*t2+R3.^2.*t2+R4.^2.*t2+t15.^2.*t56+t17.^2.*t56+t19.^2.*t56+t21.^2.*t56+R1.*t15.*t58.*cos(t60)+R2.*t17.*t58.*cos(t61)+R3.*t19.*t58.*cos(t62)+R4.*t21.*t58.*cos(t63),-R1.*t4.*t27.*t66,-R2.*t4.*t29.*t69,-R3.*t4.*t31.*t72,-R4.*t4.*t33.*t75,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t64,t65,-R1.*t4.*t27.*t66,t56,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t67,t68,-R2.*t4.*t29.*t69,0.0,t56,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t70,t71,-R3.*t4.*t31.*t72,0.0,0.0,t56,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,t73,t74,-R4.*t4.*t33.*t75,0.0,0.0,0.0,t56],[14,14]);

function H = mxyphi_10_ham(t,in2,in3,in4)
%MXYPHI_10_HAM
%    H = MXYPHI_10_HAM(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Sep-2018 09:47:19

JM33 = in3(19,:);
JP33 = in3(18,:);
R1 = in4(1,:);
R2 = in4(2,:);
R3 = in4(3,:);
R4 = in4(4,:);
R5 = in4(5,:);
R6 = in4(6,:);
R7 = in4(7,:);
R8 = in4(8,:);
R9 = in4(9,:);
R10 = in4(10,:);
dphi = in2(16,:);
dtheta1 = in2(17,:);
dtheta2 = in2(18,:);
dtheta3 = in2(19,:);
dtheta4 = in2(20,:);
dtheta5 = in2(21,:);
dtheta6 = in2(22,:);
dtheta7 = in2(23,:);
dtheta8 = in2(24,:);
dtheta9 = in2(25,:);
dtheta10 = in2(26,:);
dx = in2(14,:);
dy = in2(15,:);
g = in3(5,:);
gamma1 = in4(11,:);
gamma2 = in4(12,:);
gamma3 = in4(13,:);
gamma4 = in4(14,:);
gamma5 = in4(15,:);
gamma6 = in4(16,:);
gamma7 = in4(17,:);
gamma8 = in4(18,:);
gamma9 = in4(19,:);
gamma10 = in4(20,:);
k = in3(9,:);
kt = in3(10,:);
lb = in3(4,:);
lp = in3(3,:);
mM = in3(8,:);
mP = in3(7,:);
mb = in3(2,:);
mp = in3(1,:);
phi = in2(3,:);
psi1 = in4(21,:);
psi2 = in4(22,:);
psi3 = in4(23,:);
psi4 = in4(24,:);
psi5 = in4(25,:);
psi6 = in4(26,:);
psi7 = in4(27,:);
psi8 = in4(28,:);
psi9 = in4(29,:);
psi10 = in4(30,:);
theta1 = in2(4,:);
theta2 = in2(5,:);
theta3 = in2(6,:);
theta4 = in2(7,:);
theta5 = in2(8,:);
theta6 = in2(9,:);
theta7 = in2(10,:);
theta8 = in2(11,:);
theta9 = in2(12,:);
theta10 = in2(13,:);
x = in2(1,:);
y = in2(2,:);
t2 = sin(theta1);
t3 = dphi.^2;
t4 = sin(theta2);
t5 = sin(theta3);
t6 = sin(theta4);
t7 = sin(theta5);
t8 = sin(theta6);
t9 = sin(theta7);
t10 = sin(theta8);
t11 = sin(theta9);
t12 = sin(theta10);
t13 = cos(gamma1);
t14 = R1.^2;
t15 = sin(gamma1);
t16 = cos(gamma2);
t17 = R2.^2;
t18 = sin(gamma2);
t19 = cos(gamma3);
t20 = R3.^2;
t21 = sin(gamma3);
t22 = cos(gamma4);
t23 = R4.^2;
t24 = sin(gamma4);
t25 = cos(gamma5);
t26 = R5.^2;
t27 = sin(gamma5);
t28 = cos(gamma6);
t29 = R6.^2;
t30 = sin(gamma6);
t31 = cos(gamma7);
t32 = R7.^2;
t33 = sin(gamma7);
t34 = cos(gamma8);
t35 = R8.^2;
t36 = sin(gamma8);
t37 = cos(gamma9);
t38 = R9.^2;
t39 = sin(gamma9);
t40 = cos(gamma10);
t41 = R10.^2;
t42 = sin(gamma10);
t43 = sin(phi);
t44 = cos(phi);
t45 = phi+psi1;
t46 = cos(t45);
t47 = phi+psi2;
t48 = cos(t47);
t49 = phi+psi3;
t50 = cos(t49);
t51 = phi+psi4;
t52 = cos(t51);
t53 = phi+psi5;
t54 = cos(t53);
t55 = phi+psi6;
t56 = cos(t55);
t57 = phi+psi7;
t58 = cos(t57);
t59 = phi+psi8;
t60 = cos(t59);
t61 = phi+psi9;
t62 = cos(t61);
t63 = phi+psi10;
t64 = cos(t63);
t65 = cos(theta1);
t66 = cos(theta2);
t67 = cos(theta3);
t68 = cos(theta4);
t69 = cos(theta5);
t70 = cos(theta6);
t71 = cos(theta7);
t72 = cos(theta8);
t73 = cos(theta9);
t74 = cos(theta10);
t75 = sin(t45);
t76 = sin(t47);
t77 = sin(t49);
t78 = sin(t51);
t79 = sin(t53);
t80 = sin(t55);
t81 = sin(t57);
t82 = sin(t59);
t83 = sin(t61);
t84 = sin(t63);
t85 = sin(psi1);
t86 = cos(psi1);
t87 = sin(psi2);
t88 = cos(psi2);
t89 = sin(psi3);
t90 = cos(psi3);
t91 = sin(psi4);
t92 = cos(psi4);
t93 = sin(psi5);
t94 = cos(psi5);
t95 = sin(psi6);
t96 = cos(psi6);
t97 = sin(psi7);
t98 = cos(psi7);
t99 = sin(psi8);
t100 = cos(psi8);
t101 = sin(psi9);
t102 = cos(psi9);
t103 = sin(psi10);
t104 = cos(psi10);
t105 = lb.*mb;
t106 = t105-lp.*mp;
H = (mM+mb+mp).*(t3.*(t13.^2.*t14+t14.*t15.^2)+t3.*(t16.^2.*t17+t17.*t18.^2)+t3.*(t19.^2.*t20+t20.*t21.^2)+t3.*(t22.^2.*t23+t23.*t24.^2)+t3.*(t25.^2.*t26+t26.*t27.^2)+t3.*(t28.^2.*t29+t29.*t30.^2)+t3.*(t31.^2.*t32+t32.*t33.^2)+t3.*(t34.^2.*t35+t35.*t36.^2)+t3.*(t37.^2.*t38+t38.*t39.^2)+t3.*(t40.^2.*t41+t41.*t42.^2)-dphi.*dx.*(R1.*t13.*t43+R1.*t15.*t44).*2.0-dphi.*dx.*(R2.*t16.*t43+R2.*t18.*t44).*2.0-dphi.*dx.*(R3.*t19.*t43+R3.*t21.*t44).*2.0-dphi.*dx.*(R4.*t22.*t43+R4.*t24.*t44).*2.0-dphi.*dx.*(R5.*t25.*t43+R5.*t27.*t44).*2.0-dphi.*dx.*(R6.*t28.*t43+R6.*t30.*t44).*2.0-dphi.*dx.*(R7.*t31.*t43+R7.*t33.*t44).*2.0-dphi.*dx.*(R8.*t34.*t43+R8.*t36.*t44).*2.0-dphi.*dx.*(R9.*t37.*t43+R9.*t39.*t44).*2.0-dphi.*dx.*(R10.*t40.*t43+R10.*t42.*t44).*2.0+dphi.*dy.*(R1.*t13.*t44-R1.*t15.*t43).*2.0+dphi.*dy.*(R2.*t16.*t44-R2.*t18.*t43).*2.0+dphi.*dy.*(R3.*t19.*t44-R3.*t21.*t43).*2.0+dphi.*dy.*(R4.*t22.*t44-R4.*t24.*t43).*2.0+dphi.*dy.*(R5.*t25.*t44-R5.*t27.*t43).*2.0+dphi.*dy.*(R6.*t28.*t44-R6.*t30.*t43).*2.0+dphi.*dy.*(R7.*t31.*t44-R7.*t33.*t43).*2.0+dphi.*dy.*(R8.*t34.*t44-R8.*t36.*t43).*2.0+dphi.*dy.*(R9.*t37.*t44-R9.*t39.*t43).*2.0+dphi.*dy.*(R10.*t40.*t44-R10.*t42.*t43).*2.0).*(1.0./2.0)+t3.*(JM33.*1.0e1+JP33).*(1.0./2.0)+(dx.^2+dy.^2).*(mM.*1.0e1+mP+mb.*1.0e1+mp.*1.0e1).*(1.0./2.0)+(lb.^2.*mb+lp.^2.*mp).*(t2.^2.*t3+t3.*t4.^2+t3.*t5.^2+t3.*t6.^2+t3.*t7.^2+t3.*t8.^2+t3.*t9.^2+t3.*t10.^2+t3.*t11.^2+t3.*t12.^2+dtheta1.^2+dtheta2.^2+dtheta3.^2+dtheta4.^2+dtheta5.^2+dtheta6.^2+dtheta7.^2+dtheta8.^2+dtheta9.^2+dtheta10.^2).*(1.0./2.0)+kt.*phi.^2.*(1.0./2.0)+t106.*(dphi.*dx.*t2.*t75.*-2.0-dphi.*dx.*t4.*t76.*2.0-dphi.*dx.*t5.*t77.*2.0-dphi.*dx.*t6.*t78.*2.0-dphi.*dx.*t7.*t79.*2.0-dphi.*dx.*t8.*t80.*2.0-dphi.*dx.*t9.*t81.*2.0-dphi.*dx.*t10.*t82.*2.0-dphi.*dx.*t11.*t83.*2.0-dphi.*dx.*t12.*t84.*2.0+dphi.*dy.*t2.*t46.*2.0+dphi.*dy.*t4.*t48.*2.0+dphi.*dy.*t5.*t50.*2.0+dphi.*dy.*t6.*t52.*2.0+dphi.*dy.*t7.*t54.*2.0+dphi.*dy.*t8.*t56.*2.0+dphi.*dy.*t9.*t58.*2.0+dphi.*dy.*t10.*t60.*2.0+dphi.*dy.*t11.*t62.*2.0+dphi.*dy.*t12.*t64.*2.0+dtheta1.*dx.*t46.*t65.*2.0+dtheta2.*dx.*t48.*t66.*2.0+dtheta3.*dx.*t50.*t67.*2.0+dtheta4.*dx.*t52.*t68.*2.0+dtheta5.*dx.*t54.*t69.*2.0+dtheta6.*dx.*t56.*t70.*2.0+dtheta7.*dx.*t58.*t71.*2.0+dtheta8.*dx.*t60.*t72.*2.0+dtheta9.*dx.*t62.*t73.*2.0+dtheta10.*dx.*t64.*t74.*2.0+dtheta1.*dy.*t65.*t75.*2.0+dtheta2.*dy.*t66.*t76.*2.0+dtheta3.*dy.*t67.*t77.*2.0+dtheta4.*dy.*t68.*t78.*2.0+dtheta5.*dy.*t69.*t79.*2.0+dtheta6.*dy.*t70.*t80.*2.0+dtheta7.*dy.*t71.*t81.*2.0+dtheta8.*dy.*t72.*t82.*2.0+dtheta9.*dy.*t73.*t83.*2.0+dtheta10.*dy.*t74.*t84.*2.0+R1.*t2.*t3.*t13.*t86.*2.0+R1.*t2.*t3.*t15.*t85.*2.0+R2.*t3.*t4.*t16.*t88.*2.0+R2.*t3.*t4.*t18.*t87.*2.0+R3.*t3.*t5.*t19.*t90.*2.0+R3.*t3.*t5.*t21.*t89.*2.0+R4.*t3.*t6.*t22.*t92.*2.0+R4.*t3.*t6.*t24.*t91.*2.0+R5.*t3.*t7.*t25.*t94.*2.0+R5.*t3.*t7.*t27.*t93.*2.0+R6.*t3.*t8.*t28.*t96.*2.0+R6.*t3.*t8.*t30.*t95.*2.0+R7.*t3.*t9.*t31.*t98.*2.0+R7.*t3.*t9.*t33.*t97.*2.0+R8.*t3.*t10.*t34.*t100.*2.0+R8.*t3.*t10.*t36.*t99.*2.0+R9.*t3.*t11.*t37.*t102.*2.0+R9.*t3.*t11.*t39.*t101.*2.0+R10.*t3.*t12.*t40.*t104.*2.0+R10.*t3.*t12.*t42.*t103.*2.0+R1.*dphi.*dtheta1.*t13.*t65.*t85.*2.0-R1.*dphi.*dtheta1.*t15.*t65.*t86.*2.0+R2.*dphi.*dtheta2.*t16.*t66.*t87.*2.0-R2.*dphi.*dtheta2.*t18.*t66.*t88.*2.0+R3.*dphi.*dtheta3.*t19.*t67.*t89.*2.0-R3.*dphi.*dtheta3.*t21.*t67.*t90.*2.0+R4.*dphi.*dtheta4.*t22.*t68.*t91.*2.0-R4.*dphi.*dtheta4.*t24.*t68.*t92.*2.0+R5.*dphi.*dtheta5.*t25.*t69.*t93.*2.0-R5.*dphi.*dtheta5.*t27.*t69.*t94.*2.0+R6.*dphi.*dtheta6.*t28.*t70.*t95.*2.0-R6.*dphi.*dtheta6.*t30.*t70.*t96.*2.0+R7.*dphi.*dtheta7.*t31.*t71.*t97.*2.0-R7.*dphi.*dtheta7.*t33.*t71.*t98.*2.0+R8.*dphi.*dtheta8.*t34.*t72.*t99.*2.0-R8.*dphi.*dtheta8.*t36.*t72.*t100.*2.0+R9.*dphi.*dtheta9.*t37.*t73.*t101.*2.0-R9.*dphi.*dtheta9.*t39.*t73.*t102.*2.0+R10.*dphi.*dtheta10.*t40.*t74.*t103.*2.0-R10.*dphi.*dtheta10.*t42.*t74.*t104.*2.0).*(1.0./2.0)+k.*(x.^2+y.^2).*(1.0./2.0)+g.*t106.*(t65+t66+t67+t68+t69+t70+t71+t72+t73+t74);

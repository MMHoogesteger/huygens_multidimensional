function dxdt = mxyphi_10_row(t,in2,in3,in4,in5)
%MXYPHI_10_ROW
%    DXDT = MXYPHI_10_ROW(T,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Sep-2018 09:47:15

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
Uphi = in5(3,:);
Utheta1 = in5(4,:);
Utheta2 = in5(5,:);
Utheta3 = in5(6,:);
Utheta4 = in5(7,:);
Utheta5 = in5(8,:);
Utheta6 = in5(9,:);
Utheta7 = in5(10,:);
Utheta8 = in5(11,:);
Utheta9 = in5(12,:);
Utheta10 = in5(13,:);
Ux = in5(1,:);
Uy = in5(2,:);
c = in3(12,:);
ct = in3(13,:);
d = in3(6,:);
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
t2 = dphi.^2;
t3 = phi+psi1;
t4 = cos(t3);
t5 = sin(theta1);
t6 = phi+psi2;
t7 = cos(t6);
t8 = sin(theta2);
t9 = phi+psi3;
t10 = cos(t9);
t11 = sin(theta3);
t12 = phi+psi4;
t13 = cos(t12);
t14 = sin(theta4);
t15 = phi+psi5;
t16 = cos(t15);
t17 = sin(theta5);
t18 = phi+psi6;
t19 = cos(t18);
t20 = sin(theta6);
t21 = phi+psi7;
t22 = cos(t21);
t23 = sin(theta7);
t24 = phi+psi8;
t25 = cos(t24);
t26 = sin(theta8);
t27 = phi+psi9;
t28 = cos(t27);
t29 = sin(theta9);
t30 = phi+psi10;
t31 = cos(t30);
t32 = sin(theta10);
t33 = mM+mb+mp;
t34 = gamma1+phi;
t35 = gamma2+phi;
t36 = gamma3+phi;
t37 = gamma4+phi;
t38 = gamma5+phi;
t39 = gamma6+phi;
t40 = gamma7+phi;
t41 = gamma8+phi;
t42 = gamma9+phi;
t43 = gamma10+phi;
t44 = lb.*mb;
t76 = lp.*mp;
t45 = t44-t76;
t46 = sin(t3);
t47 = sin(t6);
t48 = sin(t9);
t49 = sin(t12);
t50 = sin(t15);
t51 = sin(t18);
t52 = sin(t21);
t53 = sin(t24);
t54 = sin(t27);
t55 = sin(t30);
t56 = dtheta1.^2;
t57 = dtheta2.^2;
t58 = dtheta3.^2;
t59 = dtheta4.^2;
t60 = dtheta5.^2;
t61 = dtheta6.^2;
t62 = dtheta7.^2;
t63 = dtheta8.^2;
t64 = dtheta9.^2;
t65 = dtheta10.^2;
t66 = cos(theta1);
t67 = cos(theta2);
t68 = cos(theta3);
t69 = cos(theta4);
t70 = cos(theta5);
t71 = cos(theta6);
t72 = cos(theta7);
t73 = cos(theta8);
t74 = cos(theta9);
t75 = cos(theta10);
t77 = gamma1-psi1;
t78 = gamma2-psi2;
t79 = gamma3-psi3;
t80 = gamma4-psi4;
t81 = gamma5-psi5;
t82 = gamma6-psi6;
t83 = gamma7-psi7;
t84 = gamma8-psi8;
t85 = gamma9-psi9;
t86 = gamma10-psi10;
t87 = lb.^2;
t88 = lp.^2;
t89 = cos(t77);
t90 = mb.*t87;
t91 = mp.*t88;
t92 = t90+t91;
t93 = cos(t78);
t94 = cos(t79);
t95 = cos(t80);
t96 = cos(t81);
t97 = cos(t82);
t98 = cos(t83);
t99 = cos(t84);
t100 = cos(t85);
t101 = cos(t86);
dxdt = [dx;dy;dphi;dtheta1;dtheta2;dtheta3;dtheta4;dtheta5;dtheta6;dtheta7;dtheta8;dtheta9;dtheta10;Ux+t33.*(R1.*t2.*cos(t34)+R2.*t2.*cos(t35)+R3.*t2.*cos(t36)+R4.*t2.*cos(t37)+R5.*t2.*cos(t38)+R6.*t2.*cos(t39)+R7.*t2.*cos(t40)+R8.*t2.*cos(t41)+R9.*t2.*cos(t42)+R10.*t2.*cos(t43))-c.*dx-k.*x+t45.*(t2.*t4.*t5+t2.*t7.*t8+t2.*t10.*t11+t2.*t13.*t14+t2.*t16.*t17+t2.*t19.*t20+t2.*t22.*t23+t2.*t25.*t26+t2.*t28.*t29+t2.*t31.*t32+t4.*t5.*t56+t7.*t8.*t57+t10.*t11.*t58+t13.*t14.*t59+t16.*t17.*t60+t19.*t20.*t61+t22.*t23.*t62+t25.*t26.*t63+t28.*t29.*t64+t31.*t32.*t65+dphi.*dtheta1.*t46.*t66.*2.0+dphi.*dtheta2.*t47.*t67.*2.0+dphi.*dtheta3.*t48.*t68.*2.0+dphi.*dtheta4.*t49.*t69.*2.0+dphi.*dtheta5.*t50.*t70.*2.0+dphi.*dtheta6.*t51.*t71.*2.0+dphi.*dtheta7.*t52.*t72.*2.0+dphi.*dtheta8.*t53.*t73.*2.0+dphi.*dtheta9.*t54.*t74.*2.0+dphi.*dtheta10.*t55.*t75.*2.0);Uy-c.*dy-k.*y+t45.*(t2.*t5.*t46+t2.*t8.*t47+t2.*t11.*t48+t2.*t14.*t49+t2.*t17.*t50+t2.*t20.*t51+t2.*t23.*t52+t2.*t26.*t53+t2.*t29.*t54+t2.*t32.*t55+t5.*t46.*t56+t8.*t47.*t57+t11.*t48.*t58+t14.*t49.*t59+t17.*t50.*t60+t20.*t51.*t61+t23.*t52.*t62+t26.*t53.*t63+t29.*t54.*t64+t32.*t55.*t65-dphi.*dtheta1.*t4.*t66.*2.0-dphi.*dtheta2.*t7.*t67.*2.0-dphi.*dtheta3.*t10.*t68.*2.0-dphi.*dtheta4.*t13.*t69.*2.0-dphi.*dtheta5.*t16.*t70.*2.0-dphi.*dtheta6.*t19.*t71.*2.0-dphi.*dtheta7.*t22.*t72.*2.0-dphi.*dtheta8.*t25.*t73.*2.0-dphi.*dtheta9.*t28.*t74.*2.0-dphi.*dtheta10.*t31.*t75.*2.0)+t33.*(R1.*t2.*sin(t34)+R2.*t2.*sin(t35)+R3.*t2.*sin(t36)+R4.*t2.*sin(t37)+R5.*t2.*sin(t38)+R6.*t2.*sin(t39)+R7.*t2.*sin(t40)+R8.*t2.*sin(t41)+R9.*t2.*sin(t42)+R10.*t2.*sin(t43));Uphi-t45.*(R1.*t5.*t56.*sin(t77)+R2.*t8.*t57.*sin(t78)+R3.*t11.*t58.*sin(t79)+R4.*t14.*t59.*sin(t80)+R5.*t17.*t60.*sin(t81)+R6.*t20.*t61.*sin(t82)+R7.*t23.*t62.*sin(t83)+R8.*t26.*t63.*sin(t84)+R9.*t29.*t64.*sin(t85)+R10.*t32.*t65.*sin(t86)+R1.*dphi.*dtheta1.*t66.*t89.*2.0+R2.*dphi.*dtheta2.*t67.*t93.*2.0+R3.*dphi.*dtheta3.*t68.*t94.*2.0+R4.*dphi.*dtheta4.*t69.*t95.*2.0+R5.*dphi.*dtheta5.*t70.*t96.*2.0+R6.*dphi.*dtheta6.*t71.*t97.*2.0+R7.*dphi.*dtheta7.*t72.*t98.*2.0+R8.*dphi.*dtheta8.*t73.*t99.*2.0+R9.*dphi.*dtheta9.*t74.*t100.*2.0+R10.*dphi.*dtheta10.*t75.*t101.*2.0)-ct.*dphi-kt.*phi-(mb.*t87.*2.0+mp.*t88.*2.0).*(dphi.*dtheta1.*t5.*t66+dphi.*dtheta2.*t8.*t67+dphi.*dtheta3.*t11.*t68+dphi.*dtheta4.*t14.*t69+dphi.*dtheta5.*t17.*t70+dphi.*dtheta6.*t20.*t71+dphi.*dtheta7.*t23.*t72+dphi.*dtheta8.*t26.*t73+dphi.*dtheta9.*t29.*t74+dphi.*dtheta10.*t32.*t75);Utheta1-d.*dtheta1+g.*t5.*t45+t2.*t5.*t66.*t92+R1.*t2.*t45.*t66.*t89;Utheta2-d.*dtheta2+g.*t8.*t45+t2.*t8.*t67.*t92+R2.*t2.*t45.*t67.*t93;Utheta3-d.*dtheta3+g.*t11.*t45+t2.*t11.*t68.*t92+R3.*t2.*t45.*t68.*t94;Utheta4-d.*dtheta4+g.*t14.*t45+t2.*t14.*t69.*t92+R4.*t2.*t45.*t69.*t95;Utheta5-d.*dtheta5+g.*t17.*t45+t2.*t17.*t70.*t92+R5.*t2.*t45.*t70.*t96;Utheta6-d.*dtheta6+g.*t20.*t45+t2.*t20.*t71.*t92+R6.*t2.*t45.*t71.*t97;Utheta7-d.*dtheta7+g.*t23.*t45+t2.*t23.*t72.*t92+R7.*t2.*t45.*t72.*t98;Utheta8-d.*dtheta8+g.*t26.*t45+t2.*t26.*t73.*t92+R8.*t2.*t45.*t73.*t99;Utheta9-d.*dtheta9+g.*t29.*t45+t2.*t29.*t74.*t92+R9.*t2.*t45.*t74.*t100;Utheta10-d.*dtheta10+g.*t32.*t45+t2.*t32.*t75.*t92+R10.*t2.*t45.*t75.*t101];

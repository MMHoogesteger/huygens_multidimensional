function dxdt = mxyphi_5_row(t,in2,in3,in4,in5)
%MXYPHI_5_ROW
%    DXDT = MXYPHI_5_ROW(T,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Sep-2018 09:50:58

R1 = in4(1,:);
R2 = in4(2,:);
R3 = in4(3,:);
R4 = in4(4,:);
R5 = in4(5,:);
Uphi = in5(3,:);
Utheta1 = in5(4,:);
Utheta2 = in5(5,:);
Utheta3 = in5(6,:);
Utheta4 = in5(7,:);
Utheta5 = in5(8,:);
Ux = in5(1,:);
Uy = in5(2,:);
c = in3(12,:);
ct = in3(13,:);
d = in3(6,:);
dphi = in2(11,:);
dtheta1 = in2(12,:);
dtheta2 = in2(13,:);
dtheta3 = in2(14,:);
dtheta4 = in2(15,:);
dtheta5 = in2(16,:);
dx = in2(9,:);
dy = in2(10,:);
g = in3(5,:);
gamma1 = in4(6,:);
gamma2 = in4(7,:);
gamma3 = in4(8,:);
gamma4 = in4(9,:);
gamma5 = in4(10,:);
k = in3(9,:);
kt = in3(10,:);
lb = in3(4,:);
lp = in3(3,:);
mM = in3(8,:);
mb = in3(2,:);
mp = in3(1,:);
phi = in2(3,:);
psi1 = in4(11,:);
psi2 = in4(12,:);
psi3 = in4(13,:);
psi4 = in4(14,:);
psi5 = in4(15,:);
theta1 = in2(4,:);
theta2 = in2(5,:);
theta3 = in2(6,:);
theta4 = in2(7,:);
theta5 = in2(8,:);
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
t18 = mM+mb+mp;
t19 = gamma1+phi;
t20 = gamma2+phi;
t21 = gamma3+phi;
t22 = gamma4+phi;
t23 = gamma5+phi;
t24 = lb.*mb;
t41 = lp.*mp;
t25 = t24-t41;
t26 = sin(t3);
t27 = sin(t6);
t28 = sin(t9);
t29 = sin(t12);
t30 = sin(t15);
t31 = dtheta1.^2;
t32 = dtheta2.^2;
t33 = dtheta3.^2;
t34 = dtheta4.^2;
t35 = dtheta5.^2;
t36 = cos(theta1);
t37 = cos(theta2);
t38 = cos(theta3);
t39 = cos(theta4);
t40 = cos(theta5);
t42 = gamma1-psi1;
t43 = gamma2-psi2;
t44 = gamma3-psi3;
t45 = gamma4-psi4;
t46 = gamma5-psi5;
t47 = lb.^2;
t48 = lp.^2;
t49 = cos(t42);
t50 = mb.*t47;
t51 = mp.*t48;
t52 = t50+t51;
t53 = cos(t43);
t54 = cos(t44);
t55 = cos(t45);
t56 = cos(t46);
dxdt = [dx;dy;dphi;dtheta1;dtheta2;dtheta3;dtheta4;dtheta5;Ux+t18.*(R1.*t2.*cos(t19)+R2.*t2.*cos(t20)+R3.*t2.*cos(t21)+R4.*t2.*cos(t22)+R5.*t2.*cos(t23))-c.*dx-k.*x+t25.*(t2.*t4.*t5+t2.*t7.*t8+t2.*t10.*t11+t2.*t13.*t14+t2.*t16.*t17+t4.*t5.*t31+t7.*t8.*t32+t10.*t11.*t33+t13.*t14.*t34+t16.*t17.*t35+dphi.*dtheta1.*t26.*t36.*2.0+dphi.*dtheta2.*t27.*t37.*2.0+dphi.*dtheta3.*t28.*t38.*2.0+dphi.*dtheta4.*t29.*t39.*2.0+dphi.*dtheta5.*t30.*t40.*2.0);Uy-c.*dy-k.*y+t18.*(R1.*t2.*sin(t19)+R2.*t2.*sin(t20)+R3.*t2.*sin(t21)+R4.*t2.*sin(t22)+R5.*t2.*sin(t23))+t25.*(t2.*t5.*t26+t2.*t8.*t27+t2.*t11.*t28+t2.*t14.*t29+t2.*t17.*t30+t5.*t26.*t31+t8.*t27.*t32+t11.*t28.*t33+t14.*t29.*t34+t17.*t30.*t35-dphi.*dtheta1.*t4.*t36.*2.0-dphi.*dtheta2.*t7.*t37.*2.0-dphi.*dtheta3.*t10.*t38.*2.0-dphi.*dtheta4.*t13.*t39.*2.0-dphi.*dtheta5.*t16.*t40.*2.0);Uphi-ct.*dphi-kt.*phi-t25.*(R1.*t5.*t31.*sin(t42)+R2.*t8.*t32.*sin(t43)+R3.*t11.*t33.*sin(t44)+R4.*t14.*t34.*sin(t45)+R5.*t17.*t35.*sin(t46)+R1.*dphi.*dtheta1.*t36.*t49.*2.0+R2.*dphi.*dtheta2.*t37.*t53.*2.0+R3.*dphi.*dtheta3.*t38.*t54.*2.0+R4.*dphi.*dtheta4.*t39.*t55.*2.0+R5.*dphi.*dtheta5.*t40.*t56.*2.0)-(mb.*t47.*2.0+mp.*t48.*2.0).*(dphi.*dtheta1.*t5.*t36+dphi.*dtheta2.*t8.*t37+dphi.*dtheta3.*t11.*t38+dphi.*dtheta4.*t14.*t39+dphi.*dtheta5.*t17.*t40);Utheta1-d.*dtheta1+g.*t5.*t25+t2.*t5.*t36.*t52+R1.*t2.*t25.*t36.*t49;Utheta2-d.*dtheta2+g.*t8.*t25+t2.*t8.*t37.*t52+R2.*t2.*t25.*t37.*t53;Utheta3-d.*dtheta3+g.*t11.*t25+t2.*t11.*t38.*t52+R3.*t2.*t25.*t38.*t54;Utheta4-d.*dtheta4+g.*t14.*t25+t2.*t14.*t39.*t52+R4.*t2.*t25.*t39.*t55;Utheta5-d.*dtheta5+g.*t17.*t25+t2.*t17.*t40.*t52+R5.*t2.*t25.*t40.*t56];

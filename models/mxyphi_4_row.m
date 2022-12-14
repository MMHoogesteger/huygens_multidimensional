function dxdt = mxyphi_4_row(t,in2,in3,in4,in5)
%MXYPHI_4_ROW
%    DXDT = MXYPHI_4_ROW(T,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Sep-2018 09:50:20

R1 = in4(1,:);
R2 = in4(2,:);
R3 = in4(3,:);
R4 = in4(4,:);
Uphi = in5(3,:);
Utheta1 = in5(4,:);
Utheta2 = in5(5,:);
Utheta3 = in5(6,:);
Utheta4 = in5(7,:);
Ux = in5(1,:);
Uy = in5(2,:);
c = in3(12,:);
ct = in3(13,:);
d = in3(6,:);
dphi = in2(10,:);
dtheta1 = in2(11,:);
dtheta2 = in2(12,:);
dtheta3 = in2(13,:);
dtheta4 = in2(14,:);
dx = in2(8,:);
dy = in2(9,:);
g = in3(5,:);
gamma1 = in4(5,:);
gamma2 = in4(6,:);
gamma3 = in4(7,:);
gamma4 = in4(8,:);
k = in3(9,:);
kt = in3(10,:);
lb = in3(4,:);
lp = in3(3,:);
mM = in3(8,:);
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
t15 = lb.*mb;
t34 = lp.*mp;
t16 = t15-t34;
t17 = sin(t3);
t18 = sin(t6);
t19 = sin(t9);
t20 = sin(t12);
t21 = dtheta1.^2;
t22 = dtheta2.^2;
t23 = dtheta3.^2;
t24 = dtheta4.^2;
t25 = cos(theta1);
t26 = cos(theta2);
t27 = cos(theta3);
t28 = cos(theta4);
t29 = mM+mb+mp;
t30 = gamma1+phi;
t31 = gamma2+phi;
t32 = gamma3+phi;
t33 = gamma4+phi;
t35 = gamma1-psi1;
t36 = gamma2-psi2;
t37 = gamma3-psi3;
t38 = gamma4-psi4;
t39 = lb.^2;
t40 = lp.^2;
t41 = cos(t35);
t42 = mb.*t39;
t43 = mp.*t40;
t44 = t42+t43;
t45 = cos(t36);
t46 = cos(t37);
t47 = cos(t38);
dxdt = [dx;dy;dphi;dtheta1;dtheta2;dtheta3;dtheta4;Ux+t29.*(R1.*t2.*cos(t30)+R2.*t2.*cos(t31)+R3.*t2.*cos(t32)+R4.*t2.*cos(t33))-c.*dx-k.*x+t16.*(t2.*t4.*t5+t2.*t7.*t8+t2.*t10.*t11+t2.*t13.*t14+t4.*t5.*t21+t7.*t8.*t22+t10.*t11.*t23+t13.*t14.*t24+dphi.*dtheta1.*t17.*t25.*2.0+dphi.*dtheta2.*t18.*t26.*2.0+dphi.*dtheta3.*t19.*t27.*2.0+dphi.*dtheta4.*t20.*t28.*2.0);Uy-c.*dy-k.*y+t16.*(t2.*t5.*t17+t2.*t8.*t18+t2.*t11.*t19+t2.*t14.*t20+t5.*t17.*t21+t8.*t18.*t22+t11.*t19.*t23+t14.*t20.*t24-dphi.*dtheta1.*t4.*t25.*2.0-dphi.*dtheta2.*t7.*t26.*2.0-dphi.*dtheta3.*t10.*t27.*2.0-dphi.*dtheta4.*t13.*t28.*2.0)+t29.*(R1.*t2.*sin(t30)+R2.*t2.*sin(t31)+R3.*t2.*sin(t32)+R4.*t2.*sin(t33));Uphi-ct.*dphi-kt.*phi-t16.*(R1.*t5.*t21.*sin(t35)+R2.*t8.*t22.*sin(t36)+R3.*t11.*t23.*sin(t37)+R4.*t14.*t24.*sin(t38)+R1.*dphi.*dtheta1.*t25.*t41.*2.0+R2.*dphi.*dtheta2.*t26.*t45.*2.0+R3.*dphi.*dtheta3.*t27.*t46.*2.0+R4.*dphi.*dtheta4.*t28.*t47.*2.0)-(mb.*t39.*2.0+mp.*t40.*2.0).*(dphi.*dtheta1.*t5.*t25+dphi.*dtheta2.*t8.*t26+dphi.*dtheta3.*t11.*t27+dphi.*dtheta4.*t14.*t28);Utheta1-d.*dtheta1+g.*t5.*t16+t2.*t5.*t25.*t44+R1.*t2.*t16.*t25.*t41;Utheta2-d.*dtheta2+g.*t8.*t16+t2.*t8.*t26.*t44+R2.*t2.*t16.*t26.*t45;Utheta3-d.*dtheta3+g.*t11.*t16+t2.*t11.*t27.*t44+R3.*t2.*t16.*t27.*t46;Utheta4-d.*dtheta4+g.*t14.*t16+t2.*t14.*t28.*t44+R4.*t2.*t16.*t28.*t47];

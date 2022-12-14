function dxdt = mxyphi_3_row(t,in2,in3,in4,in5)
%MXYPHI_3_ROW
%    DXDT = MXYPHI_3_ROW(T,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Sep-2018 09:49:43

R1 = in4(1,:);
R2 = in4(2,:);
R3 = in4(3,:);
Uphi = in5(3,:);
Utheta1 = in5(4,:);
Utheta2 = in5(5,:);
Utheta3 = in5(6,:);
Ux = in5(1,:);
Uy = in5(2,:);
c = in3(12,:);
ct = in3(13,:);
d = in3(6,:);
dphi = in2(9,:);
dtheta1 = in2(10,:);
dtheta2 = in2(11,:);
dtheta3 = in2(12,:);
dx = in2(7,:);
dy = in2(8,:);
g = in3(5,:);
gamma1 = in4(4,:);
gamma2 = in4(5,:);
gamma3 = in4(6,:);
k = in3(9,:);
kt = in3(10,:);
lb = in3(4,:);
lp = in3(3,:);
mM = in3(8,:);
mb = in3(2,:);
mp = in3(1,:);
phi = in2(3,:);
psi1 = in4(7,:);
psi2 = in4(8,:);
psi3 = in4(9,:);
theta1 = in2(4,:);
theta2 = in2(5,:);
theta3 = in2(6,:);
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
t12 = mM+mb+mp;
t13 = gamma1+phi;
t14 = gamma2+phi;
t15 = gamma3+phi;
t16 = lb.*mb;
t27 = lp.*mp;
t17 = t16-t27;
t18 = sin(t3);
t19 = sin(t6);
t20 = sin(t9);
t21 = dtheta1.^2;
t22 = dtheta2.^2;
t23 = dtheta3.^2;
t24 = cos(theta1);
t25 = cos(theta2);
t26 = cos(theta3);
t28 = gamma1-psi1;
t29 = gamma2-psi2;
t30 = gamma3-psi3;
t31 = lb.^2;
t32 = lp.^2;
t33 = cos(t28);
t34 = mb.*t31;
t35 = mp.*t32;
t36 = t34+t35;
t37 = cos(t29);
t38 = cos(t30);
dxdt = [dx;dy;dphi;dtheta1;dtheta2;dtheta3;Ux-c.*dx+t12.*(R1.*t2.*cos(t13)+R2.*t2.*cos(t14)+R3.*t2.*cos(t15))-k.*x+t17.*(t2.*t4.*t5+t2.*t7.*t8+t2.*t10.*t11+t4.*t5.*t21+t7.*t8.*t22+t10.*t11.*t23+dphi.*dtheta1.*t18.*t24.*2.0+dphi.*dtheta2.*t19.*t25.*2.0+dphi.*dtheta3.*t20.*t26.*2.0);Uy-c.*dy-k.*y+t12.*(R1.*t2.*sin(t13)+R2.*t2.*sin(t14)+R3.*t2.*sin(t15))+t17.*(t2.*t5.*t18+t2.*t8.*t19+t2.*t11.*t20+t5.*t18.*t21+t8.*t19.*t22+t11.*t20.*t23-dphi.*dtheta1.*t4.*t24.*2.0-dphi.*dtheta2.*t7.*t25.*2.0-dphi.*dtheta3.*t10.*t26.*2.0);Uphi-ct.*dphi-kt.*phi-(mb.*t31.*2.0+mp.*t32.*2.0).*(dphi.*dtheta1.*t5.*t24+dphi.*dtheta2.*t8.*t25+dphi.*dtheta3.*t11.*t26)-t17.*(R1.*t5.*t21.*sin(t28)+R2.*t8.*t22.*sin(t29)+R3.*t11.*t23.*sin(t30)+R1.*dphi.*dtheta1.*t24.*t33.*2.0+R2.*dphi.*dtheta2.*t25.*t37.*2.0+R3.*dphi.*dtheta3.*t26.*t38.*2.0);Utheta1-d.*dtheta1+g.*t5.*t17+t2.*t5.*t24.*t36+R1.*t2.*t17.*t24.*t33;Utheta2-d.*dtheta2+g.*t8.*t17+t2.*t8.*t25.*t36+R2.*t2.*t17.*t25.*t37;Utheta3-d.*dtheta3+g.*t11.*t17+t2.*t11.*t26.*t36+R3.*t2.*t17.*t26.*t38];

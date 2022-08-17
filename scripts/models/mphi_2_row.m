function dxdt = mphi_2_row(t,in2,in3,in4,in5)
%MPHI_2_ROW
%    DXDT = MPHI_2_ROW(T,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    09-Aug-2018 13:59:47

R1 = in4(1,:);
R2 = in4(2,:);
Uphi = in5(1,:);
Utheta1 = in5(2,:);
Utheta2 = in5(3,:);
ct = in3(13,:);
d = in3(6,:);
dphi = in2(4,:);
dtheta1 = in2(5,:);
dtheta2 = in2(6,:);
g = in3(5,:);
gamma1 = in4(3,:);
gamma2 = in4(4,:);
kt = in3(10,:);
lb = in3(4,:);
lp = in3(3,:);
mb = in3(2,:);
mp = in3(1,:);
phi = in2(1,:);
psi1 = in4(5,:);
psi2 = in4(6,:);
theta1 = in2(2,:);
theta2 = in2(3,:);
t2 = sin(theta1);
t3 = sin(theta2);
t4 = cos(theta1);
t5 = gamma1-psi1;
t6 = cos(theta2);
t7 = gamma2-psi2;
t8 = lb.*mb;
t14 = lp.*mp;
t9 = t8-t14;
t10 = lb.^2;
t11 = lp.^2;
t12 = dphi.^2;
t13 = cos(t5);
t15 = mb.*t10;
t16 = mp.*t11;
t17 = t15+t16;
t18 = cos(t7);
dxdt = [dphi;dtheta1;dtheta2;Uphi-ct.*dphi-kt.*phi-t9.*(R1.*dtheta1.^2.*t2.*sin(t5)+R2.*dtheta2.^2.*t3.*sin(t7)+R1.*dphi.*dtheta1.*t4.*t13.*2.0+R2.*dphi.*dtheta2.*t6.*t18.*2.0)-(dphi.*dtheta1.*t2.*t4+dphi.*dtheta2.*t3.*t6).*(mb.*t10.*2.0+mp.*t11.*2.0);Utheta1-d.*dtheta1+g.*t2.*t9+t2.*t4.*t12.*t17+R1.*t4.*t9.*t12.*t13;Utheta2-d.*dtheta2+g.*t3.*t9+t3.*t6.*t12.*t17+R2.*t6.*t9.*t12.*t18];

function M = mx_3_mass(t,in2,in3,in4)
%MX_3_MASS
%    M = MX_3_MASS(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    09-Aug-2018 14:00:24

lb = in3(4,:);
lp = in3(3,:);
mM = in3(8,:);
mP = in3(7,:);
mb = in3(2,:);
mp = in3(1,:);
psi1 = in4(7,:);
psi2 = in4(8,:);
psi3 = in4(9,:);
theta1 = in2(2,:);
theta2 = in2(3,:);
theta3 = in2(4,:);
t2 = lb.*mb;
t4 = lp.*mp;
t3 = t2-t4;
t5 = cos(psi1);
t6 = cos(theta1);
t7 = t3.*t5.*t6;
t8 = cos(psi2);
t9 = cos(theta2);
t10 = t3.*t8.*t9;
t11 = lb.^2;
t12 = mb.*t11;
t13 = lp.^2;
t14 = mp.*t13;
t15 = t12+t14;
t16 = cos(psi3);
t17 = cos(theta3);
t18 = t3.*t16.*t17;
M = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,mM.*3.0+mP+mb.*3.0+mp.*3.0,t7,t10,t18,0.0,0.0,0.0,0.0,t7,t15,0.0,0.0,0.0,0.0,0.0,0.0,t10,0.0,t15,0.0,0.0,0.0,0.0,0.0,t18,0.0,0.0,t15],[8,8]);

function M = my_2_mass(t,in2,in3,in4)
%MY_2_MASS
%    M = MY_2_MASS(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    09-Aug-2018 14:00:02

lb = in3(4,:);
lp = in3(3,:);
mM = in3(8,:);
mP = in3(7,:);
mb = in3(2,:);
mp = in3(1,:);
psi1 = in4(5,:);
psi2 = in4(6,:);
theta1 = in2(2,:);
theta2 = in2(3,:);
t2 = lb.*mb;
t6 = lp.*mp;
t3 = t2-t6;
t4 = cos(theta1);
t5 = sin(psi1);
t7 = t3.*t4.*t5;
t8 = cos(theta2);
t9 = sin(psi2);
t10 = t3.*t8.*t9;
t11 = lb.^2;
t12 = mb.*t11;
t13 = lp.^2;
t14 = mp.*t13;
t15 = t12+t14;
M = reshape([1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,0.0,0.0,mM.*2.0+mP+mb.*2.0+mp.*2.0,t7,t10,0.0,0.0,0.0,t7,t15,0.0,0.0,0.0,0.0,t10,0.0,t15],[6,6]);

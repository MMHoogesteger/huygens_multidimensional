function dxdt = mxy_3_row(t,in2,in3,in4,in5)
%MXY_3_ROW
%    DXDT = MXY_3_ROW(T,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    09-Aug-2018 14:00:13

Utheta1 = in5(3,:);
Utheta2 = in5(4,:);
Utheta3 = in5(5,:);
Ux = in5(1,:);
Uy = in5(2,:);
c = in3(12,:);
d = in3(6,:);
dtheta1 = in2(8,:);
dtheta2 = in2(9,:);
dtheta3 = in2(10,:);
dx = in2(6,:);
dy = in2(7,:);
g = in3(5,:);
k = in3(9,:);
lb = in3(4,:);
lp = in3(3,:);
mb = in3(2,:);
mp = in3(1,:);
psi1 = in4(7,:);
psi2 = in4(8,:);
psi3 = in4(9,:);
theta1 = in2(3,:);
theta2 = in2(4,:);
theta3 = in2(5,:);
x = in2(1,:);
y = in2(2,:);
t2 = lb.*mb;
t10 = lp.*mp;
t3 = t2-t10;
t4 = dtheta1.^2;
t5 = sin(theta1);
t6 = dtheta2.^2;
t7 = sin(theta2);
t8 = dtheta3.^2;
t9 = sin(theta3);
dxdt = [dx;dy;dtheta1;dtheta2;dtheta3;Ux-c.*dx-k.*x+t3.*(t4.*t5.*cos(psi1)+t6.*t7.*cos(psi2)+t8.*t9.*cos(psi3));Uy-c.*dy-k.*y+t3.*(t4.*t5.*sin(psi1)+t6.*t7.*sin(psi2)+t8.*t9.*sin(psi3));Utheta1-d.*dtheta1+g.*t3.*t5;Utheta2-d.*dtheta2+g.*t3.*t7;Utheta3-d.*dtheta3+g.*t3.*t9];
function dxdt = mx_2_row(t,in2,in3,in4,in5)
%MX_2_ROW
%    DXDT = MX_2_ROW(T,IN2,IN3,IN4,IN5)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    09-Aug-2018 14:00:00

Utheta1 = in5(2,:);
Utheta2 = in5(3,:);
Ux = in5(1,:);
c = in3(12,:);
d = in3(6,:);
dtheta1 = in2(5,:);
dtheta2 = in2(6,:);
dx = in2(4,:);
g = in3(5,:);
k = in3(9,:);
lb = in3(4,:);
lp = in3(3,:);
mb = in3(2,:);
mp = in3(1,:);
psi1 = in4(5,:);
psi2 = in4(6,:);
theta1 = in2(2,:);
theta2 = in2(3,:);
x = in2(1,:);
t2 = sin(theta1);
t3 = lb.*mb;
t6 = lp.*mp;
t4 = t3-t6;
t5 = sin(theta2);
dxdt = [dx;dtheta1;dtheta2;Ux+t4.*(dtheta1.^2.*t2.*cos(psi1)+dtheta2.^2.*t5.*cos(psi2))-c.*dx-k.*x;Utheta1-d.*dtheta1+g.*t2.*t4;Utheta2-d.*dtheta2+g.*t4.*t5];
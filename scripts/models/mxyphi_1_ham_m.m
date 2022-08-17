function HM = mxyphi_1_ham_m(t,in2,in3,in4)
%MXYPHI_1_HAM_M
%    HM = MXYPHI_1_HAM_M(T,IN2,IN3,IN4)

%    This function was generated by the Symbolic Math Toolbox version 8.1.
%    24-Sep-2018 09:48:37

JP33 = in3(18,:);
dphi = in2(7,:);
dtheta1 = in2(8,:);
dx = in2(5,:);
dy = in2(6,:);
g = in3(5,:);
k = in3(9,:);
kt = in3(10,:);
lb = in3(4,:);
lp = in3(3,:);
mP = in3(7,:);
mb = in3(2,:);
mp = in3(1,:);
phi = in2(3,:);
theta1 = in2(4,:);
x = in2(1,:);
y = in2(2,:);
HM = [dx.^2.*mP.*(1.0./2.0)+k.*x.^2.*(1.0./2.0);dy.^2.*mP.*(1.0./2.0)+k.*y.^2.*(1.0./2.0);JP33.*dphi.^2.*(1.0./2.0)+kt.*phi.^2.*(1.0./2.0);dtheta1.^2.*(lb.^2.*mb.*(1.0./2.0)+lp.^2.*mp.*(1.0./2.0))+g.*cos(theta1).*(lb.*mb-lp.*mp)];
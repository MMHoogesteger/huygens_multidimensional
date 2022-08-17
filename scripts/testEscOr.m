prepEnv
plotInit

%%
t = 0:1e-5:1;
%t = 0.14:1e-7:0.15;

theta = 0.7*sin(1.7*2*pi*t);
dtheta = diff(theta)./(t(2)-t(1));
dtheta = [dtheta(1) dtheta];

ids = 1:round(numel(t)/2000):numel(t);
figure;
subplot(4,1,1)
plot(t(ids),theta(ids),'k')
xlim([min(t) max(t)])
legend('hide')
ylabel('$ \theta$ [rad]','Interpreter','Latex')
% hold on
% plot(t(ids),dtheta(ids),'k--')

p = getparams;
e = p.epsilon;
tau = p.tau;
t_s = p.theta_s;
t_e = p.theta_e;

s = tanh(dtheta./e);
u = tau/2.*s.*(tanh((s.*theta-t_s)/e) - tanh((s.*theta-t_e)/e));
subplot(4,1,2)
plot(t(ids),10*u(ids),'k')
xlim([min(t) max(t)])
ylim([-5 5]*1e-4)
legend('hide')
ylabel('$ u$ [N]','Interpreter','Latex')
subplot(4,1,3)
dsdtd = 1/e.*sech(dtheta./e).^2;
dudtheta = tau/2.*s.*(s/e.*sech((s.*theta-t_s)/e).^2 - s/e.*sech((s.*theta-t_e)/e).^2);
dudthetadot = tau/2.*dsdtd.*(tanh((s.*theta-t_s)/e) - tanh((s.*theta-t_e)/e))...
    + tau/2.*s.*(dsdtd.*sech((s.*theta-t_s)/e).^2 - dsdtd.*sech((s.*theta-t_e)/e).^2);
plot(t(ids),dudtheta(ids),'k')
xlim([min(t) max(t)])
ylim([-5 5]*1e-4)
legend('hide')
ylabel('$\frac{\partial}{\partial \theta} u$ [N/rad]','Interpreter','Latex')
subplot(4,1,4)
plot(t(ids),dudthetadot(ids),'k')
xlim([min(t) max(t)])
ylim([-5 10]*1e-4)
legend('hide')
xlabel('$t$ [s]','Interpreter','Latex')
ylabel('$\frac{\partial}{\partial \dot{\theta}} u$ [Ns/rad]','Interpreter','Latex')
genTikz('exc_or_ders');
%genTikz('exc_or_ders_zoom');
function plotSignals3(saveName,savePlot,t,x,y,phi,theta1,theta2,theta3)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

figure;

subplot(3,1,1)
plot(t,x,'k-')
hold on
plot(t,y,'k--')
ylabel('$x,y$ [mm]','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim([-1 1])
legend('$x$','$y$')

subplot(3,1,2)
plot(t,phi,'k-')
ylabel('$\varphi$ [rad]','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim(1e-2*[-2 2])
legend('$\varphi$')

subplot(3,1,3)
plot(t,theta1,'k-')
hold on
plot(t,theta2,'k--')
plot(t,theta3,'k-.')
xlabel('$t$ [s]','Interpreter','Latex')
ylabel(['$\theta_{i}$ [rad]'],'Interpreter','Latex')
grid on
legend('$\theta_1$','$\theta_2$','$\theta_3$')
xlim([min(t) max(t)])
ylim([-1 1])

drawnow

if(savePlot)
    genTikz(saveName);
end
end


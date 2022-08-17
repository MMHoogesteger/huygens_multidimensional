function plotPhasePortrait3(saveName,savePlot,theta1,theta2,theta3)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
cGray = [0.3 0.3 0.3];
figure;
plot(theta1,theta2,'k-')
hold on
plot(theta1,theta3,'r--')
grid on
xlabel('$\theta_{1}$ [rad]','Interpreter','Latex')
ylabel('$\theta_{i}$ [rad]','Interpreter','Latex')
xlim([-1 1])
ylim([-1 1])
legend('$\theta_2$','$\theta_3$')
drawnow

if(savePlot)
    genTikz(saveName);
end
end


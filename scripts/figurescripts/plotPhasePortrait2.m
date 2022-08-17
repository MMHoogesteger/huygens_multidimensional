function plotPhasePortrait2(saveName,savePlot,theta1,theta2)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

figure;
plot(theta1,theta2,'k-')
grid on
xlabel('$\theta_{1}$ [rad]','Interpreter','Latex')
ylabel('$\theta_{2}$ [rad]','Interpreter','Latex')
xlim([-1 1])
ylim([-1 1])
legend('hide')
drawnow

if(savePlot)
    genTikz(saveName);
end
end


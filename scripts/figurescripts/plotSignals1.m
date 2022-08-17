function plotSignals1(saveName,savePlot,t,x,y,phi,theta,freeyLims)
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
if(freeyLims)
    ylim(ceil(max([max(abs(y)) max(abs(x))]))*[-1 1])
else
    ylim([-1 1])
end
legend('$x$','$y$')

subplot(3,1,2)
plot(t,phi,'k-')
ylabel('$\varphi$ [rad]','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
% if(freeyLims)
%     ylim(ceil(max(abs(phi)))*[-1 1])
% else
    ylim(1e-2*[-3 3])
% end
legend('$\varphi$')

subplot(3,1,3)
plot(t,theta,'k-')
xlabel('$t$ [s]','Interpreter','Latex')
ylabel(['$\theta$ [rad]'],'Interpreter','Latex')
grid on
legend('$\theta_1$')
xlim([min(t) max(t)])
if(freeyLims)
    ylim(ceil(max(abs(theta)))*[-1 1])
else
    ylim([-1 1])
end

drawnow

if(savePlot)
    genTikz(saveName);
end
end


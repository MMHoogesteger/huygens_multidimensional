%% Load data
loadTang3Vids
plotInit
savePlot = true;
%% Plot data
nVec = [1 2];
sVec = {'v_tang3_lin','v_tang3_rot'};
corF = [1 1];
for n = 1:numel(nVec)
    saveName = sVec{n};
    r = rv(nVec(n));

    % Plot signals of the first 50 seconds
    tS = 30;
    idxs = r.signals.t<tS;
    t = r.signals.t(idxs);

    x = r.signals.x_cor(idxs);
    y = -r.signals.y_cor(idxs);
    phi = r.signals.phi_cor(idxs);

    theta1 = corF(n)*r.signals.met(1).theta_cor(idxs);
    theta2 = corF(n)*r.signals.met(2).theta_cor(idxs);
    theta3 = corF(n)*r.signals.met(3).theta_cor(idxs);

    plotSignals3([saveName '_t_q'],savePlot,t,x,y,phi,theta1,theta2,theta3)

    % Plot signals SS
    tL = 3;
    idxs = r.signals.t>(max(r.signals.t)-tL);
    t = r.signals.t(idxs);

    x = r.signals.x_cor(idxs);
    y = -r.signals.y_cor(idxs);
    phi = r.signals.phi_cor(idxs);

    theta1 = corF(n)*r.signals.met(1).theta_cor(idxs);
    theta2 = corF(n)*r.signals.met(2).theta_cor(idxs);
    theta3 = corF(n)*r.signals.met(3).theta_cor(idxs);

    plotSignals3([saveName '_t_q_ss'],savePlot,t,x,y,phi,theta1,theta2,theta3)
    
    % Plot phase portrait
    plotPhasePortrait3([saveName '_phaseportrait_ss'],savePlot,theta1,theta2,theta3)
    
    
end


%% Plot platform phase signals

t = r.signals.t;

figure;
subplot(3,1,1)
plot(t,wrapTo2Pi(r.analysis.dphase_x),'Color',cBlack)
ylabel('$p_{x\theta_1} [rad]$','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
%ylim([-1 1])
legend('hide')

subplot(3,1,2)
plot(t,wrapTo2Pi(r.analysis.dphase_y),'Color',cBlack)
ylabel('$p_{y\theta_1} [rad]$','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
%ylim([-1 1])
legend('hide')

subplot(3,1,3)
plot(t,wrapTo2Pi(r.analysis.dphase_phi),'Color',cBlack)
xlabel('$t [s]$','Interpreter','Latex')
ylabel('$p_{\varphi\theta_1} [rad]$','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
%ylim([-0.25 0.25])
legend('hide')
drawnow

if(savePlot)
    matlab2tikz([saveName '_p_xyphi.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end


%% Plot platform freq signals

t = r.signals.t;

figure;
subplot(3,1,1)
plot(t(1:end-1),smooth(r.analysis.freq_x,70),'Color',cBlack)
ylabel('$x [Hz]$','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim([0 2.5])
legend('hide')

subplot(3,1,2)
plot(t(1:end-1),smooth(r.analysis.freq_y,70),'Color',cBlack)
ylabel('$y [Hz]$','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim([0 2.5])
legend('hide')

subplot(3,1,3)
plot(t(1:end-1),smooth(r.analysis.freq_phi,70),'Color',cBlack)
xlabel('$t [s]$','Interpreter','Latex')
ylabel('$\varphi [Hz]$','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim([0 2.5])
legend('hide')
drawnow

if(savePlot)
    matlab2tikz([saveName '_f_xyphi.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end

%% Plot metronome phase signals

t = r.signals.t;

nMets= numel(r.signals.met);
figure;
for nMet = 2:nMets
    subplot(nMets-1,1,nMet-1)
    plot(t,wrapToPi(r.analysis.dphase_mets(:,nMet)),'Color',cBlack)
    xlabel('$t [s]$','Interpreter','Latex')
    ylabel(['$p_{\theta_{' num2str(nMet) '}\theta_1} [rad]$'],'Interpreter','Latex')
    grid on
    xlim([min(t) max(t)])
    %ylim([-3 3])
    legend('hide')
end
drawnow

if(savePlot)
    matlab2tikz([saveName '_p_theta.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end


%% Plot metronome freq signals
t = r.signals.t;

nMets= numel(r.signals.met);
figure;
for nMet = 1:nMets
    subplot(nMets,1,nMet)
    plot(t(1:end-1),smooth(r.analysis.freq_mets(:,nMet),70),'Color',cBlack)
    xlabel('$t [s]$','Interpreter','Latex')
    ylabel(['$F_{\theta_{' num2str(nMet) '}} [mm]$'],'Interpreter','Latex')
    grid on
    xlim([0 max(t)])
    ylim([1.7 1.75])
    legend('hide')
end
drawnow

if(savePlot)
    matlab2tikz([saveName '_f_theta.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end

%% Plot platform fourier


figure;
semilogy(0:6,abs([r.analysis.fourier.Fz0(1) r.analysis.fourier.Fzkn(1,:)]),'o','Color',cBlack)
hold on
semilogy(0:6,abs([r.analysis.fourier.Fz0(2) r.analysis.fourier.Fzkn(2,:)]),'s','Color',cBlack)
semilogy(0:6,abs([r.analysis.fourier.Fz0(3) r.analysis.fourier.Fzkn(3,:)]),'d','Color',cBlack)
semilogy(0:6,abs([r.analysis.fourier.Fz0(4) r.analysis.fourier.Fzkn(4,:)]),'x','Color',cBlack)

xlabel('$k [-]$','Interpreter','Latex')
ylabel('$|a_k| [mm]$','Interpreter','Latex')
grid on
xlim([-0.5 6.5])
%ylim([0 2.5])

legend('x','y','\varphi','\theta')


if(savePlot)
    matlab2tikz([saveName '_fourier.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end

%% Plot platform fourier error


figure;
semilogy(1:6,abs([r.analysis.fourier.Ezk(1,:)]),'o','Color',cBlack)
hold on
semilogy(1:6,abs([r.analysis.fourier.Ezk(2,:)]),'s','Color',cBlack)
semilogy(1:6,abs([r.analysis.fourier.Ezk(3,:)]),'d','Color',cBlack)
semilogy(1:6,abs([r.analysis.fourier.Ezk(4,:)]),'x','Color',cBlack)

xlabel('$k [-]$','Interpreter','Latex')
ylabel('$|a_k| [mm]$','Interpreter','Latex')
grid on
xlim([0.5 6.5])
%ylim([0 2.5])

legend('x','y','\varphi','\theta')


if(savePlot)
    matlab2tikz([saveName '_fourier_error.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end


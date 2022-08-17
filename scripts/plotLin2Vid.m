%% Load data
load('H:\Vids\Vids_huygens\video2_results.mat')
plotInit
savePlot = true;
%% Plot data
saveName = 'v_double';
r = results;

% Plot signals of the first 150 seconds
tS = 0;
tE = 150;
idxs = r.signals.t<tE & r.signals.t>tS;
t = r.signals.t(idxs);

x = r.signals.x_cor(idxs);
y = r.signals.y_cor(idxs);
phi = r.signals.phi_cor(idxs);

theta1 = r.signals.met(1).theta_cor(idxs);
theta2 = r.signals.met(2).theta_cor(idxs);

plotSignals2([saveName '_t_q'],savePlot,t,x,y,phi,theta1,theta2)

% Plot signals SS
tL = 3;
idxs = r.signals.t>(max(r.signals.t)-tL);
t = r.signals.t(idxs);

x = r.signals.x_cor(idxs);
y = r.signals.y_cor(idxs);
phi = r.signals.phi_cor(idxs);

theta1 = -r.signals.met(1).theta_cor(idxs);
theta2 = -r.signals.met(2).theta_cor(idxs);

plotSignals2([saveName '_t_q_ss'],savePlot,t,x,y,phi,theta1,theta2)





%% Plot platform phase signals

t = r.signals.t;

figure;
subplot(3,1,1)
plot(t,wrapToPi(r.analysis.dphase_x-pi),'Color',cBlack)
ylabel('$p_{x\theta_1}$ [rad]','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
%ylim([-1 1])
legend('hide')

subplot(3,1,2)
plot(t,wrapToPi(r.analysis.dphase_y-pi),'Color',cBlack)
ylabel('$p_{y\theta_1}$ [rad]','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
%ylim([-1 1])
legend('hide')

subplot(3,1,3)
plot(t,wrapTo2Pi(r.analysis.dphase_phi-pi),'Color',cBlack)
xlabel('$t$ [s]','Interpreter','Latex')
ylabel('$p_{\varphi\theta_1}$ [rad]','Interpreter','Latex')
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
ylabel('$x$ [Hz]','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim([0 2.5])
legend('hide')

subplot(3,1,2)
plot(t(1:end-1),smooth(r.analysis.freq_y,70),'Color',cBlack)
ylabel('$y$ [Hz]','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim([0 2.5])
legend('hide')

subplot(3,1,3)
plot(t(1:end-1),smooth(r.analysis.freq_phi,70),'Color',cBlack)
xlabel('$t$ [s]','Interpreter','Latex')
ylabel('$\varphi$ [Hz]','Interpreter','Latex')
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


%% Plot pwelch signals

t = r.signals.t;
x = r.signals.x_cor;
theta = r.signals.met.theta_cor;

Fs = 60;
TSS = 100;
NOverlap = 10*Fs;
W = hanning(2*NOverlap);
figure;
subplot(2,1,1)
[Pxx,F] = pwelch(x,W,NOverlap,[],Fs);
plot(F,mag2db(Pxx),'Color',cBlack);
hold on
[Pxx,F] = pwelch(x(t>TSS),W,NOverlap,[],Fs);
plot(F,mag2db(Pxx),'Color',cRed);
ylabel('PSD $x$ [dB/Hz]','Interpreter','Latex')
grid on
xlim([0 15])
%ylim([-200 10])
legend('$t>0$','$t>100$')

subplot(2,1,2)
[Pxx,F] = pwelch(theta,W,NOverlap,[],Fs);
plot(F,mag2db(Pxx),'Color',cBlack);
hold on
[Pxx,F] = pwelch(theta(t>TSS),W,NOverlap,[],Fs);
plot(F,mag2db(Pxx),'Color',cRed);
grid on
xlim([0 15])
%ylim([-200 10])
ylabel('PSD $\theta$ [dB/Hz]','Interpreter','Latex')
xlabel('$F$ [Hz]','Interpreter','Latex')
legend('hide')

if(savePlot)
    matlab2tikz([saveName '_pwelch_xtheta.tex'],'parseStrings',false,...
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

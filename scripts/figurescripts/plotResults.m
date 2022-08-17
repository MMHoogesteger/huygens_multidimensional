r = Results;
saveplot = true;
plotInit
saveName = 'three_rot';
%% Plot platform signals

t = r.signals.t;
x = r.signals.x_cor;
y = r.signals.y_cor;
phi = r.signals.phi_cor;

figure;
subplot(3,1,1)
plot(t,x,'Color',cBlack)
ylabel('$x [mm]$','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim([-3 3])
legend('hide')

subplot(3,1,2)
plot(t,y,'Color',cBlack)
ylabel('$y [mm]$','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim([-3 3])
legend('hide')

subplot(3,1,3)
plot(t,rad2deg(phi),'Color',cBlack)
xlabel('$t [s]$','Interpreter','Latex')
ylabel('$\varphi [degree]$','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim([-0.25 0.25])
legend('hide')
drawnow

genTikz(saveplot,[saveName '_t_xy']);

%% Plot metronome signals

t = r.signals.t;
nMets= numel(r.signals.met);
figure;
for nMet = 1:nMets
subplot(nMets,1,nMet)
plot(t,r.signals.met(nMet).theta_cor,'Color',cBlack)
xlabel('$t [s]$','Interpreter','Latex')
ylabel(['$\theta_{' num2str(nMet) '} [mm]$'],'Interpreter','Latex')
grid on
xlim([min(t) max(t)])
%ylim([-3 3])
legend('hide')
end

drawnow

if(saveplot)
    matlab2tikz([saveName '_theta.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end

%% Plot platform signals SS

t = r.signals.t(r.analysis.idSS);
x = r.signals.x_cor(r.analysis.idSS);
y = r.signals.y_cor(r.analysis.idSS);
phi = r.signals.phi_cor(r.analysis.idSS);

figure;
subplot(3,1,1)
plot(t,x,'Color',cBlack)
ylabel('$x [mm]$','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim([-1 1])
legend('hide')

subplot(3,1,2)
plot(t,y,'Color',cBlack)
ylabel('$y [mm]$','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim([-1 1])
legend('hide')

subplot(3,1,3)
plot(t,rad2deg(phi),'Color',cBlack)
xlabel('$t [s]$','Interpreter','Latex')
ylabel('$\varphi [degree]$','Interpreter','Latex')
grid on
xlim([min(t) max(t)])
ylim([-1.5 1.5])
legend('hide')
drawnow

if(saveplot)
    matlab2tikz([saveName '_t_xy_ss.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end

%% Plot metronome signals SS

t = r.signals.t(r.analysis.idSS);
nMets= numel(r.signals.met);
figure;
for nMet = 1:nMets
    subplot(nMets,1,nMet)
    plot(t,r.signals.met(nMet).theta_cor(r.analysis.idSS),'Color',cBlack)
    xlabel('$t [s]$','Interpreter','Latex')
    ylabel(['$\theta_{' num2str(nMet) '} [rad]$'],'Interpreter','Latex')
    grid on
    xlim([min(t) max(t)])
    %ylim([-3 3])
    legend('hide')
end

drawnow

if(saveplot)
    matlab2tikz([saveName '_theta_ss.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end


%% Plot metronome phase portrait (SS)

%figure;
ids = r.analysis.idSS;
%ids = 1:length(r.signals.t);
t1 = r.signals.met(1).theta_cor(ids);
t2 = r.signals.met(2).theta_cor(ids);
t3 = r.signals.met(3).theta_cor(ids);
scatter3(t1,t2,t3)
xlabel(['$\theta_{1} [rad]$'],'Interpreter','Latex')
ylabel(['$\theta_{2} [rad]$'],'Interpreter','Latex')
zlabel(['$\theta_{3} [rad]$'],'Interpreter','Latex')
grid on

%ylim([-3 3])
legend('hide')

drawnow

% if(saveplot)
%     matlab2tikz([saveName '_theta_ss.tex'],'parseStrings',false,...
%         'height','\figureheight',...
%         'width','\figurewidth',...
%         'showInfo', false);
% end

%% Plot metronome phase portrait (SS)

figure;
ids = r.analysis.idSS;
%ids = 1:length(r.signals.t);
t1 = r.signals.met(1).theta_cor(ids);
t2 = r.signals.met(2).theta_cor(ids);
t3 = r.signals.met(3).theta_cor(ids);
plot(t1,t2,'-','Color',cBlack)
hold on
%plot(t2,t3,'x')
%plot(t3,t1,'x')
plot(t1,t3,'-','Color',cRed)

xlabel(['$\theta_{1} [rad]$'],'Interpreter','Latex')
ylabel(['$\theta_{j} [rad]$'],'Interpreter','Latex')
%zlabel(['$\theta_{3} [rad]$'],'Interpreter','Latex')
grid on

%ylim([-3 3])
legend('j=2','j=3','Location','NorthWest')

drawnow

if(saveplot)
    matlab2tikz([saveName '_theta_phaseportrait.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
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

if(saveplot)
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

if(saveplot)
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

if(saveplot)
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

if(saveplot)
    matlab2tikz([saveName '_f_theta.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end


%% Plot platform x pwelch signals

t = r.signals.t;

figure;
[Pxx,F] = pwelch(r.signals.y_cor,[],[],[],60);
plot(F,mag2db(Pxx),'Color',cBlack);
plot(F,mag2db(Pxx),'Color',cBlack);
xlabel('$F [\mathrm{Hz}]$','Interpreter','Latex')
ylabel('$Power /frequency [\mathrm{dB/Hz}]$','Interpreter','Latex')
grid on
xlim([0 15])
%ylim([0 2.5])

hold on
[Pxx,F] = pwelch(r.signals.y_cor(r.signals.t>60),[],[],[],60);
plot(F,mag2db(Pxx),'Color',cRed);
legend('$t>0$','$t>100$')


if(saveplot)
    matlab2tikz([saveName '_pwelch_x.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end

%% Plot metronome pwelch signals

t = r.signals.t;
metId = 1;
figure;
[Pxx,F] = pwelch(r.signals.met(metId).theta_cor,[],[],[],60);
plot(F,mag2db(Pxx),'Color',cBlack);
plot(F,mag2db(Pxx),'Color',cBlack);
xlabel('$F [\mathrm{Hz}]$','Interpreter','Latex')
ylabel('$Power /frequency [\mathrm{dB/Hz}]$','Interpreter','Latex')
grid on
xlim([0 15])
%ylim([0 2.5])

hold on
[Pxx,F] = pwelch(r.signals.met(metId).theta_cor(r.signals.t>70),[],[],[],60);
plot(F,mag2db(Pxx),'Color',cRed);
legend('$t>0$','$t>100$')


if(saveplot)
    matlab2tikz([saveName '_pwelch_theta.tex'],'parseStrings',false,...
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


if(saveplot)
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


if(saveplot)
    matlab2tikz([saveName '_fourier_error.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end



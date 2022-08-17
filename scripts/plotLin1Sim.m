%% Plot settings and initialization

dataFolder = 'H:\afstudeerdata\sim_1_xyphi_ref\';
load([dataFolder '1xyphi_original'])
plotInit
savePlot = 1;

%% Plot data 
saveName = 's_single';
r = s.results.solution;

plotFs = 60;
% Plot signals of the first 150 seconds
tS = 0;
tE = 150;
idxs = find(r.t<tE & r.t>tS);
idxs = idxs(1:round(s.fs/plotFs):end);

t = r.t(idxs);
x = 1e3*r.x(idxs);
y = 1e3*r.y(idxs);
phi = r.phi(idxs);
theta = r.theta(idxs);

plotSignals1([saveName '_t_q'],savePlot,t,x,y,phi,theta,true)

% Plot signals SS
tL = 3;
idxs = find(r.t>(max(r.t)-tL));
idxs = idxs(1:round(s.fs/plotFs):end);

t = r.t(idxs);
x = 1e3*r.x(idxs);
y = 1e3*r.y(idxs);
phi = r.phi(idxs);
theta = r.theta(idxs);

plotSignals1([saveName '_t_q_ss'],savePlot,t,x,y,phi,theta,false)

%% Plot pwelch signals

t = r.t;
x = 1e3*r.x;
theta = r.theta;

Fs = 1e3;
TSS = 900;
NOverlap = 10*Fs;
W = hanning(2*NOverlap);
figure;
subplot(2,1,1)
[Pxx,F] = pwelch(x,W,NOverlap,[],Fs);
plot(F(F<15.1),mag2db(Pxx(F<15.1)),'Color',cBlack);
hold on
[Pxx,F] = pwelch(x(t>TSS),W,NOverlap,[],Fs);
plot(F(F<15.1),mag2db(Pxx(F<15.1)),'Color',cRed);
ylabel('PSD $x$ [dB/Hz]','Interpreter','Latex')
grid on
xlim([0 15])
%ylim([-250 10])
legend('$t>0$','$t>900$')

subplot(2,1,2)
[Pxx,F] = pwelch(theta,W,NOverlap,[],Fs);
plot(F(F<15.1),mag2db(Pxx(F<15.1)),'Color',cBlack);
hold on
[Pxx,F] = pwelch(theta(t>TSS),W,NOverlap,[],Fs);
plot(F(F<15.1),mag2db(Pxx(F<15.1)),'Color',cRed);
grid on
xlim([0 15])
%ylim([-250 10])

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
semilogy(0:6,1e3*abs([s.results.analysis.fourier.Fz0(1) s.results.analysis.fourier.Fzkn(1,1:6)]),'o','Color',cBlack)
hold on
semilogy(0:6,abs([s.results.analysis.fourier.Fz0(2) s.results.analysis.fourier.Fzkn(2,1:6)]),'s','Color',cBlack)
% semilogy(0:6,abs([s.results.analysis.fourier.Fz0(3) s.results.analysis.fourier.Fzkn(3,1:6)]),'d:','Color',cBlack)
% semilogy(0:6,abs([s.results.analysis.fourier.Fz0(4) s.results.analysis.fourier.Fzkn(4,1:6)]),'x-.','Color',cBlack)

xlabel('$k [-]$','Interpreter','Latex')
ylabel('$|a_k| [mm]$','Interpreter','Latex')
grid on
xlim([-0.5 6.5])
%ylim([0 2.5])

legend('x','\theta')


genTikz(savePlot,[saveName '_fourier']);

%% Plot platform fourier error


figure;
semilogy(1:6,abs([s.results.analysis.fourier.Ezk(1,1:6)]),'o','Color',cBlack)
hold on
semilogy(1:6,abs([s.results.analysis.fourier.Ezk(2,1:6)]),'s','Color',cBlack)
% semilogy(1:6,abs([s.results.analysis.fourier.Ezk(3,1:6)]),'d','Color',cBlack)
% semilogy(1:6,abs([s.results.analysis.fourier.Ezk(4,1:6)]),'x','Color',cBlack)

xlabel('$k [-]$','Interpreter','Latex')
ylabel('$|a_k| [mm]$','Interpreter','Latex')
grid on
xlim([0.5 6.5])
%ylim([0 2.5])

legend('x','\theta')


if(savePlot)
    matlab2tikz([saveName '_fourier_error.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
end

genTikz(savePlot,[saveName '_fourier_error']);
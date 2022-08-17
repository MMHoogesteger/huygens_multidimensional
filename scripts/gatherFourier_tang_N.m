Nmetsv = [3:10];
for nMid = 1:numel(Nmetsv)
    Nmets = Nmetsv(nMid)
    prepEnv
    dataFolder = ['H:\afstudeerdata\new\per_xyphi_' num2str(Nmets) '_R\'];
    plotInit
    saveplot = 0;
    saveName = ['per_xyphi_' num2str(Nmets) '_or'];

    Nsims = 6;
    F = zeros(3+Nmets,12,Nsims);
    for simN = 1:Nsims
        simN
        load([dataFolder 'sim_' num2str(simN) '.mat'])
        for k = 1:12
            F(:,k,simN) = s.results(k).analysis.fourier.Fzkn(1:(3+Nmets),1);
        end
    end

    load([dataFolder 'simSet.mat'])
    R = cell2mat(simSet.simSettings.varVecCell{1});
    R = R(1,:);
    save([dataFolder  'gatherFourier'],'F','R');
end    
%%
Nmets = 4
dataFolder = ['H:\afstudeerdata\new\per_xyphi_' num2str(Nmets) '_R\'];
load([dataFolder 'gatherFourier']);
figure
for n = 1:Nmets
    plot(wrapTo2Pi(squeeze(angle(F(3+n,:,1)))),'ks-');
    hold on
end
%%

figure;
A =abs(F);
hold off
hy = plot(R,1e3*squeeze(A(2,:,:)),'s','Color',cBlack);
hold on
hp = plot(R,rad2deg(squeeze(A(3,:,:))),'d','Color',cBlack);
ht = plot(R,squeeze(A(5,:,:)),'+','Color',cBlack);
grid on
legend([hy(1) hp(2) ht(3)],'$y$','$\varphi$','$\theta_2$')
%xlim([0.5 10.5]*10)
xlabel('$R [mm]$','Interpreter','Latex')
ylabel('$| a_k| [mm]$','Interpreter','Latex')
genTikz(saveplot,[saveName '_amp']);
%%
P = (angle(F)+0.5*pi);

figure;
hold off
hy = plot(R(1:16),squeeze(P(2,:,1:16)),'s','Color',cBlack);
hold on
hp = plot(R(17:end),squeeze(P(3,:,17:end)),'d','Color',cBlack);
ht = plot(R,squeeze(P(5,:,:)),'+','Color',cBlack);
grid on
legend([hy(1) hp(2) ht(3)],'$y$','$\varphi$','$\theta_2$')
%xlim([0.5 10.5]*10)
xlabel('$R [mm]$','Interpreter','Latex')
ylabel('$\angle a_k [rad]$','Interpreter','Latex')
genTikz(saveplot,[saveName '_phase']);

prepEnv
dataFolder = 'D:\afstudeerdata\simulations\tang_xyphi_3_R\';
plotInit
saveplot = 1;
saveName = 'tang_xyphi_3_or';
%%

F = zeros(6,12,25);
F3 = zeros(6,12,25);
O = zeros(1,12,25);
for simN = 1:25
    simN
    load([dataFolder 'sim_' num2str(simN) '.mat'])
    for k = 1:12
        F(:,k,simN) = s.results(k).analysis.fourier.Fzkn(1:6,1);
        F3(:,k,simN) = s.results(k).analysis.fourier.Fzkn(1:6,3);
        O(:,k,simN) = s.results(k).analysis.fourier.omega;
    end
end

load([dataFolder 'simSet.mat'])
R = cell2mat(simSet.simSettings.varVecCell{1});
R = R(1,:);
save([dataFolder  'gatherFourier'],'F','F3','R','O');
    
%%
load([dataFolder 'gatherFourier']);
R = R*1e3;
A =abs(F);
P = (angle(F)+0.5*pi);
%%
figure;

hold off
hy = plot(R,1e3*squeeze(A(2,:,:)),'s','Color',cBlack);
hold on
hp = plot(R,rad2deg(squeeze(A(3,:,:))),'d','Color',cBlack);
ht1 = plot(R,squeeze(A(5,:,:)),'+','Color',cBlack);
ht2 = plot(R,squeeze(A(6,:,:)),'x','Color',cBlack);
grid on
legend([hy(1) hp(1) ht1(1) ht2(1)],'$y$','$\varphi$','$\theta_2$','$\theta_3$','Location','SouthWest')
%xlim([0.5 10.5]*10)
xlabel('$R [mm]$','Interpreter','Latex')
ylabel('$| a_k| [mm]$','Interpreter','Latex')
genTikz(saveplot,[saveName '_amp']);
%%


figure;
hold off
hy = plot(R(1:16),squeeze(P(2,:,1:16)),'s','Color',cBlack);
hold on
hp = plot(R(17:end),squeeze(P(3,:,17:end)),'d','Color',cBlack);
ht = plot(R,squeeze(P(5,:,:)),'+','Color',cBlack);
ht = plot(R,squeeze(P(6,:,:)),'x','Color',cBlack);
grid on
legend([hy(1) hp(1) ht1(1) ht2(1)],'$y$','$\varphi$','$\theta_2$','$\theta_3$','Location','SouthWest')
%xlim([0.5 10.5]*10)
xlabel('$R [mm]$','Interpreter','Latex')
ylabel('$\angle a_k [rad]$','Interpreter','Latex')
genTikz(saveplot,[saveName '_phase']);

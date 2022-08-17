dataFolder = '..\data\linper_xyphi_2_R\';
plotInit
saveplot = 1;
saveName = 'linper_2x_or';
%%

F = zeros(5,12,25);
for simN = 1:25
    simN
    load([dataFolder 'sim_' num2str(simN) '.mat'])
    for k = 1:12
        F(:,k,simN) = s.results(k).analysis.fourier.Fzkn(1:5,1);
    end
end

load([dataFolder 'simSet.mat'])
R = cell2mat(simSet.simSettings.varVecCell{1});
R = R(1,:);
save([dataFolder saveName '_gatherFourier'],'F','R');
    
%%
load([dataFolder saveName '_gatherFourier']);

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

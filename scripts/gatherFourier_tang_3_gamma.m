prepEnv
dataFolder = 'H:\afstudeerdata\tang_xyphi_3_gamma\';
plotInit
saveplot = 0;
saveName = 'tang_xyphi_3_or_gamma';
%%
simNT = 78
F = zeros(6,12,simNT);
for simN = 1:simNT
    simN
    load([dataFolder 'sim_' num2str(simN) '.mat'])
    for k = 1:12
        F(:,k,simN) = s.results(k).analysis.fourier.Fzkn(1:6,1);
    end
end

load([dataFolder 'simSet.mat'])
G = cell2mat(simSet.simSettings.varVecCell{1});
save([dataFolder  'gatherFourier'],'F','G');
    
%%
load([dataFolder 'gatherFourier']);
G = rad2deg(G);

figure;
A =abs(F);
P = rad2deg((angle(F)+0.5*pi));
P6 = squeeze(P(6,:,:));
idsTot = zeros(size(P6));

P6_id = abs(P6)<rad2deg(2e-1);
idsTot(P6_id) = 1;
P6_t= nan(size(P6));
P6_t(P6_id) = P6(P6_id);
hold on
for nRep = 1:12
hy = plot3(G(2,:),G(3,:),squeeze(P6_t(nRep,:)),'x','Color',cBlack);
end

P6_id = abs(P6-180)<rad2deg(2e-1);
idsTot(P6_id) = 1;
P6_t = nan(size(P6));
P6_t(P6_id) = P6(P6_id);
hold on
for nRep = 1:12
hy = plot3(G(2,:),G(3,:),squeeze(P6_t(nRep,:)),'+','Color',cBlack);
end

specialIds = [25 70 46 31];
P6_id = zeros(size(P6));
P6_id(:,specialIds) = 1;
P6_id = logical(P6_id);
idsTot(P6_id) = 1;
P6_t = nan(size(P6));
P6_t(P6_id) = P6(P6_id);
hold on
for nRep = 1:12
hy = plot3(G(2,:),G(3,:),squeeze(P6_t(nRep,:)),'s','Color',cRed);
end

P6_id = idsTot==0;
idsTot(P6_id) = 1;
P6_t = nan(size(P6));
P6_t(P6_id) = P6(P6_id);
hold on
for nRep = 1:12
hy = plot3(G(2,:),G(3,:),squeeze(P6_t(nRep,:)),'.','Color',cGray);
end

grid on
xlabel('$\gamma_2$','Interpreter','Latex')
ylabel('$\gamma_3$','Interpreter','Latex')

%%
load([dataFolder 'gatherFourier']);
%%
G1 = G(2,:);
G2 = G(3,:);

g1 = unique(G1);
g2 = unique(G2);


figure;
A =abs(F);
hold on
for nG = 1:numel(g1)
    xMat = repmat(G1(G1==g1(nG)),1,12);
    yMat = repmat(G2(G1==g1(nG)),1,12);
    hy = plot3(xMat,yMat,reshape(squeeze(P(6,:,G1==g1(nG))).',1,[]),'x');
end
grid on
xlabel('$\gamma_2$','Interpreter','Latex')
ylabel('$\gamma_3$','Interpreter','Latex')
%%
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


figure;
hold off
hy = scatter3(G(2,:),G(3,:),squeeze(P(5,1,:)),'s');
hold on
hp = plot(R(17:end),squeeze(P(3,:,17:end)),'d','Color',cBlack);
ht = plot(R,squeeze(P(5,:,:)),'+','Color',cBlack);
grid on
legend([hy(1) hp(2) ht(3)],'$y$','$\varphi$','$\theta_2$')
%xlim([0.5 10.5]*10)
xlabel('$R [mm]$','Interpreter','Latex')
ylabel('$\angle a_k [rad]$','Interpreter','Latex')
genTikz(saveplot,[saveName '_phase']);

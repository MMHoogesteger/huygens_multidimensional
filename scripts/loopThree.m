res = 4*10^3;
Rinc = 0.002;
%%
disp('Preparation of dist-phase R')

h = genTangOptim(3,Rinc,'h3_3_dist');
optimSet3_3_dist(h);

disp(['Number of unknowns: ' num2str(h.no_T)])

h.getTransferFunctions;
[J,Ai,Pi,Oi] = h.evalError(h.rho_i,res);
disp('init optim')
h.runOptimization(res,false);
rho_last = h.rho_o;
[J,Ai,Pi,Oi] = h.evalError(h.rho_o,res);
disp(['initial residual: ' num2str(J)])
%%
Rvec = Rinc:Rinc:0.5;

nv = length(Rvec);

Jvec = nan(1,1,nv);
Avec = nan([size(Ai) nv]);
Pvec = nan([size(Pi) nv]);
Ovec = nan([size(Oi) nv]);



for n = 1:nv
    disp(['init optim ' num2str(n) ' of ' num2str(nv)])
    R = Rvec(n);
    
    h = genTangOptim(3,R,'h3_3_dist');
    optimSet3_3_dist(h);

    h.getTransferFunctions;
    h.setRhoInit(rho_last);
    
    h.runOptimization(res,false);
%     disp('ref optim')
%     h.continueOptimization(10^4,false);
    [J,Ai,Pi,Oi] = h.evalError(h.rho_o,res);
    disp(['optim ended in ' num2str(J)])
    Jvec(:,:,n) = J;
    Avec(:,:,n) = Ai;
    Pvec(:,:,n) = Pi;
    Ovec(:,:,n) = Oi;
    rho_last = h.rho_o;
    % Generate transfer function and derivative function 
    
end

%%
R3_3_dist_cont.Rvec = Rvec;
R3_3_dist_cont.Jvec = Jvec;
R3_3_dist_cont.Avec = Avec;
R3_3_dist_cont.Pvec = Pvec;
R3_3_dist_cont.Ovec = Ovec;

save('R3_3_dist_cont','R3_3_dist_cont')


%%
disp('Preparation of in-phase R')

h = genTangOptim(3,0.5,'h3_3_in');
optimSet3_3_in(h);
disp(['Number of unknowns: ' num2str(h.no_T)])

h.getTransferFunctions;
[J,Ai,Pi,Oi] = h.evalError(h.rho_i,res);
disp('init optim')
h.runOptimization(res,false);
rho_last = h.rho_o;
[J,Ai,Pi,Oi] = h.evalError(h.rho_o,res);
disp(['initial residual: ' num2str(J)])
%%
Rvec = 0.5:-Rinc:Rinc;

nv = length(Rvec);

Jvec = nan(1,1,nv);
Avec = nan([size(Ai) nv]);
Pvec = nan([size(Pi) nv]);
Ovec = nan([size(Oi) nv]);



for n = 1:nv
    disp(['init optim ' num2str(n) ' of ' num2str(nv)])
    R = Rvec(n);
    
    h = genTangOptim(3,R,'h3_3_in');
    optimSet3_3_in(h);

    h.getTransferFunctions;
    h.setRhoInit(rho_last);
    
    h.runOptimization(res,false);
%     disp('ref optim')
%     h.continueOptimization(10^4,false);
    [J,Ai,Pi,Oi] = h.evalError(h.rho_o,res);
    disp(['optim ended in ' num2str(J)])
    Jvec(:,:,n) = J;
    Avec(:,:,n) = Ai;
    Pvec(:,:,n) = Pi;
    Ovec(:,:,n) = Oi;
    rho_last = h.rho_o;
    % Generate transfer function and derivative function 
    
end


R3_3_in_cont.Rvec = Rvec;
R3_3_in_cont.Jvec = Jvec;
R3_3_in_cont.Avec = Avec;
R3_3_in_cont.Pvec = Pvec;
R3_3_in_cont.Ovec = Ovec;

save('R3_3_in_cont','R3_3_in_cont')

%%
for r = 1:2:3
figure
for n = 1:4
subplot(5,2,2*n-1)
plot(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Avec(n,r,:))));
hold on
plot(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Avec(n,r,:))));

subplot(5,2,2*n)
plot(R3_3_in_cont.Rvec,(squeeze(R3_3_in_cont.Pvec(n,r,:))));
hold on
plot(R3_3_dist_cont.Rvec,(squeeze(R3_3_dist_cont.Pvec(n,r,:))));
end
subplot(5,2,9)
plot(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Ovec(:,:,:))));
hold on
plot(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Ovec(:,:,:))));
subplot(5,2,10)
semilogy(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Jvec(:,:,:))));
hold on
semilogy(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Jvec(:,:,:))));
ylim([10^-6 10^-1])

end

%%

figure

subplot(3,4,1)
plot(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Avec(3,1,:))));
hold on
plot(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Avec(3,1,:))));

subplot(3,4,2)
plot(R3_3_in_cont.Rvec,(squeeze(R3_3_in_cont.Pvec(3,1,:))));
hold on
plot(R3_3_dist_cont.Rvec,(squeeze(R3_3_dist_cont.Pvec(3,1,:))));

subplot(3,4,3)
plot(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Avec(3,3,:))));
hold on
plot(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Avec(3,3,:))));

subplot(3,4,4)
plot(R3_3_in_cont.Rvec,(squeeze(R3_3_in_cont.Pvec(3,3,:))));
hold on
plot(R3_3_dist_cont.Rvec,(squeeze(R3_3_dist_cont.Pvec(3,3,:))));
% theta
subplot(3,4,5)
plot(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Avec(4,1,:))));
hold on
plot(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Avec(4,1,:))));

subplot(3,4,6)
plot(R3_3_in_cont.Rvec,(squeeze(R3_3_in_cont.Pvec(4,1,:))));
hold on
plot(R3_3_dist_cont.Rvec,(squeeze(R3_3_dist_cont.Pvec(4,1,:))));

subplot(3,4,7)
plot(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Avec(4,3,:))));
hold on
plot(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Avec(4,3,:))));

subplot(3,4,8)
plot(R3_3_in_cont.Rvec,(squeeze(R3_3_in_cont.Pvec(4,3,:))));
hold on
plot(R3_3_dist_cont.Rvec,(squeeze(R3_3_dist_cont.Pvec(4,3,:))));


subplot(3,4,10)
plot(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Ovec(:,:,:))));
hold on
plot(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Ovec(:,:,:))));
subplot(3,4,11)
semilogy(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Jvec(:,:,:))));
hold on
semilogy(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Jvec(:,:,:))));
ylim([10^-6 10^-1])

matlab2tikz(['optim_3_in_dist.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);
    
 %%

figure

subplot(3,2,1)
plot(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Avec(3,1,:))));
hold on
plot(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Avec(3,1,:))));


subplot(3,2,2)
plot(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Avec(3,3,:))));
hold on
plot(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Avec(3,3,:))));


% theta
subplot(3,2,3)
plot(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Avec(4,1,:))));
hold on
plot(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Avec(4,1,:))));


subplot(3,2,4)
plot(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Avec(4,3,:))));
hold on
plot(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Avec(4,3,:))));



subplot(3,2,5)
plot(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Ovec(:,:,:))));
hold on
plot(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Ovec(:,:,:))));
subplot(3,2,6)
semilogy(R3_3_in_cont.Rvec,abs(squeeze(R3_3_in_cont.Jvec(:,:,:))));
hold on
semilogy(R3_3_dist_cont.Rvec,abs(squeeze(R3_3_dist_cont.Jvec(:,:,:))));
ylim([10^-6 10^-1])



%% Plot Fourier coefficients
dataFolder = 'D:\afstudeerdata\simulations\tang_xyphi_3_R\';
load([dataFolder 'gatherFourier'])
nSims = size(F,3);
nReps = size(F,2);
P = wrapTo2Pi(angle(F)+0.5*pi);
A = 2*abs(F);
P3 = wrapTo2Pi(angle(F3)+0.5*pi);
A3 = 2*abs(F3);


for nr = 1:12
subplot(3,2,1)
plot(R,squeeze(A(3,nr,:)),'+','Color',cBlack)
subplot(3,2,3)
    plot(R,squeeze(A(4,nr,:)),'+','Color',cBlack)
subplot(3,2,5)
    plot(R,squeeze(O(1,nr,:)),'+','Color',cBlack)

subplot(3,2,2)
plot(R,squeeze(A3(3,nr,:)),'+','Color',cBlack)
subplot(3,2,4)
    plot(R,squeeze(A3(4,nr,:)),'+','Color',cBlack)
end
%%

matlab2tikz(['optim_3_in_dist.tex'],'parseStrings',false,...
        'height','\figureheight',...
        'width','\figurewidth',...
        'showInfo', false);

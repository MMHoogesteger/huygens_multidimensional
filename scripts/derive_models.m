%clear all
clc
tic;
NMIN = 1;
NTOT = 99;
derive_matlab = 0;
%%
mdlset.m3d = 1;
mdlset.mphi = 1;
mdlset.mxy = 1;
mdlset.mxphi = 1;
mdlset.myphi = 1;
mdlset.mx = 1;
mdlset.my = 1;

%%

tloop = zeros(1,NTOT-NMIN+1);
for N = NMIN:NTOT
    mdl.N = N;
    mdl.PSI_zero = 0;
    mdl.GAMMA_zero = 0;
    if(mdlset.m3d)
        tloopstart = tic;
        disp(['Deriving models for ' num2str(N) ' metronomes.'])
        disp(['3D-model'])
        mdl.DOF_x = 1;
        mdl.DOF_y = 1;
        mdl.DOF_phi = 1;
        mdl.modelName = ['mxyphi_' num2str(N)];
        derive_model(mdl,derive_matlab);
    end
    if(mdlset.mphi)
        disp(['Rotational-model'])
        mdl.DOF_x = 0;
        mdl.DOF_y = 0;
        mdl.DOF_phi = 1;
        mdl.modelName = ['mphi_' num2str(N)];
        derive_model(mdl,derive_matlab);
    end
    if(mdlset.mxy)
        disp(['Translational-model'])
        mdl.DOF_x = 1;
        mdl.DOF_y = 1;
        mdl.DOF_phi = 0;
        mdl.modelName = ['mxy_' num2str(N)];
        derive_model(mdl,derive_matlab);
    end
    if(mdlset.mxphi)
        disp(['1D translation x plus rotation-model'])
        mdl.DOF_x = 1;
        mdl.DOF_y = 0;
        mdl.DOF_phi = 1;
        mdl.modelName = ['mxphi_' num2str(N)];
        derive_model(mdl,derive_matlab);
    end
    if(mdlset.myphi)
        disp(['1D translation y plus rotation-model'])
        mdl.DOF_x = 0;
        mdl.DOF_y = 1;
        mdl.DOF_phi = 1;
        mdl.modelName = ['myphi_' num2str(N)];
        derive_model(mdl,derive_matlab);
    end
    if(mdlset.mx)
        disp(['1D-model x'])
        mdl.DOF_x = 1;
        mdl.DOF_y = 0;
        mdl.DOF_phi = 0;
        mdl.modelName = ['mx_' num2str(N)];
        derive_model(mdl,derive_matlab);
    end
    if(mdlset.my)
        disp(['1D-model y'])
        mdl.DOF_x = 0;
        mdl.DOF_y = 1;
        mdl.DOF_phi = 0;
        mdl.modelName = ['my_' num2str(N)];
        derive_model(mdl,derive_matlab);
    end
    tloop(N-NMIN+1) = toc(tloopstart);
    ttot = toc;
    disp(['Deriving all models for ' num2str(N) ' metronomes took ' num2str(tloop(N-NMIN+1)) ' seconds.'])
    disp(['Total running time until now is ' num2str(ttot) ' seconds.'])
    
    if(N-NMIN+1<3)
        tavg = ttot/(N-NMIN+1);
        tinc = (tloop(N-NMIN+1)-tavg)/((N-NMIN+1)/2);
        test = tavg*(NTOT-N)+0.5*tinc*(NTOT-N)^2;
    else
        P = polyfit(NMIN:N,tloop(1:N-NMIN+1),1);
        test = sum(polyval(P,N+1:NTOT));
        plot(NMIN:N,tloop(1:N-NMIN+1),'x')
        hold on
        plot(NMIN:NTOT,polyval(P,NMIN:NTOT),':')
        hold off
        drawnow
    end
    disp(['Estimated time left:' num2str(test) ' seconds, which is '  num2str(test/3600) ' hours.'])
    
    
    
end

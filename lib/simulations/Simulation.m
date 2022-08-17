classdef Simulation < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name string;
        
        % Paths
        modelPath string;
        
        % Model
        N uint8;
        modelName string;
        model struct;
        enteredModel struct;
        
        % Configuration
        R;
        gamma;
        psi;
        configuration;
        parameters;
        escapement;
        T;
        fs;
        tspan;
        U_gains;
        HStar;
        
        % Booleans
        hasModel logical;
        hasParameters logical;
        hasConfiguration logical;
        
        % Repetitions
        nReps uint8;
        reps struct;
        
        % Simulation options
        simOptions;
        
        % Results
        results struct;
        
        
    end
    
    methods
        function obj = Simulation(name)
            %Simulation Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 0
                obj.name = 'default';
            else
                obj.name = name;
            end
            obj.hasParameters = false;
            obj.hasModel = false;
            obj.hasConfiguration = false;
            
            obj.reps = struct;
            
            obj.simOptions.relTol = 1e-10;
            obj.simOptions.absTol = 1e-10;
        end
        
        function name = getLatexName(obj,nq)
            name = obj.model.q_array{nq};
            if(startsWith(name,"phi"))
                name = ['\' name];
            end
            if(startsWith(name,"theta"))
                name = ['\theta_{' extractAfter(name,"theta") '}'];
            end
        end
        
        function setModel(obj,N,x,y,phi)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.N = N;
            tempName = 'm';
            if(x)
                tempName = [tempName 'x'];
            end
            if(y)
                tempName = [tempName 'y'];
            end
            if(phi)
                tempName = [tempName 'phi'];
            end
            obj.modelName = [tempName '_' num2str(N)];
            tload = load(obj.modelName);
            obj.model = tload.mdl;
            setParameters(obj,getMassParams(getparams(),N));
        end
        
        function setConfiguration(obj,R,gamma,psi)
            if(numel(R)~=obj.N)
                error('Simulation:configurationError','Incorrect length of R');
            end
            if(numel(gamma)~=obj.N)
                error('Simulation:configurationError','Incorrect length of Gamma');
            end
            if(numel(psi)~=obj.N)
                error('Simulation:configurationError','Incorrect length of Psi');
            end
            
            c = obj.model.configuration;
            
            obj.R = R;
            obj.gamma = gamma;
            obj.psi = psi;
            
            % Make an array of the correct order of configuration
            % parameters
            conf = zeros(obj.N*length(c),1);
            for i=1:length(c)
                conf((i-1)*obj.N+1:i*obj.N) = obj.(c{i});
            end
            
            obj.configuration = conf;
        end
        
        function setParameters(obj,params)
            p = obj.model.parameters;
            
            % Make an array of the correct order of parameters
            par = zeros(length(p),1);
            for i=1:length(p)
                par(i) = params.(p{i});
            end
            
            obj.parameters = par;
        end
        
        function setSimParams(obj,fs,T)
            obj.T = T;
            obj.fs = fs;
            obj.tspan = 0:1/obj.fs:obj.T;
        end
        
        function setEscapement(obj,escapement,HStar,U_gains)
            obj.escapement = escapement;
            if(isequal(escapement,"Hamiltonian"))
                if(size(U_gains,1)~=obj.model.N)
                    error('Simulation:escapementError','Incorrect length of U_gains');
                end
                if(size(HStar,1)~=obj.model.N)
                    error('Simulation:escapementError','Incorrect length of HStar');
                end
                obj.U_gains = U_gains;
                obj.HStar = HStar;
            end
        end
        
        function setInitialConditions(obj,nReps,xyphi,theta,dxyphi,dtheta)
            obj.nReps = nReps;
            if(size(xyphi,2)~=nReps)
                error('Simulation:initialConditionError','Incorrect number of xyphi');
            end
            if(size(dxyphi,2)~=nReps)
                error('Simulation:initialConditionError','Incorrect number of dxyphi');
            end
            if(size(theta,2)~=nReps)
                error('Simulation:initialConditionError','Incorrect number of theta');
            end
            if(size(dtheta,2)~=nReps)
                error('Simulation:initialConditionError','Incorrect number of dtheta');
            end
            
            if(size(xyphi,1)~=obj.model.NqP)
                error('Simulation:initialConditionError','Incorrect length of xyphi');
            end
            if(size(dxyphi,1)~=obj.model.NqP)
                error('Simulation:initialConditionError','Incorrect length of dxyphi');
            end
            if(size(theta,1)~=obj.model.N)
                error('Simulation:initialConditionError','Incorrect length of theta');
            end
            if(size(dtheta,1)~=obj.model.N)
                error('Simulation:initialConditionError','Incorrect length of dtheta');
            end
            
            for k = 1:nReps
                obj.reps(k).xyphi = xyphi(:,k);
                obj.reps(k).dxyphi = dxyphi(:,k);
                obj.reps(k).theta = theta(:,k);
                obj.reps(k).dtheta = dtheta(:,k);
                obj.reps(k).z0 = [xyphi(:,k);theta(:,k);dxyphi(:,k);dtheta(:,k)];
            end
        end
        
        function setInitialConditionsZeroVelocity(obj,nReps,xyphi,theta)
            obj.setInitialConditions(nReps,xyphi,theta,0*xyphi,0*theta);
        end
        
        function plotCircleFit(obj,repId,Tss,saveplot,saveName)
            
            sol = obj.results(repId).solution;
            circle = obj.results(repId).analysis.circleFit;
            
            sst = sol.t(sol.t>Tss);
            sst = sst-sst(1);
            
            sstN = sol.theta(sol.t>Tss);
            ssdtN = sol.dtheta(sol.t>Tss);
            
            
            alpha = 0:0.01:2*pi;
            x = cos(alpha).*circle.thetaMax;
            y = sin(alpha).*circle.dthetaMax;
            
            figure
            plot(sstN(1:10:end), ssdtN(1:10:end),'k')
            hold on
            
            
            plot(x,y,'r--')
            grid on
            %title('Steady-state solution and circle fit')
            xlabel('$\theta [rad]$','Interpreter','Latex')
            ylabel('$\dot{\theta} [rad/s]$','Interpreter','Latex')
            legend(['Solution for t > ' num2str(Tss)],'Circle fit')
            if(saveplot)
                matlab2tikz([saveName '_circle_fit.tex'],'parseStrings',false,...
                    'height','\figureheight',...
                    'width','\figurewidth',...
                    'showInfo', false);
            end
            
        end
        
        function plotEnergy(obj,repId,Tstart,Tend,fsplot,saveplot,saveName)
            
            fsidx = ceil(obj.fs/fsplot);
            
            sol = obj.results(repId).solution;
            env = obj.results(repId).analysis.energyEnvelopes;
            
            
            idx = Tstart*obj.fs+1:fsidx:Tend*obj.fs;
            
            tplot = sol.t(idx);
            Hplot = sol.H(idx);
            Hmplot = sol.Hm(idx,:);
            
            idxE = (env.t > Tstart & env.t < Tend);
            
            
            figure;
            plot(tplot,Hplot,'k')
            hold on
            plot(env.t(idxE),env.maxH(idxE),'r')
            plot(env.t(idxE),env.minH(idxE),'r--')
            %title('Total energy')
            xlabel('$t [s]$','Interpreter','Latex')
            ylabel('$H [J]$','Interpreter','Latex')
            xlim([Tstart Tend])
            grid on
            legend('Energy','Max','Min')
            drawnow
            
            if(saveplot)
                matlab2tikz([saveName '_total_energy.tex'],'parseStrings',false,...
                    'height','\figureheight',...
                    'width','\figurewidth',...
                    'showInfo', false);
            end
            
            
            figure;
            Nq = obj.model.Nq;
            for n = 1:Nq
                subplot(Nq,1,n)
                plot(tplot,Hmplot(:,n),'k')
                hold on
                plot(env.t(idxE),env.maxHm((idxE),n),'r')
                plot(env.t(idxE),env.minHm((idxE),n),'r--')
                ylabel(['$H_{' obj.getLatexName(n) '} [J]$'],'Interpreter','Latex')
                xlim([Tstart Tend])
                grid on
                legend('hide')
            end
            xlabel('$t [s]$','Interpreter','Latex')
            subplot(Nq,1,1)
            %title('Energy of the different DOFs')
            legend('Energy','Max','Min')
            drawnow
            
            if(saveplot)
                matlab2tikz([saveName '_dof_energy.tex'],'parseStrings',false,...
                    'height','\figureheight',...
                    'width','\figurewidth',...
                    'showInfo', false);
            end
        end
        
        
        
        function plotEnergyDeviations(obj,repId,Tstart,Tend,saveplot,saveName)
            
            
            env = obj.results(repId).analysis.energyEnvelopes;
            
            idxE = (env.t > Tstart & env.t < Tend);
            
            
            figure;
            semilogy(env.t(idxE),env.dmaxH(idxE),'k')
            hold on
            semilogy(env.t(idxE),env.dminH(idxE),'k--')
            %title('Deviation from steady-state total energy')
            xlabel('$t [s]$','Interpreter','Latex')
            ylabel('$\Delta H [J]$','Interpreter','Latex')
            xlim([Tstart Tend])
            grid on
            legend('Max','Min')
            drawnow
            
            if(saveplot)
                matlab2tikz([saveName '_total_energy_deviation.tex'],'parseStrings',false,...
                    'height','\figureheight',...
                    'width','\figurewidth',...
                    'showInfo', false);
            end
            
            figure;
            Nq = obj.model.Nq;
            for n = 1:Nq
                subplot(Nq,1,n)
                semilogy(env.t(idxE),env.dmaxHm((idxE),n),'k')
                hold on
                
                semilogy(env.t(idxE),env.dminHm((idxE),n),'k--')
                ylabel(['$\Delta H_{' obj.getLatexName(n) '} [J]$'],'Interpreter','Latex')
                xlim([Tstart Tend])
                grid on
                legend('hide')
            end
            xlabel('$t [s]$','Interpreter','Latex')
            subplot(Nq,1,1)
            %title('Deviation from steady-state energy of the different DOFs ')
            legend('Max','Min')
            drawnow
            
            if(saveplot)
                matlab2tikz([saveName '_dof_energy_deviation.tex'],'parseStrings',false,...
                    'height','\figureheight',...
                    'width','\figurewidth',...
                    'showInfo', false);
            end
        end
        function plotFourierAnalysis(obj,repId,saveplot,saveName)
            
            fourier = obj.results(repId).analysis.fourier;
            figure;
            Nq = obj.model.Nq;
            for n = 1:Nq
                
                subplot(Nq,1,n)
                semilogy(abs(fourier.Fzk(n,:)),'kx')
                ylabel(['$|a_{k' obj.getLatexName(n) '}| [rad]$'],'Interpreter','Latex')
                legend('hide')
                
                grid on
                
                
            end
            
            xlabel('$k [-]$','Interpreter','Latex')
            subplot(Nq,1,1)
            %title('Fourier Series coefficients','Interpreter','Latex')
            drawnow
            
            if(saveplot)
                matlab2tikz([saveName '_fourier_coeff.tex'],'parseStrings',false,...
                    'height','\figureheight',...
                    'width','\figurewidth',...
                    'showInfo', false);
            end
            
            figure;
            for n = 1:Nq
                
                subplot(Nq,1,n)
                semilogy(abs(fourier.Ezk(n,:)),'ro')
                ylabel(['Error $' obj.getLatexName(n) ' [rad]$'],'Interpreter','Latex')
                legend('hide')
                
                grid on
                
                
            end
            
            xlabel('$r [-]$','Interpreter','Latex')
            subplot(Nq,1,1)
            %title('Errors of Fourier Series approximations','Interpreter','Latex')
            drawnow
            
            if(saveplot)
                matlab2tikz([saveName '_fourier_coeff_err.tex'],'parseStrings',false,...
                    'height','\figureheight',...
                    'width','\figurewidth',...
                    'showInfo', false);
            end
        end
        
        function simulate(obj,verbose)
            reps = obj.reps;
            emdl = obj.enteredModel;
            opt = odeset('Mass',@(t,z) emdl.fMM(t,z),...
                'MStateDependence','strong','RelTol',obj.simOptions.relTol,'AbsTol',obj.simOptions.absTol);
            tVec = obj.tspan;
            U_gains = obj.U_gains;
            HStar = obj.HStar; %#ok<*PROP>
            escapement = obj.escapement;
            simSolutions = cell(1,obj.nReps);
            parfor n = 1:obj.nReps
                tstart = tic;
                rep = reps(n);
                if(verbose)
                fprintf(['Starting simulation of repetition ' num2str(n) '.\n'])
                end
                sol = Simulation.simulateRepetition(rep,emdl,opt,tVec,U_gains,HStar,escapement);
                if(verbose)
                fprintf(['Finished simulation of repetition ' num2str(n) ' in ' num2str(toc(tstart)) ' seconds.\n'])
                end
                simSolutions(1,n) = {sol};
            end
            obj.results = cell2struct(simSolutions,{'solution'},1);
        end
        
        function analyze(obj,Tss,verbose)
            res = obj.results;
            emdl = obj.enteredModel;
            fs = obj.fs;
            simSolutions = cell(2,obj.nReps);
            parfor n = 1:obj.nReps
                tstart = tic;
                sol = res(n).solution;
                if(verbose)
                fprintf(['Starting analysis of repetition ' num2str(n) '.\n'])
                end
                ana = Simulation.analyzeRepetition(sol,emdl,Tss,fs);
                if(verbose)
                fprintf(['Finished analysis of repetition ' num2str(n) ' in ' num2str(toc(tstart)) ' seconds.\n'])
                end
                simSolutions(:,n) = {sol;ana};
            end
            obj.results = cell2struct(simSolutions,{'solution','analysis'},1);
        end
        
        
        
        function generateModel(obj)
            mdl = obj.model;
            par = obj.parameters;
            conf = obj.configuration;
            
            emdl.modelName = ['emdl_' obj.name];
            
            emdl.mdl = mdl;
            emdl.par = par;
            emdl.conf = conf;
            
            s = mdl.symbolics;
            
            U = s.U;
            z = s.z;
            t = s.t;
            U_gains = s.U_gains;
            HStar = s.HStar;
            
            % Equations
            MM = mdl.fMM(t,z,par,conf);
            RHS = mdl.fRHS(t,z,par,conf,U);
            U_OR = mdl.fU_OR(t,z,par,conf);
            U_HAM = mdl.fU_HAM(t,z,par,conf,U_gains,HStar);
            H = mdl.fH(t,z,par,conf);
            HM = mdl.fHM(t,z,par,conf);
            
            % Write equations to emdl struct
            emdl.fMM = matlabFunction(MM,'vars',{t,z},'outputs',{'M'});
            emdl.fRHS = matlabFunction(RHS,'vars',{t,z,U},'outputs',{'dxdt'});
            emdl.fU_OR = matlabFunction(U_OR,'vars',{t,z},'outputs',{'u'});
            emdl.fU_HAM = matlabFunction(U_HAM,'vars',{t,z,U_gains,HStar},'outputs',{'u'});
            emdl.fH = matlabFunction(H,'vars',{t,z},'outputs',{'H'});
            emdl.fHM = matlabFunction(HM,'vars',{t,z},'outputs',{'HM'});
            
            obj.enteredModel = emdl;
            
        end
        
        function save(obj,folder,fPath)
            
            if(nargin<3)
                fPath = obj.name;
            end
            if(~(exist(folder, 'dir')==7))
                mkdir(folder);
            end
            if(exist(strcat(folder,fPath), 'file')==6)
                error('File exists');
            end
            s = obj; %#ok<NASGU>
            save(strcat(folder,fPath),'s','-v7.3');
        end
        
        
        
    end
    
    
    methods(Static)
        function simAnalysis = analyzeRepetition(sol,emdl,Tss,fs)
            simAnalysis = struct;
            simAnalysis.energyEnvelopes = Simulation.getEnergyEnvelopes(sol,emdl.mdl.Nq);
            simAnalysis.circleFit = Simulation.getCircleFit(sol,Tss,fs);
            simAnalysis.fourier = Simulation.getFourierAnalysis(sol,Tss,simAnalysis.circleFit,emdl.mdl.NqP+1);
        end
        
        function circle = getCircleFit(sol,Tss,fs)
            %% Fit circle
            
            sst = sol.t(sol.t>Tss);
            sst = sst-sst(1);
            
            sstN = sol.theta(sol.t>Tss);
            ssdtN = sol.dtheta(sol.t>Tss);
            
            circle.thetaMax = mean(abs(sstN(abs(ssdtN)<0.01*max(ssdtN))));
            circle.dthetaMax = mean(abs(ssdtN(abs(sstN)<0.01*max(sstN))));
            
            % Get frequency
            % total angle
            pos = [sstN./circle.thetaMax,ssdtN./circle.dthetaMax];
            pos1 = pos(2:end,:);
            pos2 = pos(1:end-1,:);
            dpos = pos1-pos2;
            dpos = sqrt(dpos(:,1).^2+dpos(:,2).^2);
            totangle = sum(dpos);
            circle.freq = totangle/2/pi/length(sst)*fs;
        end
        
        function env = getEnergyEnvelopes(sol,Nq)
            T = max(sol.t);
            tinc = 2; % seconds
            
            % Number of intervals
            L = floor(T/tinc);
            
            
            % Calculate in each interval the maximum, mean and minimum value of H and
            % H - Hm.
            
            t = zeros(L,1);
            
            maxH = zeros(L,1);
            minH = zeros(L,1);
            
            maxHm = zeros(L,Nq);
            minHm = zeros(L,Nq);
            
            for l = 1:L
                idx = (sol.t>((l-1)*tinc+1e-8) & sol.t <=(l*tinc));
                t(l) = mean(sol.t(idx));
                maxH(l) = max(sol.H(idx));
                minH(l) = min(sol.H(idx));
                
                [maxHm(l,:)] = max(sol.Hm(idx,:));
                [minHm(l,:)] = min(sol.Hm(idx,:));
                
                
            end
            
            % Calculate deviations from end value
            dmaxH = abs(maxH-maxH(end));
            dminH = abs(minH-minH(end));
            
            dmaxHm = abs(maxHm-maxHm(end,:));
            dminHm = abs(minHm-minHm(end,:));
            
            env.t = t;
            
            env.maxH = maxH;
            env.minH = minH;
            
            env.maxHm = maxHm;
            env.minHm = minHm;
            
            env.dmaxH = dmaxH;
            env.dminH = dminH;
            
            env.dmaxHm = dmaxHm;
            env.dminHm = dminHm;
            
        end
        
        function fourier = getFourierAnalysis(sol,Tss,circle,normId)
            % Initial guess for omega from circle
            f0 = circle.freq;
            omega0 = 2*pi*f0;
            
            R = 15;
            % idx = steady state data
            idx = sol.t>Tss;
            
            % First use the first metronome to find frequency
            t = sol.t(idx);
            T = t(end)-t(1);
            theta1 = sol.theta(idx,1);
            
            fc0 = @(omega) fourierDC(t,theta1,omega);
            fc = @(omega) fourierIntegral(t,theta1,omega,R);
            fJ = @(omega) sum((theta1.'-fourierReconstruct(t,fc0(omega),fc(omega),omega)).^2);
            options = optimset('TolFun',1e-10,'TolX',1e-10);
            omega = fminsearch(fJ,omega0,options);
            
            
            Nq = size(sol.z,2);
            
            Fz0 = zeros(Nq,1);
            Fzk = zeros(Nq,R);
            Ezk = zeros(Nq,R);
            
            for n = 1:Nq
                f = sol.z(idx,n);
                c0 = fourierDC(t,f,omega);
                c = fourierIntegral(t,f,omega,15);
                
                Fz0(n,:) = c0;
                Fzk(n,:) = c;
                
                for r = 1:R
                    fr = fourierReconstruct(t,c0,c(1:r),omega);
                    Ezk(n,r) = sqrt(trapz(t,(f-fr.').^2)./T);
                end
            end
            Fzkn = normalizeFourierCoefficients(Fzk,Fzk(normId,1));
            
            fourier.R = R;
            fourier.T = T;
            fourier.omega = omega;
            fourier.Fz0 = Fz0;
            fourier.Fzk = Fzk;
            fourier.Fzkn = Fzkn;
            fourier.Ezk = Ezk;
            
        end
        
        
        
        function simRes = simulateRepetition(rep,emdl,opt,tVec,U_gains,HStar,escapement)
            
            fU = []; %#ok<NASGU>
            switch escapement
                case "Original"
                    fU = @(t,z) emdl.fU_OR(t,z);
                case "Hamiltonian"
                    if(isempty(U_gains))
                        error('no escapement gain set');
                    elseif(isempty(HStar))
                        error('no escapement target set');
                    end
                    fU = @(t,z) emdl.fU_HAM(t,z,U_gains,HStar);
                otherwise
                    error('not a valid escapement');
            end
            fRHS = @(t,z) emdl.fRHS(t,z,fU(t,z));
            
            
            
            simRes = struct;
            [simRes.t, simRes.z] = ode15s(fRHS,tVec,rep.z0,opt);
            
            NqP = emdl.mdl.NqP;
            Nq = emdl.mdl.Nq;
            
            simRes.xyphi = simRes.z(:,1:NqP);
            simRes.theta = simRes.z(:,NqP+1:Nq);
            simRes.dxyphi = simRes.z(:,Nq+1:Nq+NqP);
            simRes.dtheta = simRes.z(:,Nq+NqP+1:2*Nq);
            
            k = 0;
            if(emdl.mdl.DOF_x)
                k = k+1;
                simRes.x = simRes.xyphi(:,k);
                simRes.dx = simRes.dxyphi(:,k);
            end
            if(emdl.mdl.DOF_y)
                k = k+1;
                simRes.y = simRes.xyphi(:,k);
                simRes.dy = simRes.dxyphi(:,k);
            end
            if(emdl.mdl.DOF_phi)
                k = k+1;
                simRes.phi = simRes.xyphi(:,k);
                simRes.dphi = simRes.dxyphi(:,k);
            end
            if(k~=size(simRes.xyphi,2))
                error('Error during extraction of results');
            end
            
            
            t = simRes.t;
            z = simRes.z;
            
            
            U = zeros(numel(t),emdl.mdl.Nq);
            H = zeros(numel(t),1);
            Hm = zeros(numel(t),emdl.mdl.Nq);
            
            for l = 1:numel(t)
                U(l,:) = fU(t(l),z(l,:).').';
                H(l) = emdl.fH(t(l),z(l,:).');
                Hm(l,:)= emdl.fHM(t(l),z(l,:).').';
                
            end
            
            simRes.U = U;
            simRes.H = H;
            simRes.Hm = Hm;
            
        end
        
    end
end


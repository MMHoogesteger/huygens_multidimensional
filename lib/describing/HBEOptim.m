classdef HBEOptim < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name string;
        
        % Model
        N uint8;
        r %order
        
        % Configuration
        R;
        gamma;
        psi;
        U_gains;
        HStar;
        params;
        
        % Booleans
        hasConfiguration logical;
        hasEscapement logical;
        hasModel logical;
        
        % Function handles
        G_func;
        dGdomega_func;
        
        % vars
        A_contr
        P_contr
        O_contr
        
        A_fixed
        P_fixed
        O_fixed
        
        % contr =0, fixed =0 --> optimize
        % contr =0, fixed =1 <<< Impossible
        % contr =1, fixed =0 --> refer other metronome
        % contr =1, fixed =1 --> Fixed value
        
        % contr is wether value is dependent or independent (dof/unknown)
        % fixed is wether value is fixed or dependent
        
        A_ref
        P_ref
        
        A
        P
        Pshift
        O
        
        As
        Ps
        Os
        
        A_scale
        P_scale
        O_scale
        
        no_A
        no_P
        no_O
        no_T
        
        rho_i
        rho_o
        
        % Results
        results struct;
        
        
    end
    
    methods
        function obj = HBEOptim(name)
            %Simulation Construct an instance of this class
            %   Detailed explanation goes here
            if nargin == 0
                obj.name = 'default';
            else
                obj.name = name;
            end
            obj.hasConfiguration = false;
            obj.hasEscapement = false;
            obj.hasModel = false;
        end
                      
        function setModel(obj,N)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.N = N;
            obj.params = getMassParams(getparams(),N);
            obj.hasModel = true;
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
                        
            obj.R = R;
            obj.gamma = gamma;
            obj.psi = psi;
            obj.hasConfiguration = true;
        end
                        
        function setEscapement(obj,HStar,U_gains)
            if(size(U_gains,1)~=obj.N)
                error('Simulation:escapementError','Incorrect length of U_gains');
            end
            if(size(HStar,1)~=obj.N)
                error('Simulation:escapementError','Incorrect length of HStar');
            end
            obj.U_gains = U_gains;
            obj.HStar = HStar;
            obj.hasEscapement = true;
        end
        
        function getTransferFunctions(obj)
            % Generate transfer function and derivative function 
            if(~obj.hasConfiguration)
                error('Simulation:configurationError','No configuration set');
            end
            if(~obj.hasModel)
                error('Simulation:configurationError','No model set');
            end
            syms symomega
            [~,Gs] = getTransfer_n_xyphi(obj.N,symomega,1,obj.params,obj.R,obj.gamma,obj.psi);
            Gds = diff(Gs,symomega);
            obj.G_func = matlabFunction(Gs,'file',['temp/G_func'],'vars',symomega,'outputs',{'G'});
            obj.dGdomega_func = matlabFunction(Gds,'file',['temp/dGdomega_func'],'vars',symomega,'outputs',{'dGdomega'});
            
        end
        
        function setOrder(obj,r)
            obj.r = r;
            obj.A_contr = false(3+obj.N,r);
            obj.P_contr = false(3+obj.N,r);
            obj.O_contr = false;
            
            obj.A_fixed = false(3+obj.N,r);
            obj.P_fixed = false(3+obj.N,r);
            obj.O_fixed = false;
            
            obj.A_ref = nan(3+obj.N,r);
            obj.P_ref = nan(3+obj.N,r);
                        
            obj.A = nan(3+obj.N,r);
            obj.P = nan(3+obj.N,r);
            obj.Pshift = nan(3+obj.N,r);
            obj.O = NaN;
            
            obj.A_scale = ones(3+obj.N,r);
            obj.P_scale = ones(3+obj.N,r);
            obj.O_scale = 10;
            
            obj.A_scale(1:2,:) = 1/1000;
            obj.A_scale(3,:) = 1/100;
            if(r>2)
                obj.A_scale(:,3) = obj.A_scale(:,3);
            end
            if(r>5)
                obj.A_scale(:,6) = obj.A_scale(:,6);
            end
            obj.setMetPhaseToZero(1,1);
            
        end
        
        function assembleInitial(obj)
            optim_A_ids = ~obj.A_contr;
            optim_P_ids = ~obj.P_contr;
            optim_O_ids = ~obj.O_contr;
            
            obj.no_A = sum(sum(optim_A_ids));            
            obj.no_P = sum(sum(optim_P_ids));
            obj.no_O = sum(sum(optim_O_ids));
            
            obj.no_T = obj.no_A + obj.no_P + obj.no_O;
            
            obj.rho_i = nan(obj.no_T,1);
            
            obj.As = obj.A./obj.A_scale;
            obj.Ps = obj.P./obj.P_scale;
            obj.Os = obj.O./obj.O_scale;
            
            obj.rho_i(1:obj.no_A) = obj.As(optim_A_ids);
            obj.rho_i(obj.no_A+1:obj.no_A+obj.no_P) = obj.Ps(optim_P_ids);
            obj.rho_i(obj.no_A+obj.no_P+1:obj.no_T) = obj.Os(optim_O_ids);
            
        end
        
        
        
        
        function [J,Ai,Pi,Oi] = evalError(obj,rho,res)
            Ai = nan(3+obj.N,obj.r);
            Pi = nan(3+obj.N,obj.r);
            Oi = NaN;
            
            optim_A_ids = ~obj.A_contr;
            optim_P_ids = ~obj.P_contr;
            optim_O_ids = ~obj.O_contr;
            
            ref_A_ids = obj.A_contr & (~obj.A_fixed);
            ref_P_ids = obj.P_contr & (~obj.P_fixed);
            
            % Enter all fixed vars            
            Ai(obj.A_fixed) = obj.As(obj.A_fixed);
            Pi(obj.P_fixed) = obj.Ps(obj.P_fixed);
            Oi(obj.O_fixed) = obj.Os(obj.O_fixed);
            
            % Enter all optim vars
            Ai(optim_A_ids) = rho(1:obj.no_A);
            Pi(optim_P_ids) = rho(obj.no_A+1:obj.no_A+obj.no_P);
            Oi(optim_O_ids) = rho(obj.no_A+obj.no_P+1:obj.no_T);
            
            % Scale
            Ai = Ai.*obj.A_scale;
            Pi = Pi.*obj.P_scale;
            Oi = Oi.*obj.O_scale;
            Nq = 3+obj.N;
            W = nan(Nq*obj.r,1);
            % Weight
            for k = 1:obj.r
                nids = (k-1)*Nq+(1:Nq).';
                W(nids,1) = 1./obj.A_scale(:,k);
            end
            
            % Replace all met ref vars
            for n = 4:3+obj.N
                for k = 1:obj.r
                    if (ref_A_ids(n,k))
                        Ai(n,k) = Ai(obj.A_ref(n,k),k);
                    end
                    if (ref_P_ids(n,k))
                        Pi(n,k) = Pi(obj.P_ref(n,k),k)+obj.Pshift(n,k)*k;
                    end
                end
            end 
            
            % Transform to complex vars:
            ak = convertAmpPhaseToComplex(Ai,Pi);
            omega = Oi;


            res = 2^nextpow2(res)-1;

            [Jk] = harmonic_n_xyphi_eval_r_fft(obj.N,ak,omega,...
                obj.G_func,obj.dGdomega_func,obj.params,obj.U_gains,...
                obj.HStar,obj.R,obj.gamma,obj.psi,res);

            Jwk = W.*Jk;
            %Jwk = Jk;
            J = sqrt(sum(abs(Jwk).^2));
        end
        
        function runOptimization(obj,res,verbose)
            if(nargin<3)
                verbose = true;
            end
            % Numerical optimization
            if(verbose)
                fmuoptionsn = optimoptions(@fminunc,...
                          'TolFun',1e-10,'TolX',1e-10,...
                          'Display','iter-detailed','Algorithm','quasi-newton',...
                          'StepTolerance',1e-10,...
                          'MaxFunctionEvaluations',3500,...
                          'SpecifyObjectiveGradient',false,...
                          'DerivativeCheck','off');
            else
                fmuoptionsn = optimoptions(@fminunc,...
                          'TolFun',1e-10,'TolX',1e-10,...
                          'Display','off','Algorithm','quasi-newton',...
                          'StepTolerance',1e-10,...
                          'MaxFunctionEvaluations',3500,...
                          'SpecifyObjectiveGradient',false,...
                          'DerivativeCheck','off');
            end
            
            fJ = @(rho) obj.evalError(rho,res);
            obj.rho_o = fminunc(fJ,obj.rho_i,fmuoptionsn);
            if(verbose)
                [J,Ai,Pi,Oi] = fJ(obj.rho_o)
            end
        end
        
        function continueOptimization(obj,res,verbose)
            if(nargin<3)
                verbose = true;
            end
            % Numerical optimization
            if(verbose)
                fmuoptionsn = optimoptions(@fminunc,...
                          'TolFun',1e-10,'TolX',1e-10,...
                          'Display','iter-detailed','Algorithm','quasi-newton',...
                          'StepTolerance',1e-10,...
                          'MaxFunctionEvaluations',3500,...
                          'SpecifyObjectiveGradient',false,...
                          'DerivativeCheck','off');
            else
                fmuoptionsn = optimoptions(@fminunc,...
                          'TolFun',1e-10,'TolX',1e-10,...
                          'Display','off','Algorithm','quasi-newton',...
                          'StepTolerance',1e-10,...
                          'MaxFunctionEvaluations',3500,...
                          'SpecifyObjectiveGradient',false,...
                          'DerivativeCheck','off');
            end
            
            fJ = @(rho) obj.evalError(rho,res);
            obj.rho_o = fminunc(fJ,obj.rho_o,fmuoptionsn);
            if(verbose)
                [J,Ai,Pi,Oi] = fJ(obj.rho_o)
            end
        end
        
        function setOrderToZero(obj,r)
            obj.A_contr(:,r) = true;
            obj.A_fixed(:,r) = true;
            obj.A(:,r) = 0;
            
            obj.P_contr(:,r) = true;
            obj.P_fixed(:,r) = true;
            obj.P(:,r) = 0;
        end
        
        function setRhoInit(obj,rho_i)
            obj.rho_i = rho_i;
        end
        
        function setMetSync(obj,n,nref,pdiff)
            for ri = 1:obj.r
                obj.setMetAmpRef(n,ri,nref);
                obj.setMetPhaseRef(n,ri,nref,pdiff);
            end
        end
        
        function setPlatAmpToInit(obj,n,r,a)
            if(n>3)
                error('Simulation:configurationError','Only 3 platform DOFs');
            end
            obj.A_contr(n,r) = false;
            obj.A_fixed(n,r) = false;
            obj.A(n,r) = a;
        end
        
        function setPlatPhaseToInit(obj,n,r,p)
            if(n>3)
                error('Simulation:configurationError','Only 3 platform DOFs');
            end
            obj.P_contr(n,r) = false;
            obj.P_fixed(n,r) = false;
            obj.P(n,r) = p;
        end
        
        function setMetAmpToInit(obj,n,r,a)
            if(n>obj.N)
                error('Simulation:configurationError','Only N mets');
            end
            obj.A_contr(3+n,r) = false;
            obj.A_fixed(3+n,r) = false;
            obj.A(3+n,r) = a;
        end
        
        function setOToInit(obj,o)
            obj.O_contr = false;
            obj.O_fixed = false;
            obj.O = o;
        end
        
        function setMetPhaseToInit(obj,n,r,p)
            if(n>obj.N)
                error('Simulation:configurationError','Only N mets');
            end
            obj.P_contr(3+n,r) = false;
            obj.P_fixed(3+n,r) = false;
            obj.P(3+n,r) = p;
        end
        
        function setPlatToZero(obj,n,r)
            obj.setPlatAmpToZero(n,r);
            obj.setPlatPhaseToZero(n,r);
        end
        
        function setMetToZero(obj,n,r)
            obj.setMetAmpToZero(n,r);
            obj.setMetPhaseToZero(n,r);
        end
        
%         function setPlatAmpRef(obj,n,r,nref)
%             if(n>3)
%                 error('Simulation:configurationError','Only 3 platform DOFs');
%             end
%             if((obj.A_contr(n,r)== true) && (obj.A_fixed== false))
%                 error('Simulation:configurationError','Reference already a reference');
%             end
%             obj.A_contr(n,r) = true;
%             obj.A_fixed(n,r) = false;
%             obj.A(n,r) = nref;
%         end
%         
%         function setPlatPhaseRef(obj,n,r,nref,pdiff)
%             if(n>3)
%                 error('Simulation:configurationError','Only 3 platform DOFs');
%             end
%             if((obj.P_contr(n,r)== true) && (obj.P_fixed== false))
%                 error('Simulation:configurationError','Reference already a reference');
%             end
%             obj.P_contr(n,r) = true;
%             obj.P_fixed(n,r) = false;
%             obj.P(n,r) = nref;
%             obj.Pshift(n,r) = pdiff;
%         end
        
        function setMetAmpRef(obj,n,r,nref)
            if(n>obj.N)
                error('Simulation:configurationError','Only N mets');
            end
            obj.A_contr(3+n,r) = true;
            obj.A_fixed(3+n,r) = false;
            obj.A_ref(3+n,r) = 3+nref;
            obj.A(3+n,r) = NaN;
        end
        
        function setMetPhaseRef(obj,n,r,nref,pdiff)
            if(n>obj.N)
                error('Simulation:configurationError','Only N mets');
            end
            obj.P_contr(3+n,r) = true;
            obj.P_fixed(3+n,r) = false;
            obj.P(3+n,r) = NaN;
            obj.P_ref(3+n,r) = 3+nref;
            obj.Pshift(3+n,r) = pdiff;
        end
        
        function setOToFixed(obj,o)
            obj.O_contr = true;
            obj.O_fixed = true;
            obj.O = o;
        end
        
        function setPlatAmpToFixed(obj,n,r,a)
            if(n>3)
                error('Simulation:configurationError','Only 3 platform DOFs');
            end
            obj.A_contr(n,r) = true;
            obj.A_fixed(n,r) = true;
            obj.A(n,r) = a;
        end
        
        function setPlatPhaseToFixed(obj,n,r,p)
            if(n>3)
                error('Simulation:configurationError','Only 3 platform DOFs');
            end
            obj.P_contr(n,r) = true;
            obj.P_fixed(n,r) = true;
            obj.P(n,r) = p;
        end
        
        function setMetAmpToFixed(obj,n,r,a)
            if(n>obj.N)
                error('Simulation:configurationError','Only N mets');
            end
            obj.A_contr(3+n,r) = true;
            obj.A_fixed(3+n,r) = true;
            obj.A(3+n,r) = a;
        end
        
        function setMetPhaseToFixed(obj,n,r,p)
            if(n>obj.N)
                error('Simulation:configurationError','Only N mets');
            end
            obj.P_contr(3+n,r) = true;
            obj.P_fixed(3+n,r) = true;
            obj.P(3+n,r) = p;
        end
        
        function setPlatAmpToZero(obj,n,r)
            obj.setPlatAmpToFixed(n,r,0);
        end
        
        function setPlatPhaseToZero(obj,n,r)
            obj.setPlatPhaseToFixed(n,r,0);
        end
        
        function setMetAmpToZero(obj,n,r)
            obj.setMetAmpToFixed(n,r,0);
        end
        
        function setMetPhaseToZero(obj,n,r)
            obj.setMetPhaseToFixed(n,r,0);
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
            s = obj; 
            save(strcat(folder,fPath),'s','-v7.3');
        end
        
        
        
    end
    
    
    methods(Static)
        
        
                
    end
end


%simulate_1_x - Simulation to find reference and explore initial describing
% function analysis
%
% Author: M.M. Hoogesteger - Msc student of Mechanical Engineering
% Eindhoven University of Technology, Mechanical Engineering, Dynamics and
% control
% email address: m.m.hoogesteger@student.tue.nl
% May 2017; Last revision: 16-July-2018
function [s] = auto_simulation(simSettings)
fprintf(['Starting job ' simSettings.name '.\n'])
% Generate loop mat
if(simSettings.hasVars)
    hasVars = true;
    varNames = simSettings.varNameCell;
    varVecs = simSettings.varVecCell;
    
    n_vars = numel(varNames);
    n_els = zeros(n_vars,1);
    for nidx = 1:n_vars
        n_els(nidx) = numel(varVecs{nidx});
    end
    n_els_tot = prod(n_els);
    
    idxsMat = zeros(n_els_tot,n_vars+1);
    idxsMat(:,1) = 1:n_els_tot;
    for nidx = 1:n_vars-1
        idxsMat(:,nidx+1) = mod(floor((idxsMat(:,1)-1)/prod(n_els(nidx+1:end))),n_els(nidx))+1;
    end
    idxsMat(:,n_vars+1) = mod((idxsMat(:,1)-1),n_els(n_vars))+1;
    nSims = n_els_tot;
else
    hasVars = false;    
    nSims = 1;
end




% Prepare  and sim


saveFolder = ['../data/' simSettings.name '/'];


Svec = cell(1,nSims);
for simIdx = 1:nSims
    fprintf(['Starting simulation ' num2str(simIdx) ' of ' num2str(nSims) '.\n'])
    % Load variable value into simSettings
    if(hasVars)
        for nidx = 1: n_vars
            simSettings.(varNames{nidx}) = varVecs{nidx}{idxsMat(simIdx,nidx+1)};
        end
    end
    if(isfield(simSettings,'overridePsi'))
        simSettings.psi = simSettings.gamma + simSettings.overridePsi;
    end

    s = Simulation(['sim_' num2str(simIdx)]);
    s.setModel(simSettings.('N'),simSettings.('x'),simSettings.('y'),simSettings.('phi'));

    % Prepare configuration
    R = simSettings.('R');
    gamma = simSettings.('gamma');
    psi = simSettings.('psi');

    s.setConfiguration(R,gamma,psi);
    s.generateModel;
    %%
    % Initial conditions and other settings
    
    s.setSimParams(simSettings.('fs'),simSettings.('Tend'))
    if(simSettings.escapement=="Original")
        s.setEscapement(simSettings.('escapement'))
    elseif(simSettings.('escapement')=="Hamiltonian")
        s.setEscapement(simSettings.('escapement'),...
                        simSettings.('Hstars'),...
                        simSettings.('Ugains'));
    end
    s.setInitialConditionsZeroVelocity(simSettings.('nReps'),simSettings.('xyphi0'),simSettings.('Ntheta0'));

    s.simulate(false);
    %%
    s.analyze(simSettings.('Tss'),false);

    s.save(saveFolder);

    Svec{simIdx} = s.name;
    
end
%%
simSet.Svec = Svec;
if(simSettings.hasVars)
    simSet.idxsMat = idxsMat;
    simSet.n_vars = n_vars;
end
simSet.simSettings = simSettings;


save([saveFolder 'simSet.mat'],'simSet');


    



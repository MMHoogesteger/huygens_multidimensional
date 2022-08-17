disp('Adding library folders')
if(~exist('bPrepEnv','var'))
    addpath(genpath('../lib/simulations/'));
    addpath(genpath('../lib/image_app/'));
    addpath(genpath('../lib/describing/'));
    addpath(genpath('../lib/matlab2tikz/'));
    addpath(genpath('../lib/model_generation/'));
    addpath(genpath('../lib/mtimesx/'));
    addpath(genpath('../models/'));
    addpath(genpath('../scripts/figurescripts'));
end
bPrepEnv = true;
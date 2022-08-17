%%
prepEnv;
listofJobs = dir('jobs/*.mat');
for jobIdx =  1:numel(listofJobs)
    job = listofJobs(jobIdx);
    load([job.folder '/' job.name]);
    
    auto_simulation(simSettings);
end



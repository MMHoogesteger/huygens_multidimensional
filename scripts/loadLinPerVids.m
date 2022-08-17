clear rv
disp('Loading Linper vid results')
vidFolder = 'H:\Vids\Vids_linper_two\analyzed\';
load([vidFolder '\video' num2str(0) '_results.mat']);
rv(1) = results;
load([vidFolder '\video' num2str(1) '_results.mat']);
rv(2) = results;
load([vidFolder '\video' num2str(2) '_results.mat']);
rv(3) = results;
load([vidFolder '\video' num2str(4) '_results.mat']);
rv(4) = results;
load([vidFolder '\video' num2str(3) '_results.mat']);
rv(5) = results;
load([vidFolder '\video' num2str(5) '_results.mat']);
rv(6) = results;
load([vidFolder '\video' num2str(6) '_results.mat']);
rv(7) = results;
disp('Finished loading')
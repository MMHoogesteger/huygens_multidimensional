vids = {'video0.avi','video1.avi','video2.avi','video3.avi','video4.avi','video5.avi','video6.avi'};

path = '..\vids';
%path = 'vids';
delete(gcp('nocreate'))



%%
nWorkers = 12;
nVids    = numel(vids);
nZeros   = nWorkers - mod(nVids,nWorkers);
ids = [1:nVids zeros(1,nZeros)];
ids = reshape(ids,nWorkers,[]).';
nLoops = size(ids,1);
pFunc = @(f,nf,tEst) fprintf('Progress for: frame %d of %d\nTime Remaining(min): %d\n',f,nf,tEst)
%%
parpool(nWorkers)
pctRunOnAll warning('off','vision:calibrate:boardShouldBeAsymmetric')
spmd(nWorkers) 
    for id = 1:nLoops
        idx = ids(id,labindex);
        if(idx~=0)
            vidName = vids{idx};
            fprintf('Starting %s\n',vidName)

            v = VideoProcessor();
            v.updateVideoFile(vidName,path);

            a= v.analysis;

            
            a.prepareAnalysis();
            %a.stepAnalysis();
            a.runAnalysis(pFunc)
            v.saveAnalysis();
            v.saveResults();
            fprintf('Finished %s\n',vidName)
        end
    end
end
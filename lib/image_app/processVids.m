vids = {'video49.avi','video52.avi','video55.avi','video60.avi','video61.avi'};
path = '/mnt/hit/mhoogesteger/Matlab/vids';
path = 'vids';
delete(gcp('nocreate'))
parpool;
pctRunOnAll warning('off','vision:calibrate:boardShouldBeAsymmetric')
parfor vId = 1:numel(vids)
    vidName = vids{vId};
    fprintf('Starting %s\n',vidName)
    
    v = VideoProcessor();
    v.updateVideoFile(vidName,path);

    a= v.analysis;

    pFunc = @(f,nf,tEst) fprintf('Progress for %s: frame %d of %d\nTime Remaining(min): %d\n',vids{vId},f,nf,tEst)
    a.prepareAnalysis();
    %a.stepAnalysis();
    a.runAnalysis(pFunc)
    v.saveAnalysis();
    v.saveResults();
    fprintf('Finished %s\n',vids{vId})
end
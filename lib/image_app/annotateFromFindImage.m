function [annotatedIm] = annotateFromFindImage(im,metaDataCam)
%annotateFromFindImage(im,metaDataCam) Annotate image based on detection limits in metaDataCam
% Without relying on other data, annotate the image based on metaDataCam
% using the relevant functions.
%
% See also: findCheckers, findMarkers
colLims = metaDataCam.colLims;
blobLims = metaDataCam.blobLims;
pcolLims = metaDataCam.pcolLims;
pblobLims = metaDataCam.pblobLims;


checkers = findCheckers(im);

for idx = 1:checkers(1).n
    checker = checkers(idx);
    px = min(checker.imagePoints(:,1));
    py = min(checker.imagePoints(:,2));
    wx = max(checker.imagePoints(:,1)-px);
    wy = max(checker.imagePoints(:,2)-py);
    im = insertObjectAnnotation(im, 'rectangle',...
        [px,py,wx,wy],['C' num2str(idx)]);
end

mMarkers = findMarkers(im,colLims,blobLims);

for idx = 1:mMarkers(1).n
    marker = mMarkers(idx);
    im = insertObjectAnnotation(im, 'rectangle',...
        marker.box.', ['M' num2str(idx)]);
end

pMarkers = findMarkers(im,pcolLims,pblobLims);

for idx = 1:pMarkers(1).n
    marker = pMarkers(idx);
    im = insertObjectAnnotation(im, 'rectangle',...
        marker.box.', ['P' num2str(idx)]);
end
annotatedIm = im;
end


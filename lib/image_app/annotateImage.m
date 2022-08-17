function [annotatedIm] = annotateImage(im,mCam)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
for idx = 1:mCam.checkers(1).n
    checker = mCam.checkers(idx);
    px = min(checker.imagePoints(:,1));
    py = min(checker.imagePoints(:,2));
    wx = max(checker.imagePoints(:,1)-px);
    wy = max(checker.imagePoints(:,2)-py);
    im = insertObjectAnnotation(im, 'rectangle',...
        [px,py,wx,wy],['C' num2str(idx)]);
end

for idx = 1:mCam.mMarkers(1).n
    marker = mCam.mMarkers(idx);
    im = insertObjectAnnotation(im, 'rectangle',...
        marker.box.', ['M' num2str(idx)]);
end

for idx = 1:mCam.pMarkers(1).n
    marker = mCam.pMarkers(idx);
    im = insertObjectAnnotation(im, 'rectangle',...
        marker.box.', ['P' num2str(idx)]);
end
annotatedIm = im;
end


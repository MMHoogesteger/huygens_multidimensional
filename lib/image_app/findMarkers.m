function [blobs,markedImage] = findMarkers(image,colLims,blobLims)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

ArMin = blobLims(1);
ArMax = blobLims(2);
PiMin = blobLims(3);
PiMax = blobLims(4);
MajMin = blobLims(5);
MajMax = blobLims(6);


[imBW,markedImage] = maskImage(image,colLims);
% Find connected components.
blobAnalysis = vision.BlobAnalysis('AreaOutputPort', true,...
    'CentroidOutputPort', true,...
    'BoundingBoxOutputPort', true,...
    'MajorAxisLengthOutputPort', true,...
    'MinimumBlobArea', 3, 'ExcludeBorderBlobs', false);
[areas, centers, boxes, majors] = step(blobAnalysis, imBW);

areas = double(areas);
% Sort
[areas, idx] = sort(areas, 'Descend');
centers = centers(idx,:).';
boxes = boxes(idx,:).';
majors = majors(idx,:);

% Pi
pis = 4*areas./(majors.^2);

%Limits
idxs = find(areas>ArMin & areas<ArMax &...
            pis>PiMin & pis<PiMax &...
            majors>MajMin & majors< MajMax);

% Struct
blobs = struct;
blobs.n = length(idxs);
for idx = 1:length(idxs)
    blobs(idx).id = idx;
    blobs(idx).area = areas(idxs(idx));
    blobs(idx).x= centers(1,idxs(idx));
    blobs(idx).y= centers(2,idxs(idx));
    blobs(idx).center = [blobs(idx).x blobs(idx).y];
    blobs(idx).box = boxes(:,idxs(idx));
    blobs(idx).major = majors(idxs(idx));
    blobs(idx).pi = pis(idxs(idx));
end


% Insert labels for the coins.
markedImage = insertObjectAnnotation(markedImage, 'rectangle',...
    boxes(:,idxs).', 'marker');



end


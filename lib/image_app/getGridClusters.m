function [] = getGridClusters(points,gridResolution,gridOverlap)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
g = struct;
g.resolution = gridResolution;
g.overlap = gridOverlap;

g.se = strel('disk',10);
g.sed = strel('disk',3);

g.xList = points(:,1);
g.yList = points(:,2);

% Grid boundaries
g.minX = min(g.xList);
g.minY = min(g.yList);
g.maxX = max(g.xList);
g.maxY = max(g.yList);

g.gridX = g.minX:g.resolution:g.maxX+g.resolution;
g.gridY = g.minY:g.resolution:g.maxY+g.resolution;


g.gridIm = zeros(length(g.gridY)-1-g.overlap,...
                 length(g.gridX)-1-g.overlap);
g.gridPointList = cell(length(g.gridY)-1-g.overlap,...
                       length(g.gridX)-1-g.overlap);
for idX = 1:length(g.gridX)-1-g.overlap
    xMin = g.gridX(idX);
    xMax = g.gridX(idX+1+g.overlap);
    for idY = 1:length(g.gridY)-1-g.overlap
        yMin = g.gridY(idY);
        yMax = g.gridY(idY+1+g.overlap);

        idxs = find(g.xList>=xMin & g.xList<xMax...
                  & g.yList>=yMin & g.yList<yMax);
        g.gridIm(idY,idX) = length(idxs);
        g.gridPointList{idY,idX} = idxs;
    end
end

gridIm = g.gridIm./max(max(g.gridIm));

figure;
imshow(gridIm)
level = graythresh(gridIm);
        BW = imbinarize(gridIm,level);
        figure, imshow(BW)


closeBW = imdilate(imclose(gridIm,g.se),g.sed);
figure, imshow(closeBW)
D = bwdist(~closeBW);
figure, imshow(D,[])
D = -D;
D(~closeBW) = Inf;
figure, imshow(D,[])
D = imhmin(D,20);
figure, imshow(D,[])
g.L = watershed(D);
g.L(~closeBW) = 0;
rgb = label2rgb(g.L,'jet',[.5 .5 .5]);
figure
imshow(rgb,'InitialMagnification','fit')
title('Watershed transform of D')
end


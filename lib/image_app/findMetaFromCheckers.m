function [metaFromCheckers] = findMetaFromCheckers(camIm,cameraParams)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here


checkers = findCheckers(camIm);

metaFromCheckers = struct;


%%
c = numel(checkers);
if(c<1)
    disp('Cannot use camera, no checkerboard');
else
    squareSize = 9.85; % in millimeters
    
    for idx = 1:c
        worldPoints = generateCheckerboardPoints(checkers(idx).boardSize, squareSize);
        [R, t] = extrinsics(checkers(idx).imagePoints, worldPoints, cameraParams);
        checkers(idx).R = R;
        checkers(idx).t = t;
        checkers(idx).PosW = -t*R.';
        checkers(idx).n = c;
        
    end
    metaFromCheckers.R = checkers(1).R;
    metaFromCheckers.t = checkers(1).t;
    
    camPosW = [];
    camPosZ = [];
    for idx = 1:c
        checkers(idx).worldPoints = pointsToWorld(cameraParams,...
            metaFromCheckers.R,metaFromCheckers.t,checkers(idx).imagePoints);
        checkers(idx).xVec = checkers(idx).worldPoints(checkers(idx).boardSize(1),:)...
            -checkers(idx).worldPoints(1,:);
        checkers(idx).yVec = checkers(idx).worldPoints(2,:)...
            -checkers(idx).worldPoints(1,:);
        checkers(idx).Or = atan2(checkers(idx).yVec(2),checkers(idx).yVec(1));
        checkers(idx).xVec = checkers(idx).xVec./norm(checkers(idx).xVec);
        checkers(idx).yVec = checkers(idx).yVec./norm(checkers(idx).yVec);
        checkers(idx).camPosW = checkers(idx).PosW(1:2)*[checkers(idx).xVec; checkers(idx).yVec]...
            +checkers(idx).worldPoints(1,:);
        checkers(idx).camPosZ = checkers(idx).PosW(3);
        camPosW = [camPosW;checkers(idx).camPosW];
        camPosZ = [camPosZ;checkers(idx).camPosZ];
    end
    camPosWMean = mean(camPosW,1);
    camPosWDiff = max(abs(camPosW-camPosWMean),1);
    metaFromCheckers.camPosWMean = camPosWMean;
    metaFromCheckers.camPosWDiff = camPosWDiff;
    
    camPosZMean = mean(camPosZ);
    camPosZDiff = mean(abs(camPosZ-camPosZMean));
    metaFromCheckers.camPosZMean = camPosZMean;
    metaFromCheckers.camPosZDiff = camPosZDiff;
    
end
metaFromCheckers.checkers = checkers;
end


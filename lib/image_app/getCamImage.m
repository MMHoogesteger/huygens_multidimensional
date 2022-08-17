function [camIm,newOrigin] = getCamImage(im,m,camId)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Split image
horVec = 0:m.video.nCamsHor-1;
verVec = 0:m.video.nCamsVer-1;

[camHorPos,camVerPos] = meshgrid(horVec,verVec);

horPos = camHorPos(camId);
verPos = camVerPos(camId);
camIm = im(verPos*480+1:(verPos+1)*480,horPos*640+1:(horPos+1)*640,:);
   
[camIm, newOrigin] = undistortImage(camIm, m.cam(camId).parameters, 'OutputView', 'full');
end


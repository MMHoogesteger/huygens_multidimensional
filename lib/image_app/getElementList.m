function [elementList] = getElementList(m,exCam,noCheckers)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

camList = setdiff(1:m.nCams,exCam);
elementList = [];
for idxCam = camList
    mCam = m.cam(idxCam);
    if(~noCheckers)
        for idx = 1:mCam.checkers(1).n
            elementList = [elementList;[num2str(idxCam) 'C' num2str(idx)]];
        end
    end

    for idx = 1:mCam.mMarkers(1).n
        elementList = [elementList;[num2str(idxCam) 'M' num2str(idx)]];
    end

    for idx = 1:mCam.pMarkers(1).n
        elementList = [elementList;[num2str(idxCam) 'P' num2str(idx)]];
    end
end
end


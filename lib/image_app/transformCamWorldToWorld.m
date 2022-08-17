function [pointsW] = transformCamWorldToWorld(metaData,pointsC,camID)

pointsW = metaData.cam(camID).tW+pointsC*metaData.cam(camID).RW;
end


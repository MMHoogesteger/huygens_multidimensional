function [pointsW] = transformCamToWorld(metaData,pointsC,camID)
pointsCW = transformCamToCamWorld(metaData,pointsC,camID);
pointsW = transformCamWorldToWorld(metaData,pointsCW,camID);
end


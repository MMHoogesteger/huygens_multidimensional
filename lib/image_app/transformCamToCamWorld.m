function [pointsW] = transformCamToCamWorld(metaData,pointsC,camID)
mCam = metaData.cam(camID);
pointsW = pointsToWorld(mCam.parameters, mCam.R, mCam.t, pointsC);
end


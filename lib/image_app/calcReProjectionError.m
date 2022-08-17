function [E,Mfx,Mfy,Mfrpx,Mfrpy] = calcReProjectionError(framesList,projectedPoints,RStor,tStor,camCenter,psi,H,L,camZ,Mpos)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here


[nIdx,~] = find(~isnan(projectedPoints));

tVec = repmat((-70:0.5:70)*2*pi/360,length(framesList),1);
rVec = L*sin(tVec);
zVec = H + L*cos(tVec);

cCx = zeros(length(framesList),1);
cCy = cCx;
Mfx = zeros(length(framesList),1);
Mfy = zeros(length(framesList),1);
Mfrpx = zeros(length(framesList),1);
Mfrpy = zeros(length(framesList),1);

for idx = 1:length(framesList)
    R = squeeze(RStor(framesList(idx),:,:));
    t = tStor(framesList(idx),:);

    cC = (camCenter.'-t)*R.';
    cCx(idx,1) = cC(1);
    cCy(idx,1) = cC(2);
end

Mpx = projectedPoints(:,1);
Mpy = projectedPoints(:,2);



dMx = rVec.*cos(psi);
dMy = rVec.*sin(psi);

Mx = Mpos(1) + dMx;
My = Mpos(2) + dMy;

Mrpx = Mx + (Mx-cCx).*(zVec./(camZ-zVec));
Mrpy = My + (My-cCy).*(zVec./(camZ-zVec));

distx = (Mpx-Mrpx).^2 + (Mpy-Mrpy).^2;
[mdist,midx] = min(distx,[],2);
for idx = 1:length(midx)
    Mfx(idx) = Mx(idx,midx(idx));
    Mfy(idx) = My(idx,midx(idx));
    Mfrpx(idx) = Mrpx(idx,midx(idx));
    Mfrpy(idx) = Mrpy(idx,midx(idx));
end
E = sum(mdist);



end


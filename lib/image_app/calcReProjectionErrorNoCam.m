function [E,Mfx,Mfy,Mfrpx,Mfrpy,Mtheta] = calcReProjectionErrorNoCam(projectedPoints,camCenters,psi,H,L,camZs,Mpos)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

n = size(projectedPoints,1);

tVec = repmat((-70:0.1:70)*2*pi/360,n,1);
rVec = L*sin(tVec);
zVec = H + L*cos(tVec);

cCx = camCenters(:,1);
cCy = camCenters(:,2);
Mfx = zeros(n,1);
Mfy = zeros(n,1);
Mfrpx = zeros(n,1);
Mfrpy = zeros(n,1);
Mtheta = zeros(n,1);


Mpx = projectedPoints(:,1);
Mpy = projectedPoints(:,2);



dMx = rVec.*cos(psi);
dMy = rVec.*sin(psi);

Mx = Mpos(1) + dMx;
My = Mpos(2) + dMy;

Mrpx = Mx + (Mx-cCx).*(zVec./(camZs-zVec));
Mrpy = My + (My-cCy).*(zVec./(camZs-zVec));

distx = (Mpx-Mrpx).^2 + (Mpy-Mrpy).^2;
[mdist,midx] = min(distx,[],2);
for idx = 1:n
    Mfx(idx) = Mx(idx,midx(idx));
    Mfy(idx) = My(idx,midx(idx));
    Mfrpx(idx) = Mrpx(idx,midx(idx));
    Mfrpy(idx) = Mrpy(idx,midx(idx));
    Mtheta(idx) = tVec(idx,midx(idx));
end
E = sum(mdist);



end


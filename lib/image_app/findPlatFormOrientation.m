function [phi,t,R] = findPlatFormOrientation(cOrT,cPosT,pPosT,m,phi0,t0)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
options =optimset('MaxFunEvals',200000,...
                      'MaxIter',200000,...
                      'TolFun',1e-7,...
                      'TolX',1e-7);
                  
angleW = 1;
posCW = 10;
posPW = 50;

refOr = m.world.cOr;
refCPos = m.world.cPos;
refPPos = m.world.pPos;

fMin = @(x) angleW* calcAngleError(cOrT,refOr,x(1))...
    + posCW* calcCPositionError(cPosT,refCPos,x(1),x(2:3))...
    + posPW* calcPPositionError(pPosT,refPPos,x(1),x(2:3));
x = fminsearch(@(x) fMin(x),[phi0 t0],options);
phi = x(1);
R = [cos(phi) sin(phi); -sin(phi) cos(phi)];
t = x(2:3);


end

 function dphi = calcAngleError(cOrT,refOr,phi)
    dphi = 0;
    for idx = 1:numel(cOrT)
        for cIdx = 1:numel(cOrT{idx})
            dphi = dphi + min((cOrT{idx}(cIdx)-refOr{idx} -phi).^2);
        end
    end
end


function dist = calcCPositionError(cPosT,refPos,phi,t)
    R = [cos(phi) sin(phi); -sin(phi) cos(phi)];
    dist = 0;
    for idx = 1:numel(cPosT)
        cPosTX = cPosT{idx}(:,1);
        cPosTY = cPosT{idx}(:,2);
        for cIdx = 1:numel(cPosTX)
            dist = dist + min(sum((([cPosTX(cIdx) cPosTY(cIdx)]-t)*R-refPos{idx}).^2,2));
        end
    end
end

function dist = calcPPositionError(pPosT,refPPos,phi,t)
    R = [cos(phi) sin(phi); -sin(phi) cos(phi)];
    dist = 0;
    for idx = 1:numel(pPosT)
        pPosTX = pPosT{idx}(:,1);
        pPosTY = pPosT{idx}(:,2);
        for pIdx = 1:numel(pPosTX)
            dist = dist + min(sum((([pPosTX(pIdx) pPosTY(pIdx)]-t)*R-refPPos{idx}).^2,2));
        end
    end
end

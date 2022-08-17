function [world] = assembleWorld(metaData,H)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%%
figure
plot(0,0,'kx')
hold on

world.checkers = {};
world.pMarkers = {};
world.mMarkers = {};
world.camCenters = {};
world.camDiffs = {};
world.camZs = {};



for camIdx = 1:metaData.nCams
    world.checkers{camIdx} = {};
    world.pMarkers{camIdx} = {};
    world.mMarkers{camIdx} = {};
    
    meta = metaData.cam(camIdx);
    camZ = meta.camPosZMean;
    camCenter = transformCamWorldToWorld(metaData,meta.camPosWMean,camIdx);
    world.camCenters = [world.camCenters {camCenter}];
    world.camDiffs = [world.camDiffs {meta.camPosWDiff}];
    world.camZs = [world.camZs {camZ}];
    
    
    for cIdx = 1:meta.checkers(1).n
        c = struct;
        checker = meta.checkers(cIdx);
        c.boardSize = checker.boardSize;
        c.points = transformCamToWorld(metaData,checker.imagePoints,camIdx);
        world.checkers{camIdx} = [world.checkers{camIdx} {c}];
    end
            
    for pIdx = 1:meta.pMarkers(1).n
        m = struct;
        marker = meta.pMarkers(pIdx);
        m.center = transformCamToWorld(metaData,marker.center,camIdx);
        world.pMarkers{camIdx} = [world.pMarkers{camIdx} {m}];
    end
    
    for mIdx = 1:meta.mMarkers(1).n
        m = struct;
        marker = meta.mMarkers(mIdx);
        m.pcenter = transformCamToWorld(metaData,marker.center,camIdx);
        camRay = camCenter-m.pcenter;
        m.rcenter = m.pcenter + camRay*(-H/camZ);
        world.mMarkers{camIdx} = [world.mMarkers{camIdx} {m}];
    end
    
    world.nMMarkers{camIdx} = length(world.mMarkers{camIdx});
    world.nPMarkers{camIdx} = length(world.pMarkers{camIdx});
    world.nCheckers{camIdx} = length(world.checkers{camIdx});

    %%
    world.mPos{camIdx} = zeros(world.nMMarkers{camIdx},2);
    world.pPos{camIdx} = zeros(world.nPMarkers{camIdx},2);
    world.cPos{camIdx} = zeros(world.nCheckers{camIdx},2);
    world.cOr{camIdx} = zeros(world.nCheckers{camIdx},1);

    for cid = 1:world.nCheckers{camIdx}
        world.cPos{camIdx}(cid,:) = mean(world.checkers{camIdx}{cid}.points);
        yVec = world.checkers{camIdx}{cid}.points(2,:)-world.checkers{camIdx}{cid}.points(1,:);
        world.cOr{camIdx}(cid,:) = atan2(yVec(2),yVec(1));
    end

    for mid = 1:world.nMMarkers{camIdx}
        world.mPos{camIdx}(mid,:) = world.mMarkers{camIdx}{mid}.pcenter;
    end

    for pid = 1:world.nPMarkers{camIdx}
        world.pPos{camIdx}(pid,:) = world.pMarkers{camIdx}{pid}.center;
    end
    
    for f = 1:world.nCheckers{camIdx}
        plot(world.checkers{camIdx}{f}.points(:,1),world.checkers{camIdx}{f}.points(:,2),'rx')

    end
    for f = 1:world.nMMarkers{camIdx}
        plot(world.mMarkers{camIdx}{f}.pcenter(1),world.mMarkers{camIdx}{f}.pcenter(2),'bo')
        plot(world.mMarkers{camIdx}{f}.rcenter(1),world.mMarkers{camIdx}{f}.rcenter(2),'bx')
    end
    
    ax = gca();
    for f = 1:length(world.camCenters)

        ax.ColorOrderIndex = f;
        errorbar(world.camCenters{f}(1),world.camCenters{f}(2),...
        world.camDiffs{f}(2),world.camDiffs{f}(2),...
        world.camDiffs{f}(1),world.camDiffs{f}(1),'s')
        %plot(world.camCenters{f}(1),world.camCenters{f}(2),'s')
    end
    for f = 1:length(world.camCenters)
        ax.ColorOrderIndex = f;
        p =transformCamToWorld(metaData,[0,0;640 0; 640 480 ; 0 480; 0 0],f);
        plot(p(:,1),p(:,2),':')
    end
    
end





end


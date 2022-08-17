function [E] = calcErrorConstraints(metaData,constraints,X0)




phis = X0(:,1);
ts = X0(:,[2 3])*1e2;

for phiidx = 1:length(phis)
    phi = phis(phiidx);
    Rs{phiidx} = [cos(phi) sin(phi); -sin(phi) cos(phi)];
end



d=zeros(1,length(constraints));
for cidx = 1:length(constraints)
    con = constraints{cidx};
    switch con{1}
        case 'p'
            pt = [con{3};con{4}];
            [cam,pos] = parseConstraint(con{2});
            posW = transformCamToCamWorld(metaData,pos,cam);
            pr = ts(cam,:).'+Rs{cam}.'*mean(posW,1).';
            distance = pr-pt;            
            d(cidx) = sum(distance.^2)*con{5};
            
        case 'm'
            [cam,pos] = parseConstraint(con{2});
            posW = transformCamToCamWorld(metaData,pos,cam);
            p1 = ts(cam,:).'+Rs{cam}.'*posW.';
            [cam,pos] = parseConstraint(con{3});
            posW = transformCamToCamWorld(metaData,pos,cam);
            p2 = ts(cam,:).'+Rs{cam}.'*posW.';
            distance = p2-p1;            
            d(cidx) = sum(sum(distance.^2))*con{4};
            
        case 'l'
            [cam,pos] = parseConstraint(con{2});
            posW = transformCamToCamWorld(metaData,pos,cam);
            p1 = ts(cam,:).'+Rs{cam}.'*posW.';
            [cam,pos] = parseConstraint(con{3});
            posW = transformCamToCamWorld(metaData,pos,cam);
            p2 = ts(cam,:).'+Rs{cam}.'*posW.';
            distance = p2-p1;
            mag = sqrt(sum(distance.^2));
            switch con{4}
                case 'X+'
                    eunit = mag*[1; 0];
                case 'X-'
                    eunit = mag*[-1; 0];
                case 'Y+'
                    eunit = mag*[0; 1];
                case 'Y-'
                    eunit = mag*[0; -1];
            end
                    
            d(cidx) = sum(sum((distance-eunit).^2))*con{5};
    end
            
            
end

E = sum(d);

    function [camI,pos] = parseConstraint(consStr)
    camI = str2double(consStr(1));
    typeI = consStr(2);
    objI = str2double(consStr(3));
        switch typeI
            case 'C'
                pos = metaData.cam(camI).checkers(objI).imagePoints;
            case 'M'
                pos = [metaData.cam(camI).mMarkers(objI).x,...
                    metaData.cam(camI).mMarkers(objI).y];
            case 'P'
                pos = [metaData.cam(camI).pMarkers(objI).x,...
                    metaData.cam(camI).pMarkers(objI).y];
        end
    end

    



end


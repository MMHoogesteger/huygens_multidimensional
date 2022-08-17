function [h] = optimSet3_3_dist(h)
if(~(h.N==3))
    error('not three');
end

%% Optimization search with in phase
h.setOrder(3);
h.setOrderToZero(2);

h.setPlatToZero(3,1);

h.setMetSync(2,1,2*pi/3*1)
h.setMetSync(3,1,2*pi/3*2)

% h.setMetSync(4,1,0)
% h.setMetSync(5,1,0)
% h.setMetSync(6,1,0)

h.setPlatAmpToInit(1,1,0.001);
h.setPlatAmpToInit(2,1,0.001);
h.setPlatAmpToInit(1,3,0.00002);
h.setPlatAmpToInit(2,3,0.00002);
h.setPlatAmpToInit(3,3,10^(-4));

h.setPlatPhaseToInit(1,1,0.5*pi);
h.setPlatPhaseToInit(2,1,pi);
h.setPlatPhaseToInit(1,3,pi);
h.setPlatPhaseToInit(2,3,pi);
h.setPlatPhaseToInit(3,3,0.2);

h.setMetAmpToInit(1,1,0.72);
h.setMetAmpToInit(1,3,0.01);

h.setMetPhaseToInit(1,3,0.2);
h.setOToInit(10.8);

h.assembleInitial();


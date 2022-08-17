function [h] = optimSet3_3_in(h)
if(~(h.N==3))
    error('not three');
end

%% Optimization search with in phase
h.setOrder(3);
h.setOrderToZero(2);

h.setPlatToZero(1,1);
h.setPlatToZero(2,1);
h.setPlatToZero(1,3);
h.setPlatToZero(2,3);
h.setMetSync(2,1,0)
h.setMetSync(3,1,0)
% h.setMetSync(4,1,0)
% h.setMetSync(5,1,0)
% h.setMetSync(6,1,0)

h.setPlatAmpToInit(3,1,0.01);
h.setPlatAmpToInit(3,3,0.002);
h.setPlatPhaseToInit(3,1,pi);
h.setPlatPhaseToInit(3,3,pi);

h.setMetAmpToInit(1,1,0.72);
h.setMetAmpToInit(1,3,0.01);

h.setMetPhaseToInit(1,3,0.2);
h.setOToInit(10.8);

h.assembleInitial();


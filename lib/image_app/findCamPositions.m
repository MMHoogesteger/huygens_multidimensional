function [X,E] = findCamPositions(metaData,constraints,X0)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

options =optimset('MaxFunEvals',2000*3*metaData.nCams,...
                  'MaxIter',2000*3*metaData.nCams,...
                  'Display','final',...
                  'TolFun',20,...
                  'TolX',1e-5);
X = fminsearch(@(x) calcErrorConstraints(metaData,constraints,x),X0,options);
E = calcErrorConstraints(metaData,constraints,X);
end


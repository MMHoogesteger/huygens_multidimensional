classdef CamData < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        camId uint8;
        parameters;
        paramName char;
        hasParams logical;
        
        calImage uint8;
        newOrigin;
        
        mMarkers;
        pMarkers;
        checkers;
        
        camPosWMean;
        camPosWDiff;
        camPosZMean;
        camPosZDiff;
        R;
        t;
        
        RW;
        tW;
        
        colLims;
        pcolLims;
        blobLims;
        pblobLims;
    end
    
    methods
        function obj = CamData(camId,baseCam)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.hasParams = 0;
            obj.camId = camId;
            obj.parameters = cameraParameters;
            
            obj.colLims = baseCam.colLims;
            obj.pcolLims = baseCam.pcolLims;
            obj.blobLims = baseCam.blobLims;
            obj.pblobLims = baseCam.pblobLims;
        end
        
        function setParameters(obj,parameters,name)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.parameters = parameters;
            obj.paramName = name;
            obj.hasParams = 1;
        end
        
        function setCalImage(obj,image,newOrigin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.calImage = image;
            obj.newOrigin = newOrigin;
        end
        
        function calibrateCam(obj)
            m = findMetaFromCheckers(obj.calImage,obj.parameters);
            
            obj.camPosWMean = m.camPosWMean;
            obj.camPosWDiff = m.camPosWDiff;
            obj.camPosZMean = m.camPosZMean;
            obj.camPosZDiff = m.camPosZDiff;
            
            obj.R = m.R;
            obj.t = m.t;
            
            obj.checkers = m.checkers;
            obj.mMarkers = findMarkers(obj.calImage,obj.colLims,obj.blobLims);
            obj.pMarkers = findMarkers(obj.calImage,obj.pcolLims,obj.pblobLims);
            
        end
        
        function setNewColLims(obj,colLims)
            obj.colLims = colLims;
            obj.calibrateCam();
        end
        
        function setNewPColLims(obj,pcolLims)
            obj.pcolLims = pcolLims;
            obj.calibrateCam();
        end
        
        function setNewBlobLims(obj,blobLims)
            obj.blobLims = blobLims;
            obj.calibrateCam();
        end
        
        function setNewPBlobLims(obj,pblobLims)
            obj.pblobLims = pblobLims;
            obj.calibrateCam();
        end
        
        function updateWorldPos(obj,phi,tW)
            obj.tW = tW;
            obj.RW = [cos(phi) sin(phi); -sin(phi) cos(phi)];
        end
        
    end
end


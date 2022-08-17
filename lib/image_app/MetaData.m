classdef MetaData < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        video;
        calImage;
        
        nCams;
        base;
        
        
        hasLoadedCameraParams;
        loadedCameraParams;
        
        cam CamData;
        vName;
        
        calibration;
        
        camParFinished
        calibrationFinished
        
        world;
    end
    
    methods(Static)
        function basePars = baseMetaParams()
            baseCam = struct;
            basePars = struct;
            baseCam.colLims = [0.278;0.455;0.2;1.0;0.174;1];
            baseCam.pcolLims = [0.911;1.000;0.187;1.0;0.190;1];
            baseCam.blobLims = [50;650;0.7;5;4;60];
            baseCam.pblobLims = [50;650;0.7;5;4;60];
            
            basePars.cam = baseCam;
            basePars.H = 36;
            basePars.L = 70;
        end
    end
    
    methods
        function obj = MetaData(videoProcessor)
            obj.base = obj.baseMetaParams();
            
            % Get properties from parent (only saved here)
            obj.video = videoProcessor.videoProperties;
            obj.calImage = videoProcessor.calImage;
            
            obj.nCams = obj.video.nCams;
            obj.vName = videoProcessor.vName;
            
            for camId  = 1: obj.nCams
                obj.cam(camId) = CamData(camId,obj.base.cam);
            end
            
            obj.hasLoadedCameraParams = 0;
            obj.camParFinished = 0;
            
            
            obj.calibration.X = zeros(obj.nCams,3);
            obj.calibration.constraints = {};
            obj.calibration.constraintsStr = {};
            obj.calibrationFinished = 0;
        end
        
        function state = isFinished(obj)
            %isFinished(obj) Determines if metadata is completed
            %(calibration)
            state = (obj.camParFinished & obj.calibrationFinished);
        end
        
        function loadCameraParameters(obj,cFile,cPath)
            %isFinished(obj) Load camera parameters
            cFullFile = fullfile(cPath,cFile);
            
            if(exist(cFullFile,'file')~=2 )
                error('MetaData:fileNotFound','Camera Parameters file %s not found',cFullFile);
            else                
                c = load(cFullFile);
                obj.hasLoadedCameraParams = 1;
                obj.loadedCameraParams = c.camPam;                
            end
        end
        
        function cParNames = getLoadedParameterNames(obj)
            %isFinished(obj) Returns loaded parameter names in a cell
            sc = struct2cell(obj.loadedCameraParams);
            cParNames = squeeze(sc(1,:,:));
        end
        
        function selectCameraParameters(obj,cameraNum,parNum)
            %isFinished(obj) Selects the parNum parameter for the camNum
            %camera
            p = obj.loadedCameraParams(parNum);
            obj.cam(cameraNum).setParameters(p.cameraParams,p.name);
            
        end
        
        function camParTable = getCamParTable(obj)
            %isFinished(obj) Returns combinations in a table
            camParTable = [{obj.cam.camId}.',{obj.cam.paramName}.'];
        end
        
        function finishCamParameters(obj)
            if(~all([obj.cam.hasParams]))
                error('MetaData:camParNotSelected','Not all cameras have parameters selected');
            end
            
            obj.camParFinished = 1;
            
            for camId = 1:obj.nCams
                [image,newOrigin] = getCamImage(obj.calImage,obj,camId);
                obj.cam(camId).setCalImage(image,newOrigin);
                obj.cam(camId).calibrateCam();
            end
        end
        
        
        
        function setNewColLims(obj,camId,colLims)
            obj.cam(camId).setNewColLims(colLims);
        end
        
        function setNewPColLims(obj,camId,pcolLims)
            obj.cam(camId).setNewPColLims(pcolLims);
        end
        
        function setNewBlobLims(obj,camId,blobLims)
            obj.cam(camId).setNewBlobLims(blobLims);
        end
        
        function setNewPBlobLims(obj,camId,pblobLims)
            obj.cam(camId).setNewPBlobLims(pblobLims);
        end
        
        function setNewCalibration(obj,calibration)
            obj.calibration = calibration;
            obj.calibrationFinished = 1;
            
            for camId = 1:obj.nCams
                obj.cam(camId).updateWorldPos(calibration.X(camId,1),...
                    calibration.X(camId,[2 3])*1e2);
            end
            
            obj.world = assembleWorld(obj,obj.base.H+obj.base.L);
        end
        
    end
end


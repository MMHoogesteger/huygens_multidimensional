classdef VideoProcessor < handle
    % VideoProcessor Controls the processing of a video
    % Detailed explanation goes here
    
    properties
        vName char;  % Video Name
        path char;   % Video Path
        
        % Filenames
        vFile char;  % Video File
        iFile char;  % Image File
        lFile char;  % Log File
        mFile char;  % Metadata File
        aFile char;  % Analysis Data File
        rFile char;  % Results File
        
        % Data
        videoReader VideoReader;    % Video Reader
        calImage uint8;       % Calibration Image
        logData;
        metaData MetaData;       % MetaData
        analysis Analysis;   % Analysis results
        results Results;    % Results
        
        videoProperties;
        
        
        % Booleans
        hasVideo logical;
        hasLog logical;
        hasImage logical;       
        hasMeta logical;
        hasAnalysis logical;
        hasResults logical;
    end
    
    
    
    methods
        function obj = VideoProcessor()
            %UNTITLED6 Construct an instance of this class
            %   Detailed explanation goes here
            obj.hasImage = 0;
            obj.hasMeta = 0;
            obj.hasAnalysis = 0;
            obj.hasResults = 0;
        end
        
        
        function updateVideoFile(obj,vFile,vPath)
            % updateVideoFile(obj,vFile,vPath) Set new video file and change all related file paths
                obj.vFile = vFile;
                obj.path = vPath;
                [~,obj.vName] = fileparts(vFile);
                
                obj.generateDefaultFileNames();
                
                % Create a Vido Reader
                obj.createVideoReader();
                obj.getVideoProperties();
                
                obj.loadLogFile();
                
                
                % If image exists, load it, otherwise, generate it
                if(obj.existImageFile)
                    obj.loadCalibrationImage();
                else
                    obj.generateCalibrationImage();
                end 
                
                % Load metadata if existing or start generating
                if(obj.existMetaFile)
                    obj.loadMetaData();
                else
                    obj.metaData = MetaData(obj);
                    obj.hasMeta = 1;
                end
                
                % Load results if existing
                if(obj.existResultsFile)
                    obj.loadResults();
                else
                    obj.results = Results();
                    obj.hasResults = 1;
                end
                
                
                % Load analysis if existing
                if(obj.existAnalysisFile)
                    obj.loadAnalysis();
                    obj.reassociateAnalysisObjects();
                else
                    obj.analysis = Analysis(obj);
                    obj.hasAnalysis = 1;
                end
                
                
                    
        end
    
        %% Create functions
        % TO BE CHANGED << RESULTS IS A PROPERTY OF ANALYSIS NOT OF
        % VIDEOPROCESSOR
        % reassociateAnalysisObjects(obj)
        % If we load the analysis, the handles within the Analysis class
        % and the VideoProcessor are separate handles to separate classes.
        % The videoReader of the processor class should be used (newly
        % created every time this class is loaded, as files can move. The
        % metaData classes should not have been changed since they have 
        % been saved and thus equal eachother. Similarly, results and
        % analysis have been saved at the same time, so again, they should
        % equal each other. For handle classes: == returns true only if
        % they are the same object, isequal returns true also if objects
        % are duplicates.
        function reassociateAnalysisObjects(obj)
            if(~isequal(obj.metaData,obj.analysis.metaData))
                error('VideoProcessor:reassociateObjects',...
                    'MetaData objects are not equal');
            end
            if(~isequal(obj.results,obj.analysis.results))
                error('VideoProcessor:reassociateObjects',...
                    'Results objects are not equal');
            end
            
            obj.analysis.metaData = obj.metaData;
            obj.analysis.results = obj.results;
            obj.analysis.videoReader = obj.videoReader;
                
        end
        
        function generateCalibrationImage(obj)
            obj.iFile = [obj.vName '.png'];
            v = obj.videoReader;
            v.CurrentTime = 0;
            obj.calImage = v.readFrame();
            v.CurrentTime = 0;
            obj.hasImage = 1;
            obj.saveCalibrationImage();
        end
            
        function getVideoProperties(obj)
            vR = obj.videoReader;
            obj.videoProperties.Duration = vR.Duration;
            obj.videoProperties.Width = vR.Width;
            obj.videoProperties.Height = vR.Height;
            obj.videoProperties.FrameRate = vR.FrameRate;
            obj.videoProperties.nCamsHor = floor(vR.Width/640);
            obj.videoProperties.nCamsVer = floor(vR.Height/480);
            obj.videoProperties.nCams = obj.videoProperties.nCamsHor*obj.videoProperties.nCamsVer;
            
            vR.CurrentTime = 0;
            vR.readFrame();
            
            obj.videoProperties.startFrame = round(vR.currentTime*vR.FrameRate);
            obj.videoProperties.nFrames = round(vR.Duration*vR.FrameRate)...
                -obj.videoProperties.startFrame+1;
        end
        
        
        
        
        %% Image and frame functions
        function annIm = getAnnotatedCamCalIm(obj,camId)
            camIm = obj.getCamCalIm(camId);
            annIm = annotateFromFindImage(camIm,obj.metaData.cam(camId));
            
        end
        
        function camIm = getCamCalIm(obj,camId)
            camIm = getCamImage(obj.calImage,obj.metaData,camId);
        end
        
        function annIm = getAnnotatedCamIm(obj,time,camId)
            camIm = obj.getCamIm(time,camId);
            annIm = annotateFromFindImage(camIm,obj.metaData.cam(camId));
            
        end
        
        function camIm = getCamIm(obj,time,camId)
            image  = obj.readFrame(time);
            camIm = getCamImage(image,obj.metaData,camId);
        end
        
        function image = readFrame(obj,time)
            obj.videoReader.CurrentTime = time;
            if(~obj.videoReader.hasFrame())
                error('VideoProcessor:readFrame','No Frame');
            end
            image = obj.videoReader.readFrame();
        end
        
        function time = getCurrentTime(obj)
            time = obj.videoReader.CurrentTime;
        end
        
        function image = readNextFrame(obj)
            if(~obj.videoReader.hasFrame())
                error('VideoProcessor:readFrame','No Frame');
            end
            image = obj.videoReader.readFrame();
        end
        
        %% Save data file functions
        function saveCalibrationImage(obj)
            % saveCalibrationImage(obj) Saves calibration image
            if(~obj.hasImage)
                error('VideoProcessor:saveNoData','No Image');
            end
            imwrite(obj.calImage,obj.getImageFile());
        end
        
        function saveMetaData(obj)
            % saveMetaData(obj) Saves MetaData
            if(~obj.hasMeta)
                error('VideoProcessor:saveNoData','No MetaData');
            end
            if(~obj.metaData.isFinished())
                error('VideoProcessor:notReady','MetaData not yet completed');
            end
            metaData = obj.metaData; %#ok<NASGU,PROP>
            save(obj.getMetaFile(),'metaData');
        end
        
        function saveAnalysis(obj)
            % saveAnalysis(obj) Saves Analysis
            if(~obj.hasAnalysis)
                error('VideoProcessor:saveNoData','No Analysis');
            end
            analysis = obj.analysis; %#ok<NASGU,PROP>
            save(obj.getAnalysisFile(),'analysis');
        end
        
        function saveResults(obj)
            % saveResults(obj) Saves Results
            if(~obj.hasResults)
                error('VideoProcessor:saveNoData','No Results');
            end
            results = obj.results; %#ok<NASGU,PROP>
            save(obj.getResultsFile(),'results');
        end
        
        %% Load data file functions
        
        function createVideoReader(obj)
            % createVideoReader(obj) Loads the video
            if(~obj.existVideoFile)
                error('VideoProcessor:fileNotFound','Video file %s not found',obj.getVideoFile());
            end
            obj.videoReader = VideoReader(obj.getVideoFile());
            obj.hasVideo = 1;
        end
        
        function loadLogFile(obj)
            % loadLogFile(obj) Loads the video log
            if(~obj.existLogFile)
                error('VideoProcessor:fileNotFound','Log file %s not found',obj.getLogFile());
            end
            
            % Open the text file.
            fileID = fopen(obj.getLogFile(),'r');

            % Read columns of data according to the format.
            % Format for each line of text:
            %   column1: text (%s)
            %	column2: double (%f)
            dataArray = textscan(fileID, '%s%f%[^\n\r]', 'Delimiter',...
                ' ', 'MultipleDelimsAsOne', true, 'TextType', 'string',...
                'EmptyValue', NaN, 'HeaderLines' ,2-1,...
                'ReturnOnError', false, 'EndOfLine', '\r\n');

            % Close the text file.
            fclose(fileID);

            % Create output variable
            text = dataArray{1};
            data = dataArray{2};
            idxN = 1:2:length(data);
            idxT = 2:2:length(data);
            if(~(all(text(idxN)=="n:")&&all(text(idxT)=="t:")))
                error('VideoProcessor:wrongData','Log file %s corrupt',obj.getLogFile());
            end
            lD = struct();
            lD.startFrame = data(idxN(1));
            lD.endFrame  = data(idxN(end));
            lD.frameVector = lD.startFrame :lD.endFrame;
            lD.numberOfFrames = length(lD.frameVector);
            if(~(lD.numberOfFrames==length(idxN)))
                error('VideoProcessor:wrongData','Log file %s corrupt',obj.getLogFile());
            end
            if(~all(data(idxN)==lD.frameVector.'))
                error('VideoProcessor:wrongData','Log file %s corrupt',obj.getLogFile());
            end
            
            lD.startTime = data(idxT(1));
            lD.endTime = data(idxT(end));
            lD.timeVector = data(idxT);
            lD.totalTime = lD.endTime -lD.startTime;
            
            lD.meanFPS = lD.numberOfFrames/(lD.totalTime/1000);
            lD.frameTimes = (lD.timeVector(2:end)-...
                lD.timeVector(1:end-1));
            lD.devFrameTimes = (lD.timeVector(2:end)-...
                lD.timeVector(1:end-1))-1000/lD.meanFPS;
            obj.hasLog = 1;
            obj.logData = lD;
        end
        
        function loadCalibrationImage(obj)
            % loadCalibrationImage(obj) Loads calibration image
            if(~obj.existImageFile)
                error('VideoProcessor:fileNotFound','Image file %s not found',obj.getImageFile());
            end
            obj.calImage = imread(obj.getImageFile());
            obj.hasImage = 1;
        end
        
        function loadMetaData(obj)
            % loadMetaData(obj) Loads MetaData
            if(~obj.existMetaFile)
                error('VideoProcessor:fileNotFound','MetaData file %s not found',obj.getMetaFile());
            end
            m = load(obj.getMetaFile());
            obj.metaData = m.metaData;
            obj.hasMeta = 1;
        end
        
        function loadAnalysis(obj)
            % loadCalibrationImage(obj) Loads calibration image
            if(~obj.existAnalysisFile)
                error('VideoProcessor:fileNotFound','Analysis file %s not found',obj.getAnalysisFile());
            end
            a = load(obj.getAnalysisFile());
            obj.analysis = a.analysis;
            obj.hasAnalysis = 1;
        end
        
        function loadResults(obj)
            % loadMetaData(obj) Loads MetaData
            if(~obj.existResultsFile)
                error('VideoProcessor:fileNotFound','MetaData file %s not found',obj.getMetaFile());
            end
            r = load(obj.getResultsFile());
            obj.results = r.results;
            obj.hasResults = 1;
        end
        
        %% Get file functions
        function fileName = getVideoFile(obj)
            % getVidoFile(obj) Returns the full video file path;
            fileName = fullfile(obj.path,obj.vFile);
        end
        
        function fileName = getLogFile(obj)
            % getLogFile(obj) Returns the full log file path;
            fileName = fullfile(obj.path,obj.lFile);
        end
        
        function fileName = getImageFile(obj)
            % getVidoFile(obj) Returns the full image file path;
            fileName = fullfile(obj.path,obj.iFile);
        end
        
        function fileName = getMetaFile(obj)
            % getVidoFile(obj) Returns the full metadata file path;
            fileName = fullfile(obj.path,obj.mFile);
        end
        
        function fileName = getAnalysisFile(obj)
            % getVidoFile(obj) Returns the full analysis data file path;
            fileName = fullfile(obj.path,obj.aFile);
        end
        
        function fileName = getResultsFile(obj)
            % getVidoFile(obj) Returns the full results file path;
            fileName = fullfile(obj.path,obj.rFile);
        end
        
        %% Exist file functions
        function fileExists = existVideoFile(obj)
            % existVideoFile(obj) Returns 1 if video file exists;
            fileExists = (exist(getVideoFile(obj),'file')==2) ;
        end
        
        function fileExists = existLogFile(obj)
            % existLogFile(obj) Returns 1 if log file exists;
            fileExists = (exist(getLogFile(obj),'file')==2) ;
        end
        
        function fileExists = existImageFile(obj)
            % existImageFile(obj) Returns 1 if image file exists;
            fileExists = (exist(getImageFile(obj),'file')==2) ;
        end
        
        function fileExists = existMetaFile(obj)
            % existMetaFile(obj) Returns 1 if video file exists;
            fileExists = (exist(getMetaFile(obj),'file')==2) ;
        end
        
        function fileExists = existAnalysisFile(obj)
            % existAnalysisFile(obj) Returns 1 if analysis file exists;
            fileExists = (exist(getAnalysisFile(obj),'file')==2) ;
        end
        
        function fileExists = existResultsFile(obj)
            % existResultsFile(obj) Returns 1 if results file exists;
            fileExists = (exist(getResultsFile(obj),'file')==2) ;
        end
        
        %%
        % Generate default file names
        function generateDefaultFileNames(obj)
            obj.iFile = [obj.vName '.jpg'];
            obj.mFile = [obj.vName '_meta.mat'];
            obj.lFile = [obj.vName '.log'];
            obj.aFile = [obj.vName '_analysis.mat'];
            obj.rFile = [obj.vName '_results.mat'];
        end
        
    end
end


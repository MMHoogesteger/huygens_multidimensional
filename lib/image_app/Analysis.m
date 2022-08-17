classdef Analysis < handle
    %Analysis - Class governing analysis and processing of videos.
    % This class contains all logic and function necessary to analyze a
    % certain video. Meant to be a member of a VideoProcessor object and
    % uses a Results object to store analysis results.
    %
    % Analysis is used of several steps, which are in order:
    % Video processing
    % Platform movements
    % Clustering
    % Deprojection
    % Signal analysis
    %
    % See also: VideoProcessor, MetaData, Results
    
    % Author: M.M. Hoogesteger - Msc student of Mechanical Engineering
    % Eindhoven University of Technology, Mechanical Engineering, Dynamics and
    % control
    % email address: m.m.hoogesteger@student.tue.nl
    % February 2018; Last revision: 8-October-2018
    
    properties
        videoProperties             % Struct containing video properties
        videoReader VideoReader;    % VideoReader object with current video
        metaData MetaData;          % MetaData object with world information
        results Results;            % Results object to store data
        
        nFrames;                    % Number of frames
        curFrame;                   % Current frame
        preparationFinished logical;% Boolean: preparation is finished
        params;                     % Containig some parameters for analysis
        
    end
    
    methods
        %% Constructor
        
        function obj = Analysis(videoProcessor)
            %Analysis(videoProcessor) Construct an instance of this class
            %   The videoProcessor is an instance of the class
            %   VideoProcessor and needs to be fully initialized.
            obj.videoProperties = videoProcessor.videoProperties;
            obj.videoReader = videoProcessor.videoReader;
            obj.metaData = videoProcessor.metaData;
            obj.results = videoProcessor.results;
            
            obj.nFrames = obj.videoProperties.nFrames;
            obj.preparationFinished = 0;
            
            obj.params.maxRotation = 0.25*pi; % 45 degrees
            obj.params.maxTranslation = 100; % 10cm
        end
        
        
        %% Video processing functions
        
        function prepareAnalysis(obj)
            %prepareAnalysis(obj) Prepare the analysis
            % Prepares the results object with empty arrays of the correct
            % size and sets all counters to the correct value.
            v = obj.videoReader;
            m = obj.metaData;
            
            % Generate some cell arrays
            obj.results.mPos = cell(obj.nFrames,m.nCams);
            obj.results.pPos = cell(obj.nFrames,m.nCams);
            obj.results.cPos = cell(obj.nFrames,m.nCams);
            obj.results.cOr  = cell(obj.nFrames,m.nCams);
            obj.results.ns   = cell(obj.nFrames,m.nCams);
            
            
            % Set time at zero
            obj.curFrame = 0;
            v.CurrentTime = 0;
            obj.preparationFinished = 1;
        end
        
        function stepAnalysis(obj)
            %stepAnalysis(obj) Take one analysis step
            % After preparing the analysis takes one step to walk through
            % the video and processes one frame.
            
            if(~obj.preparationFinished)
                error('Analysis:PreparationNotFinished',...
                    'Analysis.prepareAnalysis has not yet been run.');
            end
            v = obj.videoReader;
            m = obj.metaData;
            
            % Increase frame number and read a frame
            cF = obj.curFrame+1;
            obj.curFrame = cF;
            im = readFrame(v);
            
            
            % Loop over all cameras
            for camIdx = 1:m.nCams
                
                % Find all features
                camIm = getCamImage(im,m,camIdx);
                checkers = findCheckers(camIm);
                mMarkers = findMarkers(camIm,m.cam(camIdx).colLims,m.cam(camIdx).blobLims);
                pMarkers = findMarkers(camIm,m.cam(camIdx).pcolLims,m.cam(camIdx).pblobLims);
                
                % Store numbers
                nC = checkers(1).n;
                nM = mMarkers(1).n;
                nP = pMarkers(1).n;
                
                
                % Allocate temp variables
                tCOr = zeros(nC,1);
                tCPos = zeros(nC,2);
                tMPos = zeros(nM,2);
                tPPos = zeros(nP,2);
                
                % For every feature, determine point coordinates in the
                % world frame and for checkerboars orientation of the
                % principal axis.
                for cid = 1:nC
                    cwP = transformCamToWorld(m,checkers(cid).imagePoints,camIdx);
                    tCPos(cid,:) = cwP(1,:);
                    yVec = cwP(2,:)-cwP(1,:);
                    tCOr(cid) = atan2(yVec(2),yVec(1));
                end
                
                for mid = 1:nM
                    tMPos(mid,:) = transformCamToWorld(m,mMarkers(mid).center,camIdx);
                end
                
                for pid = 1:nP
                    tPPos(pid,:) = transformCamToWorld(m,pMarkers(pid).center,camIdx);
                end
                
                % Store data
                obj.results.ns{cF,camIdx}   = [nC nM nP];
                obj.results.cOr{cF,camIdx}  = tCOr;
                obj.results.cPos{cF,camIdx} = tCPos;
                obj.results.mPos{cF,camIdx} = tMPos;
                obj.results.pPos{cF,camIdx} = tPPos;
                
                
            end
            
            
        end
        
        function runAnalysis(obj,progressFunctionCallBack)
            %runAnalysis(obj,progressFunctionCallBack) Run analysis
            % Runs analysis for the rest of the video until it has no more
            % frames.
            % progressFunctionCallBack will be called as follows after
            % every ten frames :
            % progressFunctionCallBack(curFrame,nFrames,tRemMin), where
            % curFrame is the current Frame, nFrames is the total number of
            % frames and tRemMin is the expected remaining time and
            % minutes.
            %
            % See also: stepAnalysis(obj)
            if(~obj.preparationFinished)
                error('Analysis:PreparationNotFinished',...
                    'Analysis.prepareAnalysis has not yet been run.');
            end
            v = obj.videoReader;
            t1 = tic;
            
            % Set counter for callback to zero
            pCounter = 0;
            while(v.hasFrame())
                
                % Analyze a frame.
                obj.stepAnalysis();
                pCounter = pCounter + 1;
                
                % Call progressCallback if tenth frame.
                if(pCounter == 10)
                    t2 = toc(t1);
                    tTotEst = t2*obj.nFrames/obj.curFrame;
                    tRemMin = (tTotEst-t2)/60;
                    progressFunctionCallBack(obj.curFrame,obj.nFrames,round(tRemMin));
                    pCounter = 0;
                end
            end
            % Finalize analysis.
            obj.finishAnalysis();
        end
        
        function finishAnalysis(obj)
            %finishAnalysis(obj) Set anlysis to finished
            % Set pointsGathered to true
            obj.results.pointsGathered = true;
        end
        
        
        %% Functions used in after analysis
        
        function solvePlatformMovements(obj)
            %solvePlatformMovements(obj) Infers platform movements from data.
            % Takes raw data and fits a rotation and translation to infer
            % platform movement. Tolerances as defined in the params
            % property are used here.
            
            nF = obj.nFrames;
            r = obj.results;
            
            % Initialize data storage and initial guesses.
            r.phiStor = zeros(nF,1);
            r.tStor = zeros(nF,2);
            r.RStor = zeros(nF,2,2);
            
            phi0 = 0;
            t0 = [0 0];
            for n = 1:nF
                
                % Select platform and checkoerboard features from frame n
                cOrT  = obj.results.cOr(n,:);
                cPosT = obj.results.cPos(n,:);
                pPosT = obj.results.pPos(n,:);
                
                % Total number of features in the frame (sum over all
                % cameras)
                nT = sum(cellfun(@(c)sum(c([1,3])),...
                    obj.results.ns(n,:),'Uni',1));
                nC = sum(cellfun(@(c)c(1),...
                    obj.results.ns(n,:),'Uni',1));
                
                % Find optimal angle and translation based on checkerboard
                % and platform markers. Need at least two features or
                % one checkerboard.
                if(nT>1 || nC>0)
                    [phi,t,R] = findPlatFormOrientation(cOrT,cPosT,pPosT,obj.metaData,phi0,t0);
                    
                    % Don't allow wrong results: check max rotation and
                    % translation.
                    if(abs(phi)>obj.params.maxRotation || norm(t,2)>obj.params.maxTranslation)
                        phi = NaN;
                        t = NaN;
                        R = nan(4,4);
                    else
                        phi0 = phi;
                        t0 = t;
                    end
                else
                    phi = NaN;
                    t = NaN;
                    R = nan(4,4);
                end
                
                % Store inferred rotation, translation and the
                % rotationMatrix
                r.phiStor(n,1) = phi;
                r.tStor(n,1:2) = t;
                r.RStor(n,:,:) = R;
                
                
            end
            r.solvedMovement = true;
            
        end
        
        function correctForPlatformMovements(obj)
            %correctForPlatformMovements(obj) correct metronome and camera
            % Takes metronome and camera positions and put them in the
            % platform coordinate frame, ready for deprojection.
            m = obj.metaData;
            r = obj.results;
            
            
            r.cCCor   = cell(obj.nFrames,m.nCams);
            r.mPosCor = cell(obj.nFrames,m.nCams);
            
            for n = 1:obj.nFrames
                
                R = squeeze(r.RStor(n,:,:));
                t = r.tStor(n,:);
                
                for camIdx = 1:m.nCams
                    r.cCCor{n,camIdx}  = (m.world.camCenters{camIdx}-t)*R.';
                    r.mPosCor{n,camIdx} = (obj.results.mPos{n,camIdx}-t)*R.';
                    
                end
                
            end
            r.pointsCorrected = true;
            
            
        end
        
        
        %% Image clustering
        
        function initClustering(obj)
            %initClustering(obj) Initialize clustering
            % Create lists and gather points etc.
            m = obj.metaData;
            r = obj.results;
            r.clusters = struct;
            r.clusters.autoClustered = false;
            r.metronomes = struct([]);
            for camId = 1:m.nCams
                % Store all points in a list
                r.clusters.cam(camId).mPosList = cell2mat(r.mPosCor(:,camId));
                r.clusters.cam(camId).mPosXList = r.clusters.cam(camId).mPosList(:,1);
                r.clusters.cam(camId).mPosYList = r.clusters.cam(camId).mPosList(:,2);
                d = num2cell((1:obj.nFrames).');
                r.clusters.cam(camId).nFrameList = cell2mat(cellfun(...
                    @(c,d) ones(size(c,1),1).*d,r.mPosCor(:,camId),d,'Uni',0));
                r.clusters.cam(camId).nList = cell2mat(cellfun(...
                    @(c,d) (1:size(c,1)).',r.mPosCor(:,camId),d,'Uni',0));
                
                r.clusters.cam(camId).hasPolCluster = false;
            end
        end
        
        function clusterPoints(obj)
            %clusterPoints(obj) Deprecated - cluster points with grid
            m = obj.metaData;
            r = obj.results;
            
            r.grid.resolution = 0.25; % mm
            r.grid.overlap = 10;
            r.grid.se = strel('disk',10);
            for camId = 1:m.nCams
                g = struct;
                g.resolution = r.grid.resolution;
                g.overlap = r.grid.overlap;
                
                
                
                % Gather valid results
                vPoints = ~isnan(g.mPosXList);
                
                g.mPosXList = g.mPosXList(vPoints);
                g.mPosYList = g.mPosYList(vPoints);
                g.nFrameList = g.nFrameList(vPoints);
                g.nList = g.nList(vPoints);
                
                % Grid boundaries
                g.minX = min(g.mPosXList);
                g.minY = min(g.mPosYList);
                g.maxX = max(g.mPosXList);
                g.maxY = max(g.mPosYList);
                
                g.gridX = g.minX:g.resolution:g.maxX+g.resolution;
                g.gridY = g.minY:g.resolution:g.maxY+g.resolution;
                
                
                g.gridIm = zeros(length(g.gridY)-1-g.overlap,...
                    length(g.gridX)-1-g.overlap);
                g.gridPointList = cell(length(g.gridY)-1-g.overlap,...
                    length(g.gridX)-1-g.overlap);
                for idX = 1:length(g.gridX)-1-g.overlap
                    xMin = g.gridX(idX);
                    xMax = g.gridX(idX+1+g.overlap);
                    for idY = 1:length(g.gridY)-1-g.overlap
                        yMin = g.gridY(idY);
                        yMax = g.gridY(idY+1+g.overlap);
                        
                        idxs = find(g.mPosXList>=xMin & g.mPosXList<xMax...
                            & g.mPosYList>=yMin & g.mPosYList<yMax);
                        g.gridIm(idY,idX) = length(idxs);
                        g.gridPointList{idY,idX} = idxs;
                    end
                end
                
                gridIm = g.gridIm~=0;
                figure;
                imshow(gridIm)
                
                
                
                closeBW = imclose(gridIm,r.grid.se);
                figure, imshow(closeBW)
                D = bwdist(~closeBW);
                figure, imshow(D,[])
                D = -D;
                D(~closeBW) = Inf;
                figure, imshow(D,[])
                D = imhmin(D,20);
                figure, imshow(D,[])
                g.L = watershed(D);
                g.L(~closeBW) = 0;
                rgb = label2rgb(g.L,'jet',[.5 .5 .5]);
                figure
                imshow(rgb,'InitialMagnification','fit')
                title('Watershed transform of D')
                
                r.grid.cam(camId) = g;
            end
        end
        
        function autoClusterPointsGrid(obj,maxClusts)
            %autoClusterPointsGrid(obj,maxClusts) Euclidian clustering
            % maxClusts is the number of clusters to create
            m = obj.metaData;
            r = obj.results;
            r.clusters.autoClusters = struct;
            for camId = 1:m.nCams
                % Store all points in a list
                r.clusters.autoClusters(camId).mPosList = r.clusters.cam(camId).mPosList;
                r.clusters.autoClusters(camId).linkage = linkage(r.clusters.autoClusters(camId).mPosList,'single');
                r.clusters.autoClusters(camId).cluster = cluster(r.clusters.autoClusters(camId).linkage,'Maxclust',maxClusts(camId));
                for clId = 1:max(r.clusters.autoClusters(camId).cluster)
                    r.clusters.autoClusters(camId).clList{clId} = find(r.clusters.autoClusters(camId).cluster==clId);
                end
                
            end
            r.clusters.autoClustered = true;
            
        end
        
        function createNewMetronome(obj)
            %createNewMetronome(obj) Add new metronome (point group)
            r = obj.results;
            nCur = numel(r.metronomes);
            met = struct;
            met.posList = [];
            met.posXList = [];
            met.posYList = [];
            met.nFrameList = [];
            met.nList = [];
            met.camList = [];
            
            met.nPoints = 0;
            met.nFramesSingle = 0;
            met.nFramesMulti = 0;
            met.nFramesNot = obj.nFrames;
            
            if(nCur == 0)
                r.metronomes = met;
            else
                r.metronomes(nCur+1) = met;
            end
            
            
            
        end
        
        function addClusterToMetronome(obj,metId,clusterIds,camId)
            %addClusterToMetronome(obj,metId,clusterIds,camId)
            % metId is the metronome ID
            % clusterIds are all point IDs (not cluster number)
            % camId is the camera ID
            r = obj.results;
            
            met = r.metronomes(metId);
            nP = numel(r.clusters.cam(camId).mPosList(clusterIds));
            
            met.posList = [met.posList;r.clusters.cam(camId).mPosList(clusterIds,:)];
            met.posXList = [met.posXList;r.clusters.cam(camId).mPosXList(clusterIds)];
            met.posYList = [met.posYList;r.clusters.cam(camId).mPosYList(clusterIds)];
            met.nFrameList = [met.nFrameList;r.clusters.cam(camId).nFrameList(clusterIds)];
            met.nList = [met.nList;r.clusters.cam(camId).nList(clusterIds)];
            met.camList = [met.camList;camId*ones(nP,1)];
            
            met.nPoints = size(met.posList,1);
            met.nFramesSingle = getFramesSingle(obj,met);
            met.nFramesMulti = getFramesMulti(obj,met);
            met.nFramesNot = getFramesNone(obj,met);
            
            r.metronomes(metId) = met;
        end
        
        function removeMetronome(obj,id)
            %removeMetronome(obj,id) Remove metronome with id = id.
            r = obj.results;
            nCur = numel(r.metronomes);
            tempMet = r.metronomes;
            r.metronomes = struct;
            r.metronomes = tempMet([1:id-1,id+1:nCur]);
        end
        
        function fN = getFramesNone(obj,met)
            %getFramesNone(obj,met) Get number of frames without detection for metronome met
            % met is a struct, not an id
            fN = obj.nFrames-numel(unique(met.nFrameList));
        end
        
        function fN = getFramesSingle(obj,met)
            %getFramesNone(obj,met) Get number of frames with one detection for metronome met
            % met is a struct, not an id
            fN = numel(find(histcounts(met.nFrameList,'BinMethod','integers')==1));
        end
        
        function fN = getFramesMulti(obj,met)
            %getFramesNone(obj,met) Get number of frames with multiple detection for metronome met
            % met is a struct, not an id
            fN = numel(find(histcounts(met.nFrameList,'BinMethod','integers')>1));
        end
        
        function finishClustering(obj)
            %finishClustering(obj) Finalize clustering
            obj.results.clustered = true;
            obj.results.deprojected = false;
            obj.results.resolved = false;
        end
        
        
        %% Deprojection
        
        function deprojectMetronome(obj,metId)
            %deprojectMetronome(obj,metId) Deproject metronome with id = id
            % Calculates the optimal placement and swing direction of a
            % metronome so that the line from camera position through
            % metronome tip optimally coincides with the measured points on
            % the platform surface.
            %
            % See also: deprojectAllMetronomes
            m = obj.metaData;
            r = obj.results;
            
            % Gather all points and camera positions
            met = r.metronomes(metId);
            
            camCenters = zeros(met.nPoints,2);
            camZs = zeros(met.nPoints,1);
            for camId = 1:m.nCams
                pointsIds = find(met.camList ==camId);
                frameList = met.nFrameList(pointsIds);
                camCenters(pointsIds,:) = cell2mat(r.cCCor(frameList,camId));
                camZs(pointsIds) = -m.world.camZs{camId};
            end
            r.metronomes(metId).camPositions = camCenters;
            
                       
            % Get an approximate direction by calculating the main
            % direction based on a svd.
            meanPoint0 = mean([met.posList]);
            [~,~,V] = svd(met.posList-meanPoint0,'econ');
            Vmain = V(:,1);
            psi0 = atan2(Vmain(2),Vmain(1));
            
            % Optimally reconstruct metronome position
            camEst = mean(camCenters);
            H = 36;
            L = 70;
            camZEst = mean(camZs);
            Mest = (camEst - meanPoint0) * (H+L)/camZEst + meanPoint0;
            X0 = [psi0*10;Mest.'];
            
            options =optimset('MaxFunEvals',1e5,...
                'MaxIter',1e5,...
                'Display','final',...
                'TolFun',1e-10,...
                'TolX',1e-10);
            X = fminsearch(@(x) calcReProjectionErrorNoCam(met.posList,camCenters,x(1)/10,H,L,camZs,x(2:3)),X0,options);
            [E,Mx,My,Mrpx,Mrpy,Mtheta] = calcReProjectionErrorNoCam(met.posList,camCenters,X(1)/10,H,L,camZs,X(2:3));
            r.metronomes(metId).Mx = Mx;
            r.metronomes(metId).My = My;
            r.metronomes(metId).Mrpx = Mrpx;
            r.metronomes(metId).Mrpy = Mrpy;
            r.metronomes(metId).Mtheta = Mtheta;
            
            r.metronomes(metId).psi = X(1)/10;
            r.metronomes(metId).centerPos = X(2:3);
            
        end
        
        function deprojectAllMetronomes(obj)
            %deprojectAllMetronomes(obj) Deproject all metronomed
            %
            % See also: deprojectMetronome
            r = obj.results;
            nMets = numel(r.metronomes);
            for metId = 1:nMets
                obj.deprojectMetronome(metId);
            end
            r.deprojected = true;
        end
        
        function flipMetronome(obj,metId)
            %flipMetronome(obj,metId) Flips a metronome sign
            % metId is the id of the metronome to be flipped
            r = obj.results;
            r.metronomes(metId).psi = wrapToPi(r.metronomes(metId).psi + pi);
            r.metronomes(metId).Mtheta = -r.metronomes(metId).Mtheta;
        end
        
        %% Signals
        function resolveSignals(obj,timeVector)
            %resolveSignals(obj,timeVector)
            % Transfer optimally reconstructed metronome positions to
            % signals with corresponding times as timeVector
            if(numel(timeVector)>obj.nFrames)
                warning('Analysis:resolveSignalsTime',...
                    'Time vector to long, using last part');
                timeVector = timeVector(end-obj.nFrames+1:end);
            elseif(numel(timeVector)<obj.nFrames)
                error('Analysis:resolveSignalsTime',...
                    'Time vector to short, aborting');
            end
            r = obj.results;
            r.signals = struct;
            r.signals.x_f = r.tStor(:,1);
            r.signals.y_f = r.tStor(:,2);
            r.signals.phi_f = r.phiStor;
            r.signals.t = 1e-3*(timeVector-timeVector(1));
            
            
            nMets = numel(r.metronomes);
            % Construct signals. Single detections have preference. Then
            % chosen multi detections. Finally interpolate.
            for metId = 1:nMets
                met = r.metronomes(metId);
                multiFrameIds = find(histcounts(met.nFrameList,'BinMethod','integers')>1);
                singleFrameIds = find(histcounts(met.nFrameList,'BinMethod','integers')==1);
                
                singleListIds = ismember(met.nFrameList,singleFrameIds);
                multiListIds = ismember(met.nFrameList,multiFrameIds);
                
                frameVecSingle = met.nFrameList(singleListIds);
                frameVecMulti = unique(met.nFrameList(multiListIds));
                thetaVecSingle = met.Mtheta(singleListIds);
                
                thetaVecSingleSmoothed = smooth(frameVecSingle,thetaVecSingle,10/obj.nFrames,'rloess');
                thetaVecMultiEst = interp1(frameVecSingle,thetaVecSingleSmoothed,frameVecMulti,'spline');
                thetaVecMulti = nan(1,numel(frameVecMulti));
                for n = 1:numel(frameVecMulti)
                    thetaMultis = met.Mtheta(met.nFrameList==frameVecMulti(n));
                    dists = abs(thetaMultis-thetaVecMultiEst(n));
                    [~,minDistId] = min(dists);
                    thetaVecMulti(n) = thetaMultis(minDistId);
                end
                thetaVecFinal = nan(obj.nFrames,1);
                thetaVecFinal(frameVecSingle) = thetaVecSingle;
                thetaVecFinal(frameVecMulti) = thetaVecMulti;
                
                frameNaNIds = find(isnan(thetaVecFinal));
                thetaVecFinal(frameNaNIds) = interp1(1:obj.nFrames,thetaVecFinal,frameNaNIds,'spline');
                
                r.signals.met(metId).theta_f = thetaVecFinal;
                
                
            end
            r.resolved = true;
            
        end
        
        function analyzeSignals(obj,Tss)
            %analyzeSignals(obj,Tss)
            % Analyze signals using the same tools as simulations: hilbert
            % transforms, fourier series, etc.
            r = obj.results;
            r.analysis = struct;
            nMets = numel(r.metronomes);
            
            f_ids = r.signals.t>Tss;
            r.signals.x_mean = mean(r.signals.x_f(f_ids));
            r.signals.x_cor = r.signals.x_f - r.signals.x_mean;
            r.signals.xhat_cor = hilbert(r.signals.x_cor);
            
            r.signals.y_mean = mean(r.signals.y_f(f_ids));
            r.signals.y_cor = r.signals.y_f - r.signals.y_mean;
            r.signals.yhat_cor = hilbert(r.signals.y_cor);
            
            r.signals.phi_mean = mean(r.signals.phi_f(f_ids));
            r.signals.phi_cor = r.signals.phi_f - r.signals.phi_mean;
            r.signals.phihat_cor = hilbert(r.signals.phi_cor);
            for metId = 1:nMets
                r.signals.met(metId).theta_mean = mean(r.signals.met(metId).theta_f(f_ids));
                r.signals.met(metId).theta_cor = r.signals.met(metId).theta_f - r.signals.met(metId).theta_mean;
                r.signals.met(metId).thetahat_cor = hilbert(r.signals.met(metId).theta_cor);
            end
            
            w = 10/1.728;
            t = r.signals.t;
            x_cor = r.signals.x_cor;
            xhat_cor = r.signals.xhat_cor;
            
            y_cor = r.signals.y_cor;
            yhat_cor = r.signals.yhat_cor;
            
            phi_cor = r.signals.phi_cor;
            phihat_cor = r.signals.phihat_cor;
            
            met = r.signals.met;
            
            
            
            
            r.analysis.phase_x = nan(numel(t),1);
            r.analysis.phase_y = nan(numel(t),1);
            r.analysis.phase_phi = nan(numel(t),1);
            r.analysis.phase_mets = nan(numel(t),nMets);
            
            r.analysis.dphase_x = nan(numel(t),1);
            r.analysis.dphase_y = nan(numel(t),1);
            r.analysis.dphase_phi = nan(numel(t),1);
            r.analysis.dphase_mets = nan(numel(t),nMets);
            
            r.analysis.amp_x = nan(numel(t),1);
            r.analysis.amp_y = nan(numel(t),1);
            r.analysis.amp_phi = nan(numel(t),1);
            r.analysis.amp_mets = nan(numel(t),nMets);
            
            
            
            for n = 1:numel(t)
                ind = find(t>=t(n) & t<(t(n)+w));
                
                met1_temp = met(1).theta_cor(ind);
                met1_imag_temp = imag(met(1).thetahat_cor(ind));
                
                temp1 = dot(x_cor(ind),met1_temp)/(norm(x_cor(ind))*norm(met1_temp));
                temp2 = dot(x_cor(ind),met1_imag_temp)/(norm(x_cor(ind))*norm(met1_imag_temp));
                r.analysis.dphase_x(n) = atan2(temp2,temp1);
                r.analysis.amp_x(n) = max(abs(xhat_cor(ind)));
                r.analysis.phase_x(n) = angle(xhat_cor(n));
                
                temp1 = dot(y_cor(ind),met1_temp)/(norm(y_cor(ind))*norm(met1_temp));
                temp2 = dot(y_cor(ind),met1_imag_temp)/(norm(y_cor(ind))*norm(met1_imag_temp));
                r.analysis.dphase_y(n) = atan2(temp2,temp1);
                r.analysis.amp_y(n) = max(abs(yhat_cor(ind)));
                r.analysis.phase_y(n) = angle(yhat_cor(n));
                
                temp1 = dot(phi_cor(ind),met1_temp)/(norm(phi_cor(ind))*norm(met1_temp));
                temp2 = dot(phi_cor(ind),met1_imag_temp)/(norm(phi_cor(ind))*norm(met1_imag_temp));
                r.analysis.dphase_phi(n) = atan2(temp2,temp1);
                r.analysis.amp_phi(n) = max(abs(phihat_cor(ind)));
                r.analysis.phase_phi(n) = angle(phihat_cor(n));
                
                for metId = 1:nMets
                    temp1 = dot(met(metId).theta_cor(ind),met1_temp)/(norm(met(metId).theta_cor(ind))*norm(met1_temp));
                    temp2 = dot(met(metId).theta_cor(ind),met1_imag_temp)/(norm(met(metId).theta_cor(ind))*norm(met1_imag_temp));
                    r.analysis.dphase_mets(n,metId) = atan2(temp2,temp1);
                    r.analysis.amp_mets(n,metId) = max(abs(met(metId).thetahat_cor(ind)));
                    r.analysis.phase_mets(n,metId) = angle(met(metId).thetahat_cor(n));
                end
                
                
            end
            
            r.analysis.freq_x = diff(smooth(unwrap(r.analysis.phase_x),8))/((max(t)-min(t))/numel(t))/2/pi;
            r.analysis.freq_y = diff(smooth(unwrap(r.analysis.phase_y),8))/((max(t)-min(t))/numel(t))/2/pi;
            r.analysis.freq_phi = diff(smooth(unwrap(r.analysis.phase_phi),8))/((max(t)-min(t))/numel(t))/2/pi;
            r.analysis.freq_mets = nan(size(r.analysis.phase_mets));
            r.analysis.freq_mets = r.analysis.freq_mets(1:end-1,:);
            for metId = 1:nMets
                r.analysis.freq_mets(:,metId) = diff(smooth(unwrap(r.analysis.phase_mets(:,metId)),100),1,1)/((max(t)-min(t))/numel(t))/2/pi;
            end
            
            tHilbertCutoff = 3;  %s
            idSS = find(t>Tss & t<(t(end)-tHilbertCutoff));
            
            tSS = t(idSS);
            r.analysis.tSS = tSS;
            omega0 = mean(r.analysis.freq_mets(idSS,1))*2*pi;
            theta1 = r.signals.met(metId).theta_f(idSS);
            T = tSS(end)-tSS(1);
            R = 6;
            fc0 = @(omega) fourierDC(tSS,theta1,omega);
            fc = @(omega) fourierIntegral(tSS,theta1,omega,R);
            fJ = @(omega) sum((theta1.'-fourierReconstruct(tSS,fc0(omega),fc(omega),omega)).^2);
            options = optimset('TolFun',1e-10,'TolX',1e-10);
            omega = fminsearch(fJ,omega0,options);
            
            
            Nq = size(3+nMets,2);
            
            
            Fz0 = zeros(Nq,1);
            Fzk = zeros(Nq,R);
            Ezk = zeros(Nq,R);
            
            Fz0(1) = fourierDC(tSS,r.signals.x_f(idSS),omega);
            Fzk(1,:) = fourierIntegral(tSS,r.signals.x_f(idSS),omega,R);
            
            Fz0(2) = fourierDC(tSS,r.signals.y_f(idSS),omega);
            Fzk(2,:)  = fourierIntegral(tSS,r.signals.y_f(idSS),omega,R);
            
            Fz0(3) = fourierDC(tSS,r.signals.phi_f(idSS),omega);
            Fzk(3,:)  = fourierIntegral(tSS,r.signals.phi_f(idSS),omega,R);
            for rm = 1:R
                fr = fourierReconstruct(tSS,Fz0(1),Fzk(1,1:rm),omega);
                Ezk(1,rm) = sqrt(trapz(tSS,(r.signals.x_f(idSS)-fr.').^2)./T);
                
                fr = fourierReconstruct(tSS,Fz0(2),Fzk(2,1:rm),omega);
                Ezk(2,rm) = sqrt(trapz(tSS,(r.signals.y_f(idSS)-fr.').^2)./T);
                
                fr = fourierReconstruct(tSS,Fz0(3),Fzk(3,1:rm),omega);
                Ezk(3,rm) = sqrt(trapz(tSS,(r.signals.phi_f(idSS)-fr.').^2)./T);
            end
            for metId = 1:nMets
                
                Fz0(3+metId) = fourierDC(tSS,r.signals.met(metId).theta_f(idSS),omega);
                Fzk(3+metId,:) = fourierIntegral(tSS,r.signals.met(metId).theta_f(idSS),omega,R);
                
                for rm = 1:R
                    fr = fourierReconstruct(tSS,Fz0(3+metId),Fzk(3+metId,1:rm),omega);
                    Ezk(3+metId,rm) = sqrt(trapz(tSS,(r.signals.met(metId).theta_f(idSS)-fr.').^2)./T);
                end
            end
            
            Fzkn = normalizeFourierCoefficients(Fzk,Fzk(4,1));
            
            fourier.R = R;
            fourier.T = T;
            fourier.omega = omega;
            fourier.Fz0 = Fz0;
            fourier.Fzk = Fzk;
            fourier.Fzkn = Fzkn;
            fourier.Ezk = Ezk;
            
            r.analysis.fourier = fourier;
            r.analysis.idSS = idSS;
            r.analyzed = true;
        end
        
    end
end


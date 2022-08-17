classdef Results < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        % Raw detections
        mPos
        pPos
        cPos
        cOr
        
        % Number of features in every frame
        ns
        
        % 
        phiStor
        tStor
        RStor
        
        pointsGathered logical;
        
        solvedMovement logical;
        
        pointsCorrected logical;
        
        clustered logical;
        
        deprojected logical;
        
        resolved logical;
        
        analyzed logical;
        
        cCCor
        mPosCor
        
        grid
        
        Mx
        My
        Mrpx
        Mrpy
        
        clusters
        metronomes
        
        signals
        analysis
    end
    
    methods
        function obj = Results()
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
        end
        
    end
end


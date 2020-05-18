classdef LocationEstimator_tdoa
    % Superclass for all Location Estimators relying on TDoA measurements
    
    properties
        Xenb % 2xN-matrix containing the source positions 
        % as column 2-vectors. 
        % The reference source is in column 1.
        inParallel = 0
    end
    
    methods
%         function obj = LocationEstimator(Xenb)
%             % Construct an instance of this class
%            
%             obj.Xenb = Xenb;
%         end
        
        [v_estimatedLocation, locUncertainty] = ...
            estimateOneLocation(obj,v_measurements)
        % INPUT:  (N-1)-vector containing the TDOA measurements
        % OUTPUT: v_estimatedLocation: 2-vector, location estimate
        %         locUncertainty:      scalar, uncertainty measure in loc
        
        function [m_estimatedLocations, v_locUncertainties] = ...
                estimateLocations(obj, m_measurements)
            % This is the method of the superclass - can be overridden
            % with code that estimates multiple locations efficiently
            % INPUT:  (N-1)xM matrix containing M sets of TDOA measurements
            % OUTPUT: 2xM matrix     containing M location estimates
            % M can be seen as the number of UEs
            
            M = size(m_measurements, 2);
            m_estimatedLocations = zeros(2, M);
            v_locUncertainties   = zeros(1, M);
            if obj.inParallel
                disp('running in parallel.')
                parfor m = 1:M
                    [m_estimatedLocations(:,m), v_locUncertainties(m)] = ...
                        obj.estimateOneLocation(m_measurements(:, m)); ...
                        %#ok<PFBNS>                    
                end
            else
                ltc = LoopTimeControl(M);
                for m = 1:M
                    [m_estimatedLocations(:,m), v_locUncertainties(m)] = ...
                        obj.estimateOneLocation(m_measurements(:, m)); 
                    ltc.go(m);
                end
            end
        end
            
    end
end


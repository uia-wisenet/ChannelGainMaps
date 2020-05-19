classdef LocationEstimator_tdoa
    % Superclass for all Location Estimators relying on TDoA measurements
    
    properties
        Xenb % 2xN-matrix containing the source positions 
        % as column 2-vectors. 
        % The reference source is in column 1.
        b_inParallel = 0
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
            % OUTPUT: Mx2 matrix     containing M location estimates
            % M can be seen as the number of UEs
            
            n_loc = size(m_measurements, 2);
            m_estimatedLocations = zeros(n_loc, 2);
            v_locUncertainties   = zeros(n_loc, 1);
            if obj.b_inParallel
                disp('running in parallel.')
                parfor m = 1:n_loc
                    [m_estimatedLocations(m,:), v_locUncertainties(m)] = ...
                        obj.estimateOneLocation(m_measurements(:, m)); ...
                        %#ok<PFBNS>                    
                end
            else
                ltc = LoopTimeControl(n_loc);
                for m = 1:n_loc
                    [m_estimatedLocations(m,:), v_locUncertainties(m)] = ...
                        obj.estimateOneLocation(m_measurements(:, m)); 
                    ltc.go(m);
                end
            end
        end
            
    end
end


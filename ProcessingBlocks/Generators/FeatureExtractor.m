classdef FeatureExtractor
    % FeatureExtractor  class.
    %  This class allows to extract the LocF and LocB features from pilot
    %  signals.
    properties (Constant)
        c = 3e8;
    end
    properties
        sampling_period
    end
    
    methods
        function  estimatedCentersOfMass  = locFreeExtract(obj,pilotSignals)
            % this function extracts the LocF features from pilot signals.
            % Accepts input pilotSignals of size S-by-K-by-N where S is the
            % number of sources (or anchor nodes), K is the number of samples
            % per signal, and N is the number of sensor locations.
            % Outputs a matrix of size M-by-N where M is the number of
            % LocF features extracted at one sensor location.
            n_ues=size(pilotSignals,3);
            n_sources=size(pilotSignals,1);
            estimatedCentersOfMass=zeros(nchoosek(n_sources,2),n_ues);
            combin_sources=combnk(1:n_sources,2);
            if n_sources < 6
                combin_sources=flipud(combin_sources);
            end
            for ind_ues=1:n_ues
                allCenterOfMass=zeros(1,nchoosek(n_sources,2));
                for ind_comb=1:size(combin_sources,1)
                    [hD_corr, lags]=xcorr(pilotSignals(combin_sources(ind_comb,2),:,ind_ues),pilotSignals(combin_sources(ind_comb,1),:,ind_ues));
                    lag_CenterOfMass=lags*(abs(hD_corr).^2)'/sum(abs(hD_corr).^2);
                    allCenterOfMass(ind_comb)=lag_CenterOfMass*obj.sampling_period*obj.c;
                end
                estimatedCentersOfMass(:,ind_ues)=allCenterOfMass;
            end
        end
        
        
        function estimatedDistances2 = locBasedExtract(obj,pilotSignals)
            % this function extracts the LocB features from pilot signals.
            % Accepts input pilotSignals of size S-by-K-by-N where S is the
            % number of sources (or anchor nodes), K is the number of samples
            % per signal, and N is the number of sensor locations.
            % Outputs a matrix of size M-by-N where M is the number of
            % LocB features extracted at one sensor location.
            n_ues=size(pilotSignals,3);
            n_sources=size(pilotSignals,1);
            estimatedDistances2=zeros(n_sources-1,n_ues);
            for ind_ue=1:n_ues
                averaDistAllSources=zeros(n_sources-1,1);
                for ind_sourceHd=2:n_sources
                    [hD_corr, lags]=xcorr(pilotSignals(ind_sourceHd,:,ind_ue),pilotSignals(1,:,ind_ue));
                    [~, Ind_max]=max(abs(hD_corr));
                    tdoaOneSource=lags(Ind_max)*obj.sampling_period;
                    averaDistAllSources(ind_sourceHd-1)=tdoaOneSource*obj.c;
                end
                estimatedDistances2(:,ind_ue)=averaDistAllSources;
            end
            
        end
        
        function distancesDiff = obtainRangeDiffer(obj,source_locs, sensor_locs)
            % this function obtains range differences from source and sensor locations.
            % Accepts inputs: (1) source locations of size 2-by-S where S is the
            % number of sources (or anchor nodes), (2) sensor locations of size 2-by-N
            % where N is the number of sensor locations.
            % Outputs a matrix of size M-by-N where M=S-1 is the number of
            % range differences at one sensor location.
            n_ues=size(sensor_locs,2);
            n_sources=size(source_locs,2);
            assert(size(source_locs, 1)==2);
            assert(size(source_locs, 2)==2);
            distancesDiff=zeros(n_sources-1,n_ues);
            for ind_ue=1:n_ues
                 distancesDiff(:,ind_ue) = vecnorm(sensor_locs(:,ind_ue)-...
                     source_locs(:, 2:end))-vecnorm(sensor_locs(:,ind_ue)-...
                     source_locs(:, 1));
               
            end    
        end
    end
    
end


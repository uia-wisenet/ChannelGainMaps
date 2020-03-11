classdef SpectrumMapSampler
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        UELocations
        UELocationsPairs
        evalGrid_x
        evalGrid_y
        Xenb
        map
        % map should be an object of class
        % MapGenerator
        maxSamplesPerPilot % property not needed anymore.
        pilotNoiseSTD
        powerNoiseSTD
    end
    
    methods
        function measurements_out = sample(obj)
            % extract the received power entries at the map points
            % corresponding to UE locations
            n_sources = size( obj.Xenb , 2);
            n_ues=length(obj.UELocations);
            dist_to_source = zeros(n_sources,n_ues);
            powerAtLocation=zeros(1,n_ues);
            
            for ind_s=1:n_sources
                s_loc=obj.Xenb(:,ind_s);
                for ind_ue = 1:n_ues
                    dist_to_source(ind_s,ind_ue) = norm(obj.UELocations(:,ind_ue)-s_loc);
                    %p_received(ind_ue)= obj.map(find(x_coord==obj.UELocations(1,ind_ue)),find(y_coord==obj.UELocations(2,ind_ue)));
                    powerAtLocation(ind_ue) = obj.map(...
                        obj.evalGrid_x==obj.UELocations(1, ind_ue) & obj.evalGrid_y==obj.UELocations(2, ind_ue));
                    % TODO: this will work OK when the UE locations are a subset
                    % of the evaluation grid. Check that it actually is
                end
                
            end
            % add noise
            noisy_dist_to_source=dist_to_source+obj.locationNoiseSTD*randn(n_sources,n_ues);
            receivedPower = powerAtLocation + obj.powerNoiseSTD*randn(1,n_ues);
            measurements_out=[noisy_dist_to_source;receivedPower]; % TODO think about changing the format of this data
        end
        
        function measurements_out=sampleGivenDistances(obj,distances)
            % distances is an NxNxNs tensor
            % extract the received power entries at the map points
            % corresponding to UE locations
            n_sources = size(distances,3);
            n_ues=length(obj.UELocations);
            powerAtLocation=zeros(1,n_ues);
            dist_to_source = zeros(n_sources,n_ues);
            
            for ind_s=1:n_sources
                for ind_ue = 1:n_ues
                    powerAtLocation(ind_ue) = obj.map(...
                        obj.evalGrid_x==obj.UELocations(1, ind_ue) & obj.evalGrid_y==obj.UELocations(2, ind_ue));
                    % TODO: this will work OK when the UE locations are a subset
                    % of the evaluation grid. Check that it actually is
                    my_distances = distances(:,:,ind_s);
                    dist_to_source(ind_s, ind_ue) = my_distances (obj.evalGrid_x==obj.UELocations(1, ind_ue) & obj.evalGrid_y==obj.UELocations(2, ind_ue));
                end
                
            end
            % add noise
            noisy_dist_to_source=dist_to_source+obj.locationNoiseSTD*randn(n_sources,n_ues);
            receivedPower = powerAtLocation + obj.powerNoiseSTD*randn(1,n_ues);
            measurements_out=[noisy_dist_to_source;receivedPower];
%             noisy_features=noisy_dist_to_source;
            % TODO think about changing the format of this data
        end
        
        function [noiseFree_pilots_sources, noisy_measurement_pilots,measurements_power]=sampleGivenPilotSignals(obj,pilots)
            % distances is an NxNxNs tensor
            % extract the received power entries at the map points
            % corresponding to UE locations
            n_sources = size(pilots,1);
            n_ues=length(obj.UELocations);
            powerAtLocation=zeros(1,n_ues);
            my_maxSamplesPerPilot = size(pilots, 2);
            noiseFree_pilots_sources = zeros(n_sources, my_maxSamplesPerPilot, n_ues);   
                for ind_ue = 1:n_ues
                    powerAtLocation(ind_ue) = obj.map(...
                        obj.evalGrid_x==obj.UELocations(1, ind_ue) & obj.evalGrid_y==obj.UELocations(2, ind_ue));
                    % TODO: this will work OK when the UE locations are a subset
                    % of the evaluation grid. Check that it actually is
                    noiseFree_pilots_sources(:,:, ind_ue) = pilots(:,:,obj.evalGrid_x==obj.UELocations(1, ind_ue) & obj.evalGrid_y==obj.UELocations(2, ind_ue));
                end
                
            % add noise
            noisy_measurement_pilots=noiseFree_pilots_sources+obj.pilotNoiseSTD/sqrt(2)*(randn(n_sources,obj.maxSamplesPerPilot,n_ues)+1i*randn(n_sources,obj.maxSamplesPerPilot,n_ues));
            measurements_power = powerAtLocation + obj.powerNoiseSTD*randn(1,n_ues);
%             noisy_features=noisy_dist_to_source;
            % TODO think about changing the format of this data
        end
        
        function [noiseFree_pilots_sources,noisy_measurement_pilots] = sampleGivenPilotSignals_CG(obj,pilots)
            % pilots is an 4D tensor
            % extract the pilot signals at the points
            % corresponding to UE location pairs
            n_sources = size(pilots,1);
            n_ues_pairs=length(obj.UELocationsPairs);
            noiseFree_pilots_sources = zeros(n_sources,obj.maxSamplesPerPilot,n_ues_pairs,2);   
            for ind_ue = 1:n_ues_pairs
                noiseFree_pilots_sources(:,:, ind_ue,1) = pilots(:,:,abs(obj.evalGrid_x-obj.UELocationsPairs(1, ind_ue))< 1e-12 & abs(obj.evalGrid_y-obj.UELocationsPairs(2, ind_ue))< 1e-12);
                noiseFree_pilots_sources(:,:, ind_ue,2) = pilots(:,:,abs(obj.evalGrid_x-obj.UELocationsPairs(3, ind_ue))< 1e-12 & abs(obj.evalGrid_y-obj.UELocationsPairs(4, ind_ue))< 1e-12);
            end
            noisy_measurement_pilots=noiseFree_pilots_sources+obj.pilotNoiseSTD/sqrt(2)*(randn(n_sources,obj.maxSamplesPerPilot,n_ues_pairs)+1i*randn(n_sources,obj.maxSamplesPerPilot,n_ues_pairs));
%             noisy_features=noisy_dist_to_source;
        end
    end
end




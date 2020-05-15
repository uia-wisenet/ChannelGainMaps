classdef ChannelGainGenerator
    properties (Constant)
        c = 3e8;
    end
    
    methods
        channel_gain = calculateCGBetween(obj, v_xy_tx, v_xy_rx)
    end
    
    methods (Static)    
        function [h_D,samplingTime] = digitalImpulseResponse(alpha, delays, sampling_period, f)
            % calculate digital impulse response given the delays and powers of
            % the analog (time-continuous) channel impulse response
            N_rays=length(alpha);  % Set N_rays to 3 while simulating with the multiwall generator; 6 with the streetcanyon generator
            maximumSamplingTime=ceil(max(delays)/sampling_period); % 10 added is for ensuring some margin, so that all rays are considered
            samplingTime=0:maximumSamplingTime;
            impulseAllRays=zeros(length(samplingTime),N_rays);
            %sampling frequency in MHz
            for ind_ray=1:N_rays
                impulseAllRays(:,ind_ray)=alpha(ind_ray)*exp(-1i*2*pi*f*delays(ind_ray))*sinc(samplingTime-delays(ind_ray)/sampling_period);
            end
            h_D=sum(impulseAllRays,2);
        end
    end
        
end
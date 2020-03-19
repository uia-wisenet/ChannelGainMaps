classdef WangLocationEstimator
    properties 
        
        Xenb % 2xN-matrix containing the source positions. The reference source is in column 1.
        rho=15 % rho: upper bound on the NLOS bias
       
        
    end
    
    methods
        
        function [estimatedLocationRobustSDR, eig_vals_ratios] = estimateLocationRobustSDR(obj, measurements)
            dim=2;
            n_ues=size(measurements,2);
            ue_loc_calc_RobustSDR=zeros(dim,n_ues);
            eig_vals_ratios=zeros(1,n_ues);
            for ind_ue=1:n_ues
                [ue_loc_calc_RobustSDR(:,ind_ue), eig_vals_ratios(ind_ue)]=obj.solve_robust_loc(measurements(:,ind_ue));
            end
            estimatedLocationRobustSDR=ue_loc_calc_RobustSDR;
        end
        
        function[x_out, eig_vals_ratio] = solve_robust_loc(obj, ...
                v_d)
            % Solves the SDP Relaxation of the Robust Location Estimation
            % by Wang et al.
            % Inputs:
            % v_d: (N-1)-vector (N=number of sources) containing the
            %      measured TDOAs (assumed corrupted with small Gaussian 
            %      noise and big n_sourcesLOS bias)
            
            n_sources = length(v_d)+1;
            b = zeros(n_sources-1, 1);
            a_= zeros(n_sources+1, n_sources-1);
            for ind_s_tdoa = 1:n_sources-1
                b(ind_s_tdoa) = -v_d(ind_s_tdoa)^2-norm(obj.Xenb(:,ind_s_tdoa+1))^2+norm(obj.Xenb(:,1))^2;
                a_(:,ind_s_tdoa)=[ 2*(obj.Xenb(:,1)-obj.Xenb(:,ind_s_tdoa+1));
                          zeros(ind_s_tdoa-1,1);
                          -2*v_d(ind_s_tdoa);
                          zeros(n_sources-ind_s_tdoa-1,1)   ];
            end
            [y, eig_vals_ratio]= sdp_ce(obj.Xenb, v_d, b, n_sources, obj.rho, a_);
            x_out = y(1:2);
            v_r_out = y(3:end);
        end
    end
end
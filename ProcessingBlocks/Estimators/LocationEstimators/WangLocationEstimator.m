classdef WangLocationEstimator < LocationEstimator_tdoa
    properties
        
        %Xenb %(in superclass)
        % 2xN-matrix containing the source positions(reference source in
        % column 1)
                
        param_rho=10 % rho: upper bound on the NLOS bias
        % This is a robust location estimator
        
    end
    
    methods
        
%         function estimatedLocationRobustSDR = estimateLocationRobustSDR(obj, measurements)
%             dim=2;
%             n_ues=size(measurements,2);
%             ue_loc_calc_RobustSDR=zeros(dim,n_ues);
%             parfor ind_ue=1:n_ues
%                 ue_loc_calc_RobustSDR(:,ind_ue)=obj.solve_robust_loc(measurements(:,ind_ue));
%             end
%             estimatedLocationRobustSDR=ue_loc_calc_RobustSDR;
%         end
        
        function v_estimatedLocation = estimateOneLocation(obj,  v_d)
            % Solves the SDP Relaxation of the Robust Location Estimation
            % by Wang et al.
            % Inputs:
            % v_d: (N-1)- vector (N=number of sources) containing the
            %      measured TDOAs (assumed corrupted with small Gaussian 
            %      noise and big NLOS bias)
            % Xenb(property of superclass LocationEstimator_tdoa): 
            %      2xN-matrix containing the source positions. The
            %      reference source is in column 1.
            % rho(property of superclass LocationEstimator_tdoa): 
            %      upper bound on the NLOS bias
            
            
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
            % Call method implemented by Gang Wang, see original file:
            % sdp_ce.m
            y= sdp_ce(obj.Xenb, v_d, b, n_sources, obj.param_rho, a_);
            v_estimatedLocation = y(1:2);
            % v_r_out = y(3:end);
            
%             % Luismi code
%             %============
%             arle = ARSRobustLocationEstimator();
%             arle.param_rho = obj.rho;
%             m_s= obj.Xenb;
%             [x_out, ~, ~] = arle.solveRelaxed(v_d, m_s);
        end
    end
end
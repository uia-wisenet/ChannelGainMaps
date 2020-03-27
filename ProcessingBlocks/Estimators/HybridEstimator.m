classdef HybridEstimator < Estimator
    % HybridEstimator class
    %   This class allows to train and query a LF and LB kernel
    %   machine desgined to estimate the channel gain from features of both types 
    %   (location-free features annd location features) received at sensors pairs.
    
    properties
        h_kernelLF % function handle defining the LF kernel
        h_kernelLB % function handle defining the LF kernel
        regularizationParameterLF % lambda for LocFree machine  (KRR)
        regularizationParameterLB % lambda for LocBased machine (KRR)
        max_itrs_alternating =40;

    end
    
    methods
        function [m_coefsOut, intercept, v_d_out] = train(obj, ...
                t_locFeatures, t_estimatedLocation, t_locErrors, v_channelForPairs)
            % Train hybrid Channel-gain estimator.
            
            % GIVEN:    locF features, estimated locations, and 
            % OPTIMIZE: the coefficients of the kernel machines used by the
            %           hybrid estimator of the channel gain (as a linear 
            %           combination of the outputs of the latter machines)
            
            % Inputs: 
            % t_locFeatures: MxNx2-tensor containing locF features
            %               (where M is the number of features at one
            %               sensor location, N is the number of pairs 
            %               for training, and 2 stands for pair)
            % t_estimatedLocation: 2xN-matrix containing the estimated
            %               locations (from a LocationEstimator) for each
            %               training pair
            % v_channelForPairs: N-vector containing channel gain for each 
            %               training pair
            
            % Outputs:
            % m_coefsOut: Nx2 matrix containing the coefficients resulting
            %             from training (N = number of training pairs, 2 
            %             stands for the two machines (location-based and
            %             location-free)
            % intercept:  scalar intercept (mean of training channel gains)
            % v_d_diag:   N-vector containing coefficients at the diagonal 
            %             of matrix D (optimal weights for the linear 
            %             combination at the training points 
            
            %%
            % Dimensionalities and coherence check
            n_ues= length(v_channelForPairs); % number of (training) UEs
            d = size(t_estimatedLocation, 1); % dim. of location vectors
            assert(d==2 || d==3);             % 2D or 3D location
            assert(size(t_estimatedLocation, 2)==n_ues);
            
            disp('Building kernel matrices row by row')
            m_Ke_p = zeros(n_ues, n_ues);
            m_Ke_l = zeros(n_ues, n_ues);           
            ltc = LoopTimeControl(n_ues);
            for ind_ue = 1:n_ues
                m_row_inputs_to_kernels = repmat([t_locFeatures(:, ind_ue, 1);...
                    t_locFeatures(:,ind_ue,2)], [1 n_ues]);
                m_row_inputs_to_kernels_2 = repmat([t_estimatedLocation(:, ind_ue, 1);...
                    t_estimatedLocation(:,ind_ue,2)], [1 n_ues]);
                
                v_my_row = feval(obj.h_kernelLF, m_row_inputs_to_kernels, ...
                    [t_locFeatures(:,:,1);t_locFeatures(:,:,2)]);
                v_my_row_2 = feval(obj.h_kernelLB, m_row_inputs_to_kernels_2, ...
                    [t_estimatedLocation(:,:,1);t_estimatedLocation(:,:,2)]);
                m_Ke_p(ind_ue,:) = v_my_row;
                m_Ke_l(ind_ue,:) = v_my_row_2;
                ltc.go(ind_ue);
            end
            % keyboard %!
            intercept = mean(v_channelForPairs);
            v_zm_channelForPairs=v_channelForPairs-intercept; % zero mean
            
            %% Alternating minimization 
            % (matrix D and the kernel weights in \gamma: see document sent to Seung Jun)  
            [~, indices_t] = sort(t_locErrors(1,:,1));
            [~, indices_r] = sort(t_locErrors(1,:,2));
            norm_dvect=zeros(1,obj.max_itrs_alternating);
            norm_dvecr=zeros(1,obj.max_itrs_alternating);
            v_d_t = rand(n_ues, 1);
            v_d_r = rand(n_ues, 1);
                     
            ltc = LoopTimeControl(obj.max_itrs_alternating);
            for ind_iters=1:obj.max_itrs_alternating  
                % Minimization 1: Given v_d_diag, solve for kernel machine coefficients
                v_coefficients_LF=((eye(n_ues)-(diag(v_d_t)+diag(v_d_r)))*m_Ke_p+((n_ues*obj.regularizationParameterLF)*eye(n_ues)))\v_zm_channelForPairs';
                v_coefficients_LB=((diag(v_d_t)+diag(v_d_r))*m_Ke_l+((n_ues*obj.regularizationParameterLB)*eye(n_ues)))\v_zm_channelForPairs';
                
                % Minimization 2: Given v_a and v_b, solve for v_d
                v_b = m_Ke_l*v_coefficients_LB - m_Ke_p*v_coefficients_LF;
                v_a = v_zm_channelForPairs'- m_Ke_p*v_coefficients_LF;
                % d_vec = diag(b_vec)\a_vec; %?
                % solve using CVX
                cvx_begin quiet
                    variable v_d_t(n_ues, 1);
                    variable v_d_r(n_ues, 1);
                    minimize(norm(v_a - (diag(v_b) * (v_d_t + v_d_r))));
                    subject to
                        0 <= v_d_t + v_d_r;
                        v_d_t + v_d_r <= 1;
                        v_d_t(indices_t(2:end)) >= v_d_t(indices_t(1:end-1));
                        v_d_r(indices_r(2:end)) >= v_d_r(indices_r(1:end-1));        
                cvx_end
                
                % keeping track of the norm of v_d_diag
                norm_dvect(ind_iters)= norm(v_d_t); %if this line 
                norm_dvecr(ind_iters)= norm(v_d_r);
                %  is needed, maybe it is better to move it to the beginning of the loop
                ltc.go(ind_iters);
            end
%%            
%             plots (should be traced outside of the processing block)
            close all
            figure(4)
            plot(norm_dvect+norm_dvecr,'b-o', 'LineWidth',2);   grid on
            xlabel('Interation index k',  'Interpreter', 'latex')
            ylabel('\boldmath $\vert\! \vert$ $\! d_t + d_r $ $\!\vert\!\vert$', 'Interpreter', 'latex')
            
            v_d_out = v_d_t + v_d_r;
            m_coefsOut = [v_coefficients_LF, v_coefficients_LB]; % \gamma
        end
        
        function predictedMeasurements = estimateGivenEstimatedDistancesAndLoc(obj, ...
                m_coefs, t_extractedLocfreeFeaturesToConsider, t_trainingLocations, v_trainingMeasurements,...
                v_funct_weights, test_locF_features, estimatedLocationVal)
            % Estimate channel gain using Hybrid method (locF + locB
            % cartography)
            
            % Inputs: (where M = no. locF features per sensor location, 
            %                N = no. pairs for evaluation;
            %                N_train = no. pairs at training 
            %                d = dimensionality of location space)
            % m_coefs:              N_train x 2- matrix containing 
            %                       coefficients gamma=[alpha, beta]
            % t_extractedLocfreeFeaturesToConsider: 
            %                       MxN_trainx2- tensor containing the
            %                       locFree features of the TRAINING data; 
            %                       2 stands for pair
            % t_trainingLocations:  dxN_trainx2- tensor, estimated locations
            %                       of the TRAINING data
            % v_trainingMeasurements: N_train- vector, measured channel gains of
            %                       TRAINING data
            % v_funct_weights:      N-vector with weights (externally
            %                       calculated for the query/test points)
            % test_locF_features:   MxNx2- array, locF features of query/test points
            % estimatedLocationVal  dxNx2- array, estimated locations of query/test points
            
            % Outputs:
            % predictedMeasurements: N-vector with predicted channel gains
            %                       for query/test points
            
            b_hybrid_option = true; 
            ed_dim = size(test_locF_features); % M x N x 2
            n_uesPair_tr = length(v_trainingMeasurements); % N_train
            predictedMeasurements = zeros(1,ed_dim(2));
            d = size(t_trainingLocations, 1);
            for ind_ue_val = 1:ed_dim(2)
                v_dist_to_source_input_to_ker_tx = ... %(such a long name could be shortened)
                    [test_locF_features(:,ind_ue_val,1); test_locF_features(:,ind_ue_val,2)];
                v_queryLocations=[estimatedLocationVal(:,ind_ue_val,1);estimatedLocationVal(:,ind_ue_val,2)];
                m_features_LF = zeros(ed_dim(1)*2, n_uesPair_tr); %2*M x N_train
                m_features_LF(1:ed_dim(1), :)     = t_extractedLocfreeFeaturesToConsider(:,:,1);
                m_features_LF(ed_dim(1)+1:end, :) = t_extractedLocfreeFeaturesToConsider(:,:,2);
                m_features_LB = zeros(d*2, n_uesPair_tr); %2*d x N_train
                m_features_LB(1:d, :)    = t_trainingLocations(:,:,1);
                m_features_LB(d+1:end,:) = t_trainingLocations(:,:,2);
                row_kernel_LF = feval(obj.h_kernelLF, v_dist_to_source_input_to_ker_tx,  m_features_LF);
                row_kernel_LB = feval(obj.h_kernelLB, v_queryLocations,                  m_features_LB);
                if b_hybrid_option==true
                    predictedMeasurements(ind_ue_val)=((1-v_funct_weights').*row_kernel_LF)*m_coefs(:,1)+...
                        (v_funct_weights'.*row_kernel_LB)*m_coefs(:,2)+mean(v_trainingMeasurements);
                else %if b_hybrid_option is false, this is just an additive model
                    predictedMeasurements(ind_ue_val)=row_kernel_LF*m_coefs(:,1)+...
                        row_kernel_LB*m_coefs(:,2)+mean(v_trainingMeasurements);
                end
                
            end
        end
        
    end %methods
end




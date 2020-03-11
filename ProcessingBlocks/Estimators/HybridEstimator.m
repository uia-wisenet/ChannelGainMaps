classdef HybridEstimator < Estimator
    % HybridEstimator class
    %   This class allows to train and query a LF and LB kernel
    %   machine desgined to estimate the channel gain from location
    %   features received at sensors pairs.
    
    properties
        kernelLF % function handle defining the LF kernel
        kernelLB % function handle defining the LF kernel
        regularizationParameterLF % lambda (ridge regression)
        regularizationParameterLB
    end
    
    methods
        function [coefficientsOut, intercept, d_vec_diag] = train(obj, LocFeatures, estimatedLocation, channelForPairs)
            % given a 3D tensor including the locF features (M-by-N-by-2 array
            % where M is the number of features at one sensor location,
            % N is the number of pairs for training, and 2 stands for pair ) and
            % channel gains, optimize the coefficients
            % of the kernel machine that estimates the channel gain.
            max_itrs =10;
            n_ues= length(channelForPairs);
            d = size(estimatedLocation, 1);
            assert(d==2 || d==3);
            assert(size(estimatedLocation, 2)==n_ues);
            
            disp('Building array row by row')
            Ke_p = zeros(n_ues, n_ues);
            Ke_l = zeros(n_ues, n_ues);
            ltc = LoopTimeControl(n_ues);
            for i = 1:n_ues
                m_row_inputs_to_kernels = repmat([LocFeatures(:, i, 1);...
                    LocFeatures(:,i,2)], [1 n_ues]);
                m_row_inputs_to_kernels_2 = repmat([estimatedLocation(:, i, 1);...
                    estimatedLocation(:,i,2)], [1 n_ues]);
                
                my_row = feval(obj.kernelLF, m_row_inputs_to_kernels, ...
                    [LocFeatures(:,:,1);LocFeatures(:,:,2)]);
                my_row_2 = feval(obj.kernelLB, m_row_inputs_to_kernels_2, ...
                    [estimatedLocation(:,:,1);estimatedLocation(:,:,2)]);
                Ke_p(i,:) = my_row;
                Ke_l(i,:) = my_row_2;
                ltc.go(i);
            end
            
            %             keyboard
            intercept = mean(channelForPairs);
            channelForPairs=channelForPairs-intercept;
            
            %Alternate minimization (matrix D and the kernel weights in the vector \gamma: see document sent to Seung Jun)
            
            norm_dvec=zeros(1,max_itrs);
            ltc = LoopTimeControl(max_itrs);
            d_vec_diag_init =rand(n_ues, 1);
            d_vec_diag = d_vec_diag_init;
            ind_iters =1;
            while (ind_iters <= max_itrs)
                coefficients_LF=((eye(n_ues)-diag(d_vec_diag))*Ke_p+((n_ues*obj.regularizationParameterLF)*eye(n_ues)))\channelForPairs';
                coefficients_LB=(diag(d_vec_diag)*Ke_l+((n_ues*obj.regularizationParameterLB)*eye(n_ues)))\channelForPairs';
                b_vec=Ke_l*coefficients_LB - Ke_p*coefficients_LF;
                a_vec=channelForPairs'- Ke_p*coefficients_LF;
                %                 d_vec = diag(b_vec)\a_vec;
                cvx_begin
                variable d_vec(n_ues, 1)
                minimize(norm(a_vec- diag(b_vec)* d_vec))
                subject to
                for ind_ue=1:n_ues
                    0<=d_vec(ind_ue);
                    d_vec(ind_ue)<= 0.1;
                end
                cvx_end
                norm_dvec(ind_iters)= norm(d_vec_diag);
                d_vec_diag = d_vec;
                ltc.go(i);
                ind_iters = ind_iters +1;
            end
            %             d_vec_diag=d_vec_diag/norm(d_vec_diag);
            %             plot(norm_dvec,'b-o', 'LineWidth',2); grid on
            %             xlabel('Interation index k',  'Interpreter', 'latex')
            %             ylabel('\boldmath $\vert\! \vert$ $\! \sigma(e)$ $\!\vert\!\vert$', 'Interpreter', 'latex')
            
            
            coefficientsOut = [coefficients_LF, coefficients_LB];
        end
        
        function predictedMeasurements = estimateGivenEstimatedDistancesAndLoc(obj, coefficients, extractedLocfreeFeaturesToConsider, trainingLocations, trainingMeasurements,...
                funct_weights,test_locF_features, estimatedLocationVal)
            % Estimate channel gain using Location Free cartography.
            % Inputs: test_locF_features M-by-N-by-2 array
            % where M is the number of features at one sensor location
            % N is the number of pairs where the kernel machine is
            % evaluated,  and 2 stands for pair
            hybrid_option=true;
            ed_dim = size(test_locF_features);
            n_uesPair_tr=length(trainingMeasurements);
            predictedMeasurements=zeros(1,ed_dim(2));
            d = size(trainingLocations, 1);
            for ind_ue_val = 1:ed_dim(2)
                dist_to_source_input_to_ker_tx=[test_locF_features(:,ind_ue_val,1); test_locF_features(:,ind_ue_val,2)];
                queryLocations=[estimatedLocationVal(:,ind_ue_val,1);estimatedLocationVal(:,ind_ue_val,2)];
                m_features_now = zeros(ed_dim(1)*2, n_uesPair_tr);
                m_features_now(1:ed_dim(1), :)     = extractedLocfreeFeaturesToConsider(:,:,1);
                m_features_now(ed_dim(1)+1:end, :) = extractedLocfreeFeaturesToConsider(:,:,2);
                m_features_now_loc = zeros(d*2, n_uesPair_tr);
                m_features_now_loc(1:d, :)    = trainingLocations(:,:,1);
                m_features_now_loc(d+1:end,:) = trainingLocations(:,:,2);
                row_kernel = feval(obj.kernelLF, dist_to_source_input_to_ker_tx, ...
                    m_features_now);
                row1_kernel_loc = feval(obj.kernelLB, queryLocations, m_features_now_loc);
                if hybrid_option==true
                    predictedMeasurements(ind_ue_val)=((1-funct_weights').*row_kernel)*coefficients(:,1)+...
                        (funct_weights'.*row1_kernel_loc)*coefficients(:,2)+mean(trainingMeasurements);
                else
                    predictedMeasurements(ind_ue_val)=row_kernel*coefficients(:,1)+...
                        row1_kernel_loc*coefficients(:,2)+mean(trainingMeasurements);
                end
                
            end
        end
        
    end %methods
end




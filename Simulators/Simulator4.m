classdef Simulator4 < Simulator3
    properties        
        n_monteCarloRuns = 10;
        n_train = 100;
        n_test  = 400;
        b_inParallel = 0; % use parfor
    end
    
    methods
        function str_NMSE = ...
                simulateMonteCarlo ( obj, str_datasetMultipleRealizations_in)
            str_m = str_datasetMultipleRealizations_in;
            n_realizationsInDataset = size(str_m.av_noisy_channelGains, 3);
            s_broadCastFields = ["m_locations", "m_pairs", "v_trueGains", ...
                "datasetGen", "m_grid_x", "m_grid_y"];
            s_stochasticFields = ["m_features_LF", "v_noisy_channelGains", ...
                "m_estimatedLocations", "v_locUncertainties"];
            
            my_randperm = randperm(n_realizationsInDataset);
            for i_r = obj.n_monteCarloRuns:-1:1
                my_str = struct;
                for i_name = 1:length(s_broadCastFields)
                    s_fieldName = s_broadCastFields(i_name);
                    my_str.(s_fieldName) = str_m.(s_fieldName);
                end
                for i_name = 1:length(s_stochasticFields)
                    s_fieldName = s_stochasticFields(i_name);
                    t_source = str_m.("a"+s_fieldName);
                    if obj.n_monteCarloRuns<=n_realizationsInDataset
                        my_str.(s_fieldName) = t_source(:,:,my_randperm(i_r));
                    else
                        v_size = size(t_source);
                        m_value = zeros(v_size(1:2));
                        m_realizationIndices = randi(...
                            n_realizationsInDataset, v_size(1:2));
                        v_indices = zeros(prod(v_size(1:2)),1);
                        for i_i=1:prod(v_size(1:2))
                            v_indices(i_i) = sub2ind(...
                                [prod(v_size(1:2)), n_realizationsInDataset], ...
                                i_i, m_realizationIndices(i_i));
                        end
                        m_value(:) = t_source(v_indices);
                        
                        my_str.(s_fieldName) = m_value;
                    end
                    
                end 
                str_oneRealization(i_r) = my_str;
            
            end
            
            n_allPairs = size(str_oneRealization(1).m_pairs, 1);
            if obj.b_inParallel
                disp 'Simulator4: Running in parallel.'
                parfor i_r = 1:obj.n_monteCarloRuns
                    v_pairs_rand  = randperm(n_allPairs, obj.n_train + obj.n_test)';
                    v_pairs_train = v_pairs_rand(1:obj.n_train);
                    v_pairs_test  = v_pairs_rand(obj.n_train + (1:obj.n_test));
                    [~,~,~,str_NMSE_num(i_r), str_NMSE_den(i_r)] = ...
                        obj.simulate(str_oneRealization(i_r), ...
                        v_pairs_train, v_pairs_test);
                end
            else
                for i_r = 1:obj.n_monteCarloRuns
                    v_pairs_rand  = randperm(n_allPairs, obj.n_train + obj.n_test)';
                    v_pairs_train = v_pairs_rand(1:obj.n_train);
                    v_pairs_test  = v_pairs_rand(obj.n_train + (1:obj.n_test));
                    [~,~,~,str_NMSE_num(i_r), str_NMSE_den(i_r)] = ... 
                        obj.simulate(str_oneRealization(i_r), ...
                        v_pairs_train, v_pairs_test); %#ok<AGROW>
                end
            end
            % calculating the NMSEs and finishing
            s_eNames = string(fieldnames(str_NMSE_num(1)));
            for i_e = 1:length(s_eNames)
                v_num = []; v_den = [];
                for ii =length(str_NMSE_num):-1:1
                    v_num(ii) = str_NMSE_num(ii).(s_eNames(i_e));
                    v_den(ii) = str_NMSE_den(ii).(s_eNames(i_e));
                end
                str_NMSE.(s_eNames(i_e)) = mean(v_num)/mean(v_den);
            end            
        end
    end
end

classdef ChannelGainExperiments < ExperimentFunctionSet
        
    methods
        %% These methods are meant to prepare datasets only        
        function F = experiment_1010(obj, niter)
            % This experiment function only creates a dataset to be used in
            % experiment_1020
            
            % Aiming to generate the same environment as in:
            % LocFCGCartogrExperiments.experiment_4020
            
            rng(1)
            
            x_wall_file_name='./modelFiles/x_coord_4walls.mat';
            y_wall_file_name='./modelFiles/y_coord_4walls.mat';
            selectedWalls=2:4; %!
                        
            gridSize = [5 5 ]; %! 20 20
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            samplingPeriod = 1/receiverBandwidth;
            maxSamplesPerPilot = 10;
            
            my_datasetGen = DatasetGenerator;
            my_datasetGen.b_syntheticLocEstimate = 1;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator2(selectedWalls, ...
                x_wall_file_name, y_wall_file_name, ...
                gridSize, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 0;% 1e-5; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 0.5;  %dB
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
            my_datasetGen.locEstimator = WangLocationEstimator;
            my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);   x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);   y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x2, x1, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            
            % the Mambo line:
            str_dataset = my_datasetGen.generate(m_locations, m_allPairs);
            
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
            
            save (['datasets' filesep 'dataset_ChannelGain_' whichExp], '-struct', 'str_dataset');
            
            figure(999); 
            my_datasetGen.generator.plot_environment;
            F = GFigure.captureCurrentFigure;

        end
        function F = experiment_1030(obj, niter)
            % This experiment function creates a dataset to be used in
            % experiment_2030, and also creates the variables needed
            % to create an animation
            
            % Aiming to generate the same environment as in:
            % LocFCGCartogrExperiments.experiment_4020
            
            rng(1)
            
            x_wall_file_name='./modelFiles/x_coord_4walls.mat';
            y_wall_file_name='./modelFiles/y_coord_4walls.mat';
            selectedWalls=2:4; %!
                        
            gridSize = [10 10];
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            samplingPeriod = 1/receiverBandwidth;
            maxSamplesPerPilot = 10;
            
            my_datasetGen = DatasetGenerator;
            my_datasetGen.b_syntheticLocEstimate = 1;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator2(selectedWalls, ...
                x_wall_file_name, y_wall_file_name, ...
                gridSize, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 1e-5; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 2;  %dB
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
            my_datasetGen.locEstimator = WangLocationEstimator;
            my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);   
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);   
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            n_allPairs = size(m_allPairs, 1);
            
            % for the training set:
            n_trainSamples = 100;
            m_pairsTrain = m_allPairs(randperm(n_allPairs, n_trainSamples),:);

            % for the animation:
            an = Animator;
            an.t_grid_xy = cat(3, m_grid_x, m_grid_y);
            v_rows_frame = 1:gridSize(1);
            an.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
                v_rows_frame, ...
                floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
            an.generator = my_datasetGen.generator;
            m_pairsMap = an.m_pairs();
             
            m_pairsGen = [m_pairsMap; m_pairsTrain];
            v_indicesMap = 1:size(m_pairsMap,1);
            v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
            
            % the Mambo line:
            str_dataset = my_datasetGen.generate(m_locations, m_pairsGen);
            %M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            F = [];
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
            str_dataset.v_indicesMap   = v_indicesMap;
            str_dataset.v_indicesTrain = v_indicesTrain;
            str_dataset.animator   = an;
            
            save (['datasets' filesep 'dataset_ChannelGain_' whichExp], '-struct', 'str_dataset');
%             
%             figure(999); 
%             my_datasetGen.generator.plot_environment;
%             F = GFigure.captureCurrentFigure;

        end
        function F = experiment_1040(obj, niter)
            % This experiment function creates a dataset to be used in
            % experiment_2040, and also creates the variables needed
            % to create an animation
                        
            rng(1)
            
            x_wall_file_name='./modelFiles/x_coord_4walls.mat';
            y_wall_file_name='./modelFiles/y_coord_4walls.mat';
            selectedWalls=2:5; %!
                        
            gridSize = [18 18];
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            samplingPeriod = 1/receiverBandwidth;
            maxSamplesPerPilot = 10;
            
            my_datasetGen = DatasetGenerator;
            my_datasetGen.b_syntheticLocEstimate = 1;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator2(selectedWalls, ...
                x_wall_file_name, y_wall_file_name, ...
                gridSize, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 1e-5; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 1;    % dB
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
            my_datasetGen.locEstimator = WangLocationEstimator;
            my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            n_allPairs = size(m_allPairs, 1);
            
            % for the training set:
            n_trainSamples = 1000;
            m_pairsTrain = m_allPairs(randperm(n_allPairs, n_trainSamples),:);

            % for the animation:
            an = Animator;
            an.t_grid_xy = cat(3, m_grid_x, m_grid_y);
            v_rows_frame = 1:2:gridSize(1);
            an.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
                v_rows_frame, ...
                floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
            an.generator = my_datasetGen.generator;
            m_pairsMap = an.m_pairs();
             
            m_pairsGen = [m_pairsMap; m_pairsTrain];
            v_indicesMap = 1:size(m_pairsMap,1);
            v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
            
            % the Mambo line:
            str_dataset = my_datasetGen.generate(m_locations, m_pairsGen);
            %M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            F = [];
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
            str_dataset.v_indicesMap   = v_indicesMap;
            str_dataset.v_indicesTrain = v_indicesTrain;
            str_dataset.animator   = an;
            
            filename = 'dataset_ChannelGain_1040';
            
            %%%
            save (['datasets' filesep filename], '-struct', 'str_dataset');
%             
%             figure(999); 
%             my_datasetGen.generator.plot_environment;
%             F = GFigure.captureCurrentFigure;

        end
        function F = experiment_1060(obj, niter)
            % This experiment function creates a dataset to be used in
            % experiment_2060, where we want to average over a number of
            % Monte Carlo realizations
                        
            rng(1)
            s_fileName_environment = 'modelFiles/env1.mat';           
            selectedWalls=2:5;
                        
            gridSize = [18 18];
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            samplingPeriod = 1/receiverBandwidth;
            maxSamplesPerPilot = 10;
            
            my_datasetGen = DatasetGenerator2;
            my_datasetGen.b_syntheticLocEstimate = 1;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator3(load(s_fileName_environment), ...
                selectedWalls, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 1e-5; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 1;    % dB
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
            my_datasetGen.locEstimator = WangLocationEstimator;
            my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            n_allPairs = size(m_allPairs, 1);
            
            % for the training set:
            n_trainSamples = 1000;
            m_pairsTrain = m_allPairs(...
                randperm(n_allPairs, n_trainSamples),:);

%             % for the animation:
%             an = Animator;
%             an.t_grid_xy = cat(3, m_grid_x, m_grid_y);
%             v_rows_frame = 1:2:gridSize(1);
%             an.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
%                 v_rows_frame, ...
%                 floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
%             an.generator = my_datasetGen.generator;
%             m_pairsMap = an.m_pairs();
%              
%             m_pairsGen = [m_pairsMap; m_pairsTrain];
%             v_indicesMap = 1:size(m_pairsMap,1);
%             v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
            
            % the Mambo line:
            my_datasetGen.n_realizations = 10;
            str_dataset = my_datasetGen.generate(m_locations, m_pairsTrain);
            %M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            F = [];
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
%             str_dataset.v_indicesMap   = v_indicesMap;
%             str_dataset.v_indicesTrain = v_indicesTrain;
%             str_dataset.animator   = an;
                        
            %%%
            save (['datasets' filesep 'dataset_ChannelGain_' whichExp], '-struct', 'str_dataset');
%             
%             figure(999); 
%             my_datasetGen.generator.plot_environment;
%             F = GFigure.captureCurrentFigure;

        end
        
        function F = experiment_1080(obj, niter)
            % This experiment function creates a dataset to be used in
            % experiment_2080, 
            % new environment with 7 walls, and half of the map has free-
            % space propagation conditions.
            % 10 Realizations to conform Monte Carlo runs
                        
            rng(1)
            
            s_fileName_environment = 'modelFiles/env2.mat';           
            selectedWalls=1:7;
                        
            gridSize = [20 14];
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            samplingPeriod = 1/receiverBandwidth;
            maxSamplesPerPilot = 12;
            
            my_datasetGen = DatasetGenerator2;
            my_datasetGen.b_syntheticLocEstimate = 0;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator3(load(s_fileName_environment), ...
                selectedWalls, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 1e-5; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 2;    % dB
            my_datasetGen.std_syntheticLocationNoise = 7;
            my_datasetGen.snr_locErrors = 20;
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
            my_datasetGen.locEstimator = ARSRobustLocationEstimator;
            my_datasetGen.locEstimator.param_rho = 8;
            my_datasetGen.locEstimator.b_inParallel = 0;
            my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            n_allPairs = size(m_allPairs, 1);
            
            % for the training set:
            n_trainSamples = 3000;
            m_pairsTrain = m_allPairs(...
                randperm(n_allPairs, n_trainSamples),:);

%             % for the animation:
%             an = Animator;
%             an.t_grid_xy = cat(3, m_grid_x, m_grid_y);
%             v_rows_frame = 1:2:gridSize(1);
%             an.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
%                 v_rows_frame, ...
%                 floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
%             an.generator = my_datasetGen.generator;
%             m_pairsMap = an.m_pairs();
%              
%             m_pairsGen = [m_pairsMap; m_pairsTrain];
%             v_indicesMap = 1:size(m_pairsMap,1);
%             v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
            
            
            my_datasetGen.n_realizations = 10;
            my_datasetGen.b_inParallel = 1;
            % the Mambo line:
            str_dataset = my_datasetGen.generate(m_locations, m_pairsTrain);
            %M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            F = [];
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
%             str_dataset.v_indicesMap   = v_indicesMap;
%             str_dataset.v_indicesTrain = v_indicesTrain;
%             str_dataset.animator   = an;
            
            save (['datasets' filesep 'dataset_ChannelGain_' whichExp], '-struct', 'str_dataset');
%             
%             figure(999); 
%             my_datasetGen.generator.plot_environment;
%             F = GFigure.captureCurrentFigure;

        end
        function F = experiment_1099(obj, niter)
            % This experiment function creates a dataset to be used in
            % experiment_2080, 
            % new environment with 7 walls, and half of the map has free-
            % space propagation conditions.
            % 10 Realizations to conform Monte Carlo runs
                        
            rng(1)
            
            s_fileName_environment = 'modelFiles/env3.mat';           
            selectedWalls=1:7;
                        
            gridSize = [20 14];
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            samplingPeriod = 1/receiverBandwidth;
            maxSamplesPerPilot = 12;
            
            my_datasetGen = DatasetGenerator2;
            my_datasetGen.b_syntheticLocEstimate = 0;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator3(load(s_fileName_environment), ...
                selectedWalls, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 1e-5; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 2;    % dB
            my_datasetGen.std_syntheticLocationNoise = 7;
            my_datasetGen.snr_locErrors = 20;
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
            my_datasetGen.locEstimator = ARSRobustLocationEstimator;
            my_datasetGen.locEstimator.param_rho = 8;
            my_datasetGen.locEstimator.b_inParallel = 0;
            my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            n_allPairs = size(m_allPairs, 1);
            
            % for the training set:
            n_trainSamples = 3000;
            m_pairsTrain = m_allPairs(...
                randperm(n_allPairs, n_trainSamples),:);

%             % for the animation:
%             an = Animator;
%             an.t_grid_xy = cat(3, m_grid_x, m_grid_y);
%             v_rows_frame = 1:2:gridSize(1);
%             an.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
%                 v_rows_frame, ...
%                 floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
%             an.generator = my_datasetGen.generator;
%             m_pairsMap = an.m_pairs();
%              
%             m_pairsGen = [m_pairsMap; m_pairsTrain];
%             v_indicesMap = 1:size(m_pairsMap,1);
%             v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
            
            
            my_datasetGen.n_realizations = 30;
            my_datasetGen.b_inParallel = 1;
            % the Mambo line:
            str_dataset = my_datasetGen.generate(m_locations, m_pairsTrain);
            %M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            F = [];
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
%             str_dataset.v_indicesMap   = v_indicesMap;
%             str_dataset.v_indicesTrain = v_indicesTrain;
%             str_dataset.animator   = an;
            
            save (['datasets' filesep 'dataset_ChannelGain_' whichExp], '-struct', 'str_dataset');
%             
%             figure(999); 
%             my_datasetGen.generator.plot_environment;
%             F = GFigure.captureCurrentFigure;

        end
        function F = experiment_1095(obj, niter)
            % This experiment function creates a dataset to be used in
            % experiment_2095, 
            % new environment with 7 walls, and half of the map has free-
            % space propagation conditions.
            % 10 Realizations to conform Monte Carlo runs
                        
            rng(1)
            
            s_fileName_environment = 'modelFiles/env2.mat';           
            selectedWalls=1:7;
                        
            gridSize = [20 14];
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            samplingPeriod = 1/receiverBandwidth;
            maxSamplesPerPilot = 12;
            
            my_datasetGen = DatasetGenerator2;
            my_datasetGen.b_syntheticLocEstimate = 0;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator3(load(s_fileName_environment), ...
                selectedWalls, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 1e-5; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 2;    % dB
            my_datasetGen.std_syntheticLocationNoise = 7;
            my_datasetGen.snr_locErrors = 20;
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
            my_datasetGen.locEstimator = ARSRobustLocationEstimator;
            my_datasetGen.locEstimator.param_rho = 8;
            my_datasetGen.locEstimator.b_inParallel = 0;
            my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            n_allPairs = size(m_allPairs, 1);
            
            % for the training set:
            n_trainSamples = 2000;
            m_pairsTrain = m_allPairs(...
                randperm(n_allPairs, n_trainSamples),:);

            % for the animation:
            an = Animator;
            an.t_grid_xy = cat(3, m_grid_x, m_grid_y);
            v_rows_frame = 1:2:gridSize(1);
            an.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
                v_rows_frame, ...
                floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
            an.generator = my_datasetGen.generator;
            m_pairsMap = an.m_pairs();
             
            m_pairsGen = [m_pairsMap; m_pairsTrain];
            v_indicesMap = 1:size(m_pairsMap,1);
            v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
            
            
            my_datasetGen.n_realizations = 1;
            my_datasetGen.b_inParallel = 1;
            % the Mambo line:
            str_dataset = my_datasetGen.generate(m_locations, m_pairsGen);
            %M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            F = [];
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
            str_dataset.v_indicesMap   = v_indicesMap;
            str_dataset.v_indicesTrain = v_indicesTrain;
            str_dataset.animator   = an;
            
            save (['datasets' filesep 'dataset_ChannelGain_' whichExp], '-struct', 'str_dataset');
%             
%             figure(999); 
%             my_datasetGen.generator.plot_environment;
%             F = GFigure.captureCurrentFigure;

        end
        function F = experiment_1100(obj, niter)
            % This experiment function creates a dataset to be used in
            % experiment_2100, 
            % new environment with 7 walls, and half of the map has free-
            % space propagation conditions.
            % 10 Realizations to conform Monte Carlo runs
                        
            rng(1)
            
            s_fileName_environment = 'modelFiles/env3.mat';           
            selectedWalls=1:7;
                        
            gridSize = [24 18];
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            samplingPeriod = 1/(10*receiverBandwidth); %! was 1/
            maxSamplesPerPilot = 100; %! was 10
            
            my_datasetGen = DatasetGenerator2;
            my_datasetGen.b_syntheticLocEstimate = 0;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator3(load(s_fileName_environment), ...
                selectedWalls, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 0; %! 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 1e-6; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 2;    % dB
            my_datasetGen.std_syntheticLocationNoise = 7;
            my_datasetGen.snr_locErrors = 20;
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
            my_datasetGen.locEstimator = ARSRobustLocationEstimator;
            my_datasetGen.locEstimator.param_rho = 12;
            my_datasetGen.locEstimator.b_inParallel = 0;
            my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            n_allPairs = size(m_allPairs, 1);
            
            % for the training set:
            n_trainSamples = 2000;
            m_pairsTrain = m_allPairs(...
                randperm(n_allPairs, n_trainSamples),:);

            % for the animation:
            an = Animator;
            an.t_grid_xy = cat(3, m_grid_x, m_grid_y);
            v_rows_frame = 1:2:gridSize(1);
            an.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
                v_rows_frame, ...
                floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
            an.generator = my_datasetGen.generator;
            m_pairsMap = an.m_pairs();
             
            m_pairsGen = [m_pairsMap; m_pairsTrain];
            v_indicesMap = 1:size(m_pairsMap,1);
            v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
            
            
            my_datasetGen.n_realizations = 1;
            my_datasetGen.b_inParallel = 1;
            % the Mambo line:
            str_dataset = my_datasetGen.generate(m_locations, m_pairsGen);
            %M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            F = [];
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
            str_dataset.v_indicesMap   = v_indicesMap;
            str_dataset.v_indicesTrain = v_indicesTrain;
            str_dataset.animator   = an;
            
            save (['datasets' filesep 'dataset_ChannelGain_' whichExp], '-struct', 'str_dataset');
%             
%             figure(999); 
%             my_datasetGen.generator.plot_environment;
%             F = GFigure.captureCurrentFigure;

        end
        
        function F = experiment_1110(obj, niter)
            % This experiment function creates a dataset to be used in
            % experiment_2110, 
            % new environment with 7 walls, and half of the map has free-
            % space propagation conditions.
            % 10 Realizations to conform Monte Carlo runs
                        
            rng(1)
            
            s_fileName_environment = 'modelFiles/env3.mat';           
            selectedWalls=1:7;
                        
            gridSize = [24 18];
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            samplingPeriod = 1/(10*receiverBandwidth);
            maxSamplesPerPilot = 100;
            
            my_datasetGen = DatasetGenerator2;
            my_datasetGen.b_syntheticLocEstimate = 0;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator3(load(s_fileName_environment), ...
                selectedWalls, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 3e-6; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 2;    % dB
            my_datasetGen.std_syntheticLocationNoise = 7;
            my_datasetGen.snr_locErrors = 20;
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
            my_datasetGen.locEstimator = ARSRobustLocationEstimator;
            my_datasetGen.locEstimator.param_rho = 12;
            my_datasetGen.locEstimator.b_inParallel = 0;
            my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            n_allPairs = size(m_allPairs, 1);
            
            % for the training set:
            n_trainSamples = 3000;
            m_pairsTrain = m_allPairs(...
                randperm(n_allPairs, n_trainSamples),:);

%             % for the animation:
%             an = Animator;
%             an.t_grid_xy = cat(3, m_grid_x, m_grid_y);
%             v_rows_frame = 1:2:gridSize(1);
%             an.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
%                 v_rows_frame, ...
%                 floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
%             an.generator = my_datasetGen.generator;
%             m_pairsMap = an.m_pairs();
%              
%             m_pairsGen = [m_pairsMap; m_pairsTrain];
%             v_indicesMap = 1:size(m_pairsMap,1);
%             v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
            
            
            my_datasetGen.n_realizations = 30;
            my_datasetGen.b_inParallel = 1;
            % the Mambo line:
            str_dataset = my_datasetGen.generate(m_locations, m_pairsTrain);
            %M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            F = [];
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
%             str_dataset.v_indicesMap   = v_indicesMap;
%             str_dataset.v_indicesTrain = v_indicesTrain;
%             str_dataset.animator   = an;
            
            save (['datasets' filesep 'dataset_ChannelGain_' whichExp], '-struct', 'str_dataset');
%             
%             figure(999); 
%             my_datasetGen.generator.plot_environment;
%             F = GFigure.captureCurrentFigure;

        end

        function F = experiment_1119(obj, niter)
            % Generates very few traning samples. Intended for representing
            % the location estimation error (pixelwise) in experiment_3119
            rng(1)
            
            s_fileName_environment = 'modelFiles/env3.mat';           
            selectedWalls=1:7;
                        
            gridSize = [20 15];
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            samplingPeriod = 1/(receiverBandwidth);
            maxSamplesPerPilot = 12;
            
            my_datasetGen = DatasetGenerator2;
            my_datasetGen.b_syntheticLocEstimate = 0;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator3(load(s_fileName_environment), ...
                selectedWalls, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 3e-6; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 2;    % dB
            my_datasetGen.std_syntheticLocationNoise = 7;
            my_datasetGen.snr_locErrors = 20;
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
            my_datasetGen.locEstimator = ARSRobustLocationEstimator;
            my_datasetGen.locEstimator.param_rho = 12;
            my_datasetGen.locEstimator.b_inParallel = 0;
            my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            n_allPairs = size(m_allPairs, 1);
            
            % for the training set:
            n_trainSamples = 2;
            m_pairsTrain = m_allPairs(...
                randperm(n_allPairs, n_trainSamples),:);

%             % for the animation:
%             an = Animator;
%             an.t_grid_xy = cat(3, m_grid_x, m_grid_y);
%             v_rows_frame = 1:2:gridSize(1);
%             an.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
%                 v_rows_frame, ...
%                 floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
%             an.generator = my_datasetGen.generator;
%             m_pairsMap = an.m_pairs();
%              
%             m_pairsGen = [m_pairsMap; m_pairsTrain];
%             v_indicesMap = 1:size(m_pairsMap,1);
%             v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
            
            
            my_datasetGen.n_realizations = 1;
            my_datasetGen.b_inParallel = 0;
            % the Mambo line:
            str_dataset = my_datasetGen.generate(m_locations, m_pairsTrain);
            %M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            F = [];
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
%             str_dataset.v_indicesMap   = v_indicesMap;
%             str_dataset.v_indicesTrain = v_indicesTrain;
%             str_dataset.animator   = an;
            
            save (['datasets' filesep 'dataset_ChannelGain_' whichExp], '-struct', 'str_dataset');

        end

        function F = experiment_1129(obj, niter)
            % Generates very few traning samples. Intended for representing
            % the location estimation error (pixelwise) in experiment_3129
            rng(1)
            
            s_fileName_environment = 'modelFiles/env3.mat';           
            selectedWalls=1:7;
                        
            gridSize = [20 15];
            
            carrier_frequency = 800e6;
            receiverBandwidth = 40e6; % Twice as the normal bandwidth
            samplingPeriod = 1/(receiverBandwidth);
            maxSamplesPerPilot = 24;
            
            my_datasetGen = DatasetGenerator2;
            my_datasetGen.b_syntheticLocEstimate = 0;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator3(load(s_fileName_environment), ...
                selectedWalls, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 3e-6; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 2;    % dB
            my_datasetGen.std_syntheticLocationNoise = 7;
            my_datasetGen.snr_locErrors = 20;
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
            my_datasetGen.locEstimator = ARSRobustLocationEstimator;
            my_datasetGen.locEstimator.param_rho = 12;
            my_datasetGen.locEstimator.b_inParallel = 0;
            my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            n_allPairs = size(m_allPairs, 1);
            
            % for the training set:
            n_trainSamples = 2;
            m_pairsTrain = m_allPairs(...
                randperm(n_allPairs, n_trainSamples),:);

%             % for the animation:
%             an = Animator;
%             an.t_grid_xy = cat(3, m_grid_x, m_grid_y);
%             v_rows_frame = 1:2:gridSize(1);
%             an.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
%                 v_rows_frame, ...
%                 floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
%             an.generator = my_datasetGen.generator;
%             m_pairsMap = an.m_pairs();
%              
%             m_pairsGen = [m_pairsMap; m_pairsTrain];
%             v_indicesMap = 1:size(m_pairsMap,1);
%             v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
            
            
            my_datasetGen.n_realizations = 1;
            my_datasetGen.b_inParallel = 1;
            % the Mambo line:
            str_dataset = my_datasetGen.generate(m_locations, m_pairsTrain);
            %M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            F = [];
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
%             str_dataset.v_indicesMap   = v_indicesMap;
%             str_dataset.v_indicesTrain = v_indicesTrain;
%             str_dataset.animator   = an;
            

            save (['datasets' filesep 'dataset_ChannelGain_' whichExp], '-struct', 'str_dataset');
%             
%             figure(999); 
%             my_datasetGen.generator.plot_environment;
%             F = GFigure.captureCurrentFigure;

        end
        
        function F = experiment_1510(obj, niter)
            % Environment 3, synthetic loc estimate, 
            % 
            rng(1)
            
            s_fileName_environment = 'modelFiles/env3.mat';           
            selectedWalls=1:7;
                        
            gridSize = [20 14];
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            samplingPeriod = 1/receiverBandwidth;
            maxSamplesPerPilot = 12;
            
            my_datasetGen = DatasetGenerator2;
            my_datasetGen.b_syntheticLocEstimate = 1;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator3(load(s_fileName_environment), ...
                selectedWalls, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 1e-5; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 2;    % dB
            my_datasetGen.std_syntheticLocationNoise = 7;
            my_datasetGen.snr_locErrors = 20;
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
%             my_datasetGen.locEstimator = ARSRobustLocationEstimator;
%             my_datasetGen.locEstimator.param_rho = 8;
%             my_datasetGen.locEstimator.b_inParallel = 0;
%             my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            n_allPairs = size(m_allPairs, 1);
            
            % for the training set:
            n_trainSamples = 3000;
            m_pairsTrain = m_allPairs(...
                randperm(n_allPairs, n_trainSamples),:);

%             % for the animation:
%             an = Animator;
%             an.t_grid_xy = cat(3, m_grid_x, m_grid_y);
%             v_rows_frame = 1:2:gridSize(1);
%             an.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
%                 v_rows_frame, ...
%                 floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
%             an.generator = my_datasetGen.generator;
%             m_pairsMap = an.m_pairs();
%              
%             m_pairsGen = [m_pairsMap; m_pairsTrain];
%             v_indicesMap = 1:size(m_pairsMap,1);
%             v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
            
            
            my_datasetGen.n_realizations = 30;
            my_datasetGen.b_inParallel = 1;
            % the Mambo line:
            str_dataset = my_datasetGen.generate(m_locations, m_pairsTrain);
            %M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            F = [];
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
%             str_dataset.v_indicesMap   = v_indicesMap;
%             str_dataset.v_indicesTrain = v_indicesTrain;
%             str_dataset.animator   = an;
            
            save (['datasets' filesep 'dataset_ChannelGain_' whichExp], '-struct', 'str_dataset');

            
        end
        function F = experiment_1515(obj, niter)
            % Same settings as 1510, but used to create the
            % animation/storyboard
            rng(1)
            
            s_fileName_environment = 'modelFiles/env3.mat';           
            selectedWalls=1:7;
                        
            gridSize = [20 14];
            
            carrier_frequency = 800e6;
            receiverBandwidth = 20e6;
            samplingPeriod = 1/receiverBandwidth;
            maxSamplesPerPilot = 12;
            
            my_datasetGen = DatasetGenerator2;
            my_datasetGen.b_syntheticLocEstimate = 1;
            my_datasetGen.b_syntheticLocError = 1;
            
            generator_tmp = MultiWallChannelGainGenerator;
            generator_tmp.f = carrier_frequency;
            generator_tmp.sampling_period = samplingPeriod;
            generator_tmp.maxSamplesPerPilot = maxSamplesPerPilot;
            [my_datasetGen.generator, m_source_loc] = ...
                baselineGenerator3(load(s_fileName_environment), ...
                selectedWalls, generator_tmp);
            my_datasetGen.generator.delay_estimation_offset = 2*rand;
            
            my_datasetGen.sampler = SpectrumMapSampler;
            my_datasetGen.sampler.pilotNoiseSTD = 1e-5; % natural units
            my_datasetGen.sampler.powerNoiseSTD = 2;    % dB
            my_datasetGen.std_syntheticLocationNoise = 7;
            my_datasetGen.snr_locErrors = 20;
            my_datasetGen.sampler.maxSamplesPerPilot = maxSamplesPerPilot;

            my_datasetGen.featureExtractor = FeatureExtractor;
            my_datasetGen.featureExtractor.sampling_period = samplingPeriod;
            
%             my_datasetGen.locEstimator = ARSRobustLocationEstimator;
%             my_datasetGen.locEstimator.param_rho = 8;
%             my_datasetGen.locEstimator.b_inParallel = 0;
%             my_datasetGen.locEstimator.Xenb = m_source_loc;
            
            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_l = numel(m_grid_x);
            m_allPairs  = combnk(1:n_l, 2);
            n_allPairs = size(m_allPairs, 1);
            
            % for the training set:
            n_trainSamples = 2000;
            m_pairsTrain = m_allPairs(...
                randperm(n_allPairs, n_trainSamples),:);

            % for the animation:
            an = Animator;
            an.t_grid_xy = cat(3, m_grid_x, m_grid_y);
            v_rows_frame = 1:2:gridSize(1);
            an.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
                v_rows_frame, ...
                floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
            an.generator = my_datasetGen.generator;
            m_pairsMap = an.m_pairs();
             
            m_pairsGen = [m_pairsMap; m_pairsTrain];
            v_indicesMap = 1:size(m_pairsMap,1);
            v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
                       
            my_datasetGen.n_realizations = 1;
            my_datasetGen.b_inParallel = 0;
            % the Mambo line:
            str_dataset = my_datasetGen.generate(m_locations, m_pairsGen);
            %M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            F = [];
            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
            str_dataset.v_indicesMap   = v_indicesMap;
            str_dataset.v_indicesTrain = v_indicesTrain;
            str_dataset.animator   = an;
            
            save (['datasets' filesep 'dataset_ChannelGain_' whichExp], '-struct', 'str_dataset');

            
        end

        %% These methods run simulations using the datasets generated above
        function F = experiment_2010(obj, niter)
            str_dataset = load ('datasets/dataset_ChannelGain_1010');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator3;
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 1;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 14);    
            mySim.v_sigmas_toTryLF  = linspace(85, 200, 15); 
            
            mySim.v_lambdas_toTryLB = logspace(-5, -2,  15);   
            mySim.v_sigmas_toTryLB  = linspace(20, 100, 12);
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator;
            mySim.hybridEstimator.b_debugPlots = 0;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

            n_allPairs = size(str_dataset.m_pairs, 1);
            v_trainTestPairs = randperm(n_allPairs, 120)';
            v_trainPairs = v_trainTestPairs(1:100);
            v_testPairs  = v_trainTestPairs(101:120);
            
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulate(...
                str_dataset, v_trainPairs, v_testPairs);
            str_NMSE
            
            save (['savedResults' filesep 'results_' whichExp]);

            F = GFigure.captureCurrentFigure();

        end
        
        function F = experiment_2030(obj, niter)
            str_dataset = load ('datasets/dataset_ChannelGain_1030');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator3;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0; %!!
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );    
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10); 
            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator;
            mySim.hybridEstimator.b_debugPlots = 0;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            v_trainPairs = str_dataset.v_indicesTrain(:);
            v_mapPairs   = str_dataset.v_indicesMap(:);
            
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulate(...
                str_dataset, v_trainPairs, v_mapPairs);

            M = str_dataset.animator.create( [str_mapEstimates.locFree'...
                 str_mapEstimates.locBased' trueGains(:)])
             
            save (['savedResults' filesep 'results_' whichExp]);
            
            F = GFigure.captureCurrentFigure();

        end

        function F = experiment_2040(obj, niter)
            str_dataset = load ('datasets/dataset_ChannelGain_1040');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator3;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0; %!!
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = 4e-4; %logspace(-4, -2, 10 );    
            mySim.v_sigmas_toTryLF  = 67;   %linspace(50, 130, 10); 
            
            mySim.v_lambdas_toTryLB = 4e-4; %logspace(-4, -2,  10);   
            mySim.v_sigmas_toTryLB  = 25;   %linspace(20, 40, 10);
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 2;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            v_trainPairs = str_dataset.v_indicesTrain(:);
            v_mapPairs   = str_dataset.v_indicesMap(:);
            
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulate(...
                str_dataset, v_trainPairs, v_mapPairs);

%             M = str_dataset.animator.create( permute(cat(3, ...
%                 [str_mapEstimates.locFree'...
%                  str_mapEstimates.locBased' ...
%                  trueGains(:)], ...
%                 [str_mapEstimates.locFree_fromHybrid'...
%                  str_mapEstimates.locBased_fromHybrid' ...
%                  str_mapEstimates.hybrid(:)]), [1 3 2]), ...
%                  ["locFree", "locBased", "true";
%                  "locFreeH", "locBasedH", "hybrid"]);
             
             M = str_dataset.animator.createGivenMatrix([...
                 str_mapEstimates.locFree'...
                 str_mapEstimates.locFree_fromHybrid'...
                 str_mapEstimates.locBased' ...
                 str_mapEstimates.locBased_fromHybrid' ...
                 trueGains(:), ...
                 str_mapEstimates.hybrid(:) ], ...
                 ["locFree", "locBased", "true";
                 "locFreeH", "locBasedH", "hybrid"]);
             
            save (['savedResults' filesep 'results_' whichExp]);
            
            F = GFigure.captureCurrentFigure();

        end

        function F = experiment_2050(obj, niter)
            str_dataset = load ('datasets/dataset_ChannelGain_1040');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator3;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0; %!!
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = 4e-4; %logspace(-4, -2, 10 );    
            mySim.v_sigmas_toTryLF  = 67;   %linspace(50, 130, 10); 
            
            mySim.v_lambdas_toTryLB = 4e-4; %logspace(-4, -2,  10);   
            mySim.v_sigmas_toTryLB  = 25;   %linspace(20, 40, 10);
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 20;
            mySim.hybridEstimator.b_tryToBalance = 1;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            v_trainPairs = str_dataset.v_indicesTrain(:);
            v_mapPairs   = str_dataset.v_indicesMap(:);
            
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulate(...
                str_dataset, v_trainPairs, v_mapPairs);

%             M = str_dataset.animator.create( permute(cat(3, ...
%                 [str_mapEstimates.locFree'...
%                  str_mapEstimates.locBased' ...
%                  trueGains(:)], ...
%                 [str_mapEstimates.locFree_fromHybrid'...
%                  str_mapEstimates.locBased_fromHybrid' ...
%                  str_mapEstimates.hybrid(:)]), [1 3 2]), ...
%                  ["locFree", "locBased", "true";
%                  "locFreeH", "locBasedH", "hybrid"]);
             
             M = str_dataset.animator.createGivenMatrix([...
                 str_mapEstimates.locFree'...
                 str_mapEstimates.locFree_fromHybrid'...
                 str_mapEstimates.locBased' ...
                 str_mapEstimates.locBased_fromHybrid' ...
                 trueGains(:), ...
                 str_mapEstimates.hybrid(:) ], ...
                 ["locFree", "locBased", "true";
                 "locFreeH", "locBasedH", "hybrid"]);
             
            save (['savedResults' filesep 'results_' whichExp]);
            
            F = GFigure.captureCurrentFigure();

        end
        
        function F = experiment_2055(obj, niter)
            str_dataset = load ('datasets/dataset_ChannelGain_1040');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator3;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0; %!!
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = 4e-4; %logspace(-4, -2, 10 );    
            mySim.v_sigmas_toTryLF  = 67;   %linspace(50, 130, 10); 
            
            mySim.v_lambdas_toTryLB = 4e-4; %logspace(-4, -2,  10);   
            mySim.v_sigmas_toTryLB  = 25;   %linspace(20, 40, 10);
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 20;
            mySim.hybridEstimator.b_tryToBalance = 1;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            v_trainPairs = str_dataset.v_indicesTrain(:);
            v_mapPairs   = str_dataset.v_indicesMap(:);
            
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulate(...
                str_dataset, v_trainPairs, v_mapPairs);

%             M = str_dataset.animator.create( permute(cat(3, ...
%                 [str_mapEstimates.locFree'...
%                  str_mapEstimates.locBased' ...
%                  trueGains(:)], ...
%                 [str_mapEstimates.locFree_fromHybrid'...
%                  str_mapEstimates.locBased_fromHybrid' ...
%                  str_mapEstimates.hybrid(:)]), [1 3 2]), ...
%                  ["locFree", "locBased", "true";
%                  "locFreeH", "locBasedH", "hybrid"]);
             ms_titles = ["locFree", "locBased", "true";
                 "locFreeH", "locBasedH", "hybrid"];
             M = str_dataset.animator.createGivenMatrix([...
                 str_mapEstimates.locFree'...
                 str_mapEstimates.locFree_fromHybrid'...
                 str_mapEstimates.locBased' ...
                 str_mapEstimates.locBased_fromHybrid' ...
                 trueGains(:), ...
                 str_mapEstimates.hybrid(:) ], ...
                 ms_titles);
             
            save (['savedResults' filesep 'results_' whichExp]);
            
            figure;
            vs_titles = ["locFree", "locBased", ...
                 "locFreeH", "locBasedH", "hybrid", "true"];
            str_dataset.animator.storyBoard([...
                 str_mapEstimates.locFree'...
                 str_mapEstimates.locBased' ...
                 str_mapEstimates.locFree_fromHybrid' ...
                 str_mapEstimates.locBased_fromHybrid' ...
                 str_mapEstimates.hybrid(:) ...
                 trueGains(:)], ...
                 vec(vs_titles), 3:3:9);
            
            F = GFigure.captureCurrentFigure();

        end

        % Monte Carlo simulations
        function F = experiment_2060(obj, niter)
            str_dataset = load ('datasets/dataset_ChannelGain_1060');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator4;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0; %!!
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = 4e-4; %logspace(-4, -2, 10 );    
            mySim.v_sigmas_toTryLF  = 67;   %linspace(50, 130, 10); 
            
            mySim.v_lambdas_toTryLB = 4e-4; %logspace(-4, -2,  10);   
            mySim.v_sigmas_toTryLB  = 25;   %linspace(20, 40, 10);
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 20;
            mySim.hybridEstimator.b_tryToBalance = 1;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            
            mySim.b_inParallel = 1;
            v_nTrains = [10 20 40 70 100 200];
            train_test_proportion = 4;
            for i_nTrain = 1:length(v_nTrains) 
                mySim.n_train = v_nTrains(i_nTrain);
                mySim.n_test  = v_nTrains(i_nTrain)*train_test_proportion;
                str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
            end
            
            tb_NMSE = struct2table(str_NMSE);
            varNames = tb_NMSE.Properties.VariableNames;
            ch_expNum = whichExp;
            save (['savedResults' filesep 'results_' ch_expNum]);

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = varNames;
            F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            F.figureNumber = str2num(ch_expNum);
            
            keyboard
        end

        function F = experiment_2070(obj, niter)
            str_dataset = load ('datasets/dataset_ChannelGain_1060');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator4;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 1;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            
            mySim.b_inParallel = 1;
            v_nTrains = [100 200 300 400 500]; %![10 20 40 70 100 200];
            %train_test_proportion = 4;
            n_test = 500;
            for i_nTrain = length(v_nTrains):-1:1
                mySim.n_train = v_nTrains(i_nTrain);
                %mySim.n_test  = v_nTrains(i_nTrain)*train_test_proportion;
                mySim.n_test = n_test;
                str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
            end
            
            tb_NMSE = struct2table(str_NMSE);
            varNames = tb_NMSE.Properties.VariableNames;
                        
            ch_expNum = whichExp;
            save (['savedResults' filesep 'results_' ch_expNum]);

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = varNames;
            F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            F.figureNumber = str2num(ch_expNum);
            
        end
        
        function F = experiment_2080(obj, niter)
            % Same as 2070, but with dataset generated in 1080
            str_dataset = load ('datasets/dataset_ChannelGain_1080');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator4;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 1;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            
            mySim.b_inParallel = 1;
            v_nTrains = [200 500];
            %train_test_proportion = 4;
            n_test = 1000;
            for i_nTrain = length(v_nTrains):-1:1
                mySim.n_train = v_nTrains(i_nTrain);
                %mySim.n_test  = v_nTrains(i_nTrain)*train_test_proportion;
                mySim.n_test = n_test;
                str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
            end
            
            tb_NMSE = struct2table(str_NMSE);
            varNames = tb_NMSE.Properties.VariableNames;
            ch_expNum = whichExp;

            save (['savedResults' filesep 'results_' whichExp]);

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = varNames;
            F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            F.figureNumber = str2num(ch_expNum);
            
        end
        
        function F = experiment_2090(obj, niter)
            % Same as 2070, but with dataset generated in 1080
            str_dataset = load ('datasets/dataset_ChannelGain_1080');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator4;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            
            mySim.b_inParallel = 1;
            v_nTrains = 2000:-300:200;
            %train_test_proportion = 4;
            n_test = 1000;
            for i_nTrain = length(v_nTrains):-1:1
                mySim.n_train = v_nTrains(i_nTrain);
                %mySim.n_test  = v_nTrains(i_nTrain)*train_test_proportion;
                mySim.n_test = n_test;
                str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
            end
            
            tb_NMSE = struct2table(str_NMSE);
            varNames = tb_NMSE.Properties.VariableNames;
            ch_expNum = whichExp;

            save (['savedResults' filesep 'results_' ch_expNum]);

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = varNames;
            F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            F.figureNumber = str2num(ch_expNum);

            
        end

        function F = experiment_2098(obj, niter)
            %just a toy-size version of 2099
            
            str_dataset = load ('datasets/dataset_ChannelGain_1080');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator4;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            
            mySim.b_inParallel = 1;
            mySim.n_monteCarloRuns = 2; %!30
            v_nTrains = [300 200]; %!2000:-300:200;
            %train_test_proportion = 4;
            n_test = 1000;
            for i_nTrain = length(v_nTrains):-1:1
                mySim.n_train = v_nTrains(i_nTrain);
                %mySim.n_test  = v_nTrains(i_nTrain)*train_test_proportion;
                mySim.n_test = n_test;
                str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
            end
            
            tb_NMSE = struct2table(str_NMSE);
            varNames = tb_NMSE.Properties.VariableNames;
            ch_expNum = whichExp;

            save (['savedResults' filesep 'results_' ch_expNum]);

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = varNames;
            F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            F.figureNumber = str2num(ch_expNum);

            
        end
        
        function F = experiment_2099(obj, niter)
            
            str_dataset = load ('datasets/dataset_ChannelGain_1099');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator4;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            
            mySim.b_inParallel = 1;
            mySim.n_monteCarloRuns = 30;
            v_nTrains = 2000:-300:200;
            %train_test_proportion = 4;
            n_test = 1000;
            for i_nTrain = length(v_nTrains):-1:1
                mySim.n_train = v_nTrains(i_nTrain);
                %mySim.n_test  = v_nTrains(i_nTrain)*train_test_proportion;
                mySim.n_test = n_test;
                str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
            end
            
            tb_NMSE = struct2table(str_NMSE);
            varNames = tb_NMSE.Properties.VariableNames;
            ch_expNum = whichExp;

            save (['savedResults' filesep 'results_' ch_expNum]);

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = varNames;
            F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            F.figureNumber = str2num(ch_expNum);

            
        end
        
        function F = experiment_2095(obj, niter)
            str_dataset = load ('datasets/dataset_ChannelGain_1095');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator3;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            
            mySim.b_cv_inParallel = 1;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            
%             mySim.b_inParallel = 1;
%             v_nTrains = 2000:-300:200;
%             %train_test_proportion = 4;
%             n_test = 1000;
%             for i_nTrain = length(v_nTrains):-1:1
%                 mySim.n_train = v_nTrains(i_nTrain);
%                 %mySim.n_test  = v_nTrains(i_nTrain)*train_test_proportion;
%                 mySim.n_test = n_test;
%                 str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
%             end
%             
%             tb_NMSE = struct2table(str_NMSE);
%             varNames = tb_NMSE.Properties.VariableNames;

            v_trainPairs = vec(str_dataset.v_indicesTrain);
            v_mapPairs   = str_dataset.v_indicesMap(:);
            
            str_dataset_sim = struct;
            s_fieldNames = fieldnames(str_dataset);
            for i_name = 1:length(s_fieldNames)
                ch_fName = char(s_fieldNames(i_name));
                if ch_fName(1)=='a'
                    str_dataset_sim.(ch_fName(2:end)) = str_dataset.(ch_fName)(:,:,1);
                else
                    str_dataset_sim.(ch_fName) = str_dataset.(ch_fName);
                end
            end
                            
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulate(...
                str_dataset_sim, v_trainPairs, v_mapPairs);
            
            ms_titles = ["locFree", "locBased", "true";
                 "locFreeH", "locBasedH", "hybrid"];
             M = str_dataset.animator.createGivenMatrix([...
                 str_mapEstimates.locFree'...
                 str_mapEstimates.locFree_fromHybrid'...
                 str_mapEstimates.locBased' ...
                 str_mapEstimates.locBased_fromHybrid' ...
                 trueGains(:), ...
                 str_mapEstimates.hybrid(:) ], ...
                 ms_titles);

            ch_expNum = whichExp;
            save (['savedResults' filesep 'results_' ch_expNum]);
            
            F = [];
        end

        function F = experiment_2100(obj, niter)
            str_dataset = load ('datasets/dataset_ChannelGain_1100');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator3;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            
            mySim.b_cv_inParallel = 1;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            
%             mySim.b_inParallel = 1;
%             v_nTrains = 2000:-300:200;
%             %train_test_proportion = 4;
%             n_test = 1000;
%             for i_nTrain = length(v_nTrains):-1:1
%                 mySim.n_train = v_nTrains(i_nTrain);
%                 %mySim.n_test  = v_nTrains(i_nTrain)*train_test_proportion;
%                 mySim.n_test = n_test;
%                 str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
%             end
%             
%             tb_NMSE = struct2table(str_NMSE);
%             varNames = tb_NMSE.Properties.VariableNames;

            v_trainPairs = vec(str_dataset.v_indicesTrain);
            v_mapPairs   = str_dataset.v_indicesMap(:);
            
            str_dataset_sim = struct;
            s_fieldNames = fieldnames(str_dataset);
            for i_name = 1:length(s_fieldNames)
                ch_fName = char(s_fieldNames(i_name));
                if ch_fName(1)=='a'
                    str_dataset_sim.(ch_fName(2:end)) = str_dataset.(ch_fName)(:,:,1);
                else
                    str_dataset_sim.(ch_fName) = str_dataset.(ch_fName);
                end
            end
                            
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulate(...
                str_dataset_sim, v_trainPairs, v_mapPairs);
            
            ms_titles = ["locFree", "locBased", "true";
                 "locFreeH", "locBasedH", "hybrid"];
             M = str_dataset.animator.createGivenMatrix([...
                 str_mapEstimates.locFree'...
                 str_mapEstimates.locFree_fromHybrid'...
                 str_mapEstimates.locBased' ...
                 str_mapEstimates.locBased_fromHybrid' ...
                 trueGains(:), ...
                 str_mapEstimates.hybrid(:) ], ...
                 ms_titles);

            ch_expNum = whichExp;

            save (['savedResults' filesep 'results_' ch_expNum]);
            
            F = [];
        end

        function F = experiment_2110(obj, niter)

            str_dataset = load ('datasets/dataset_ChannelGain_1110');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator4;
            mySim.n_monteCarloRuns = 30;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            mySim.b_cv_inParallel = 1;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            
            mySim.b_inParallel = 0;
            v_nTrains = 2000:-300:200;
            %train_test_proportion = 4;
            n_test = 1000;
            for i_nTrain = length(v_nTrains):-1:1
                mySim.n_train = v_nTrains(i_nTrain);
                disp ("starting simulation with nTrain = " + ...
                    string(mySim.n_train))
                %mySim.n_test  = v_nTrains(i_nTrain)*train_test_proportion;
                mySim.n_test = n_test;
                str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
                save results_2110_partial
            end
            
            tb_NMSE = struct2table(str_NMSE);
            varNames = tb_NMSE.Properties.VariableNames;

            ch_expNum = whichExp;

            save (['savedResults' filesep 'results_' ch_expNum]);

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = varNames;
            F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            F.figureNumber = str2num(ch_expNum);
        end
        
        function F = experiment_2510(obj, niter)
            
            str_dataset = load ('datasets/dataset_ChannelGain_1510');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator4;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.locBasedEstimator = LocationBasedEstimator;            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
            
            mySim.b_inParallel = 1;
            mySim.n_monteCarloRuns = 30;
            v_nTrains = 2000:-300:200;
            n_test = 1000;
            for i_nTrain = length(v_nTrains):-1:1
                mySim.n_train = v_nTrains(i_nTrain);
                mySim.n_test = n_test;
                str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
            end
            
            tb_NMSE = struct2table(str_NMSE);
            varNames = tb_NMSE.Properties.VariableNames;
            ch_expNum = whichExp;

            save (['savedResults' filesep 'results_' ch_expNum]);

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = {'locF (step 1)', 'locB (step 2)', 'MoE (step 8)', ...
                'MoE.locF ($\mathbf{f}_p$)', 'MoE.locB ($\mathbf{f}_l$)'};
            F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            F.ch_legendPosition = 'NorthWest';
            F.ch_interpreter = 'Latex';
            F.figureNumber = str2num(ch_expNum);
            
        end

        function F = experiment_2515(obj, niter) %produces storyboard
            
            str_dataset = load ('datasets/dataset_ChannelGain_1515');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator4;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.locBasedEstimator = LocationBasedEstimator;            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
            
            mySim.b_cv_inParallel = 1;
            mySim.b_inParallel = 0;
            
            v_trainPairs = vec(str_dataset.v_indicesTrain);
            v_mapPairs   = str_dataset.v_indicesMap(:);
            
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulateOne(...
                str_dataset, v_trainPairs, v_mapPairs);
            ch_expNum = whichExp;
            save (['savedResults' filesep 'results_' ch_expNum]);
           
            figure;
            vs_titles = ["locFree", "locBased", ...
                 "locFreeH", "locBasedH", "hybrid", "true"];
            str_dataset.animator.storyBoard([...
                 str_mapEstimates.locFree'...
                 str_mapEstimates.locBased' ...
                 str_mapEstimates.locFree_fromHybrid' ...
                 str_mapEstimates.locBased_fromHybrid' ...
                 str_mapEstimates.hybrid(:) ...
                 trueGains(:)], ...
                 vec(vs_titles), 3:3:9);
            
            F = GFigure.captureCurrentFigure();


            
%             mySim.n_monteCarloRuns = 30;
%             v_nTrains = 2000:-300:200;
%             n_test = 1000;
%             for i_nTrain = length(v_nTrains):-1:1
%                 mySim.n_train = v_nTrains(i_nTrain);
%                 mySim.n_test = n_test;
%                 str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
%             end
%             
%             tb_NMSE = struct2table(str_NMSE);
%             varNames = tb_NMSE.Properties.VariableNames;
%             ch_expNum = whichExp;
% 
%             save (['savedResults' filesep 'results_' ch_expNum]);
% 
%             F = GFigure;
%             F.m_X = v_nTrains;
%             F.m_Y = tb_NMSE.Variables';
%             F.ch_interpreter = 'none';
%             F.c_legend = varNames;
%             F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
%             F.ch_xlabel = 'Number of training samples';
%             F.ch_ylabel = 'NMSE';
%             F.figureNumber = str2num(ch_expNum);
            
        end

        function F = experiment_2520(obj, niter) %also plots corrected MoE.LocB
            
            str_dataset = load ('datasets/dataset_ChannelGain_1510');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator4;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.locBasedEstimator = LocationBasedEstimator;            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
            
            mySim.b_cv_inParallel = 1;
            mySim.b_inParallel = 1;
            
            mySim.n_monteCarloRuns = 15;
            v_nTrains = 1700:-300:800;
            n_test = 1000;
            for i_nTrain = length(v_nTrains):-1:1
                mySim.n_train = v_nTrains(i_nTrain);
                mySim.n_test = n_test;
                str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
            end
            
            tb_NMSE = struct2table(str_NMSE);
            varNames = tb_NMSE.Properties.VariableNames;
            ch_expNum = whichExp;

            save (['savedResults' filesep 'results_' ch_expNum]);

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = {'locF (step 1)', 'locB (step 2)', 'MoE (step 8)', ...
                'MoE.locF ($\mathbf{f}_p$)', 'MoE.locB ($\mathbf{f}_l$)', ...
                'MoE.locB CORRECTED'};
            F.c_styles = {'-v', '-s', '-h', '--v', '--x', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            F.ch_legendPosition = 'NorthWest';
            F.ch_interpreter = 'Latex';
            F.figureNumber = str2num(ch_expNum);
            
        end

        function F = experiment_2521(obj, niter) %also plots corrected MoE.LocB
            
            str_dataset = load ('datasets/dataset_ChannelGain_1510');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator4;
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.locBasedEstimator = LocationBasedEstimator;            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                                        
            mySim.hybridEstimator = HybridEstimator2;
            mySim.hybridEstimator.b_debugPlots = 0;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
            
            mySim.b_cv_inParallel = 1;
            mySim.b_inParallel = 1;
            
            mySim.n_monteCarloRuns = 15;
            v_nTrains = niter;
            n_test = 1000;
            for i_nTrain = length(v_nTrains):-1:1
                mySim.n_train = v_nTrains(i_nTrain);
                mySim.n_test = n_test;
                str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
            end
            
            tb_NMSE = struct2table(str_NMSE);
            varNames = tb_NMSE.Properties.VariableNames;
            ch_expNum = whichExp;

            save (['savedResults' filesep 'results_' ch_expNum]);

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = {'locF (step 1)', 'locB (step 2)', 'MoE (step 8)', ...
                'MoE.locF ($\mathbf{f}_p$)', 'MoE.locB ($\mathbf{f}_l$)', 'MoE.locB CORRECTED'};
            F.c_styles = {'-v', '-s', '-h', '--v', '--x', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            F.ch_legendPosition = 'NorthWest';
            F.ch_interpreter = 'Latex';
            F.figureNumber = str2num(ch_expNum);
            
        end

        %% Methods that show the estimated locations and the discrepancy 
        %  in three imagesc plots
        function F = experiment_3119(obj, niter)
            load datasets/dataset_ChannelGain_1119.mat ...
                m_grid_x m_grid_y am_estimatedLocations m_locations datasetGen
            i_r = 1;
            figure(1999);
            subplot(1, 3,1);
                imagesc(m_grid_x(:,1), m_grid_y(1,:), reshape(am_estimatedLocations(:,1,i_r), size(m_grid_x))')
                xlabel x
                ylabel y
                axis xy
                title 'Estimated x'
            subplot(1, 3, 2);
                imagesc(m_grid_x(:,1), m_grid_y(1,:), reshape(am_estimatedLocations(:,2,i_r), size(m_grid_x))')
                xlabel x
                ylabel y
                axis xy
                title 'Estimated y'
                
            subplot(1,3,3);
                imagesc(m_grid_x(:,1), m_grid_y(1,:), reshape(mean(vecnorm(m_locations - am_estimatedLocations(:,:,:), 2, 2), 3), size(m_grid_x))')
                xlabel x
                ylabel y
                axis xy
                title 'Norm of estimation error vector'
            
                hold on
            
                datasetGen.generator.plot_environment()
                
            F = GFigure.captureCurrentFigure();
        end
        
        function F = experiment_3129(obj, niter)
            load datasets/dataset_ChannelGain_1129.mat ...
                m_grid_x m_grid_y am_estimatedLocations m_locations datasetGen
            i_r = 1;
            figure(1999);
            subplot(1, 3,1);
                imagesc(m_grid_x(:,1), m_grid_y(1,:), reshape(am_estimatedLocations(:,1,i_r), size(m_grid_x))')
                xlabel x
                ylabel y
                axis xy
                title 'Estimated x'
            subplot(1, 3, 2);
                imagesc(m_grid_x(:,1), m_grid_y(1,:), reshape(am_estimatedLocations(:,2,i_r), size(m_grid_x))')
                xlabel x
                ylabel y
                axis xy
                title 'Estimated y'
                
            subplot(1,3,3);
                imagesc(m_grid_x(:,1), m_grid_y(1,:), reshape(mean(vecnorm(m_locations - am_estimatedLocations(:,:,:), 2, 2), 3), size(m_grid_x))')
                xlabel x
                ylabel y
                axis xy
                title 'Norm of estimation error vector'
            
                hold on
            
                datasetGen.generator.plot_environment()
                
            F = GFigure.captureCurrentFigure();
        end

        %% Experiments with HybridEstimator3
        % (the one where we attempt to estimate the gating function
        % with tanh basis functions)
        function F = experiment_4010(obj, niter)

            str_dataset = load ('datasets/dataset_ChannelGain_1110');
            % m_source_loc = str_dataset.datasetGen.locEstimator.Xenb;
            
            mySim = Simulator4;
            mySim.n_monteCarloRuns = 1; %!
            mySim.b_augmentTraining = 1;
            
            mySim.b_trainLocFree = 1;
            mySim.b_trainLocBased = 1;
            mySim.b_trainHybrid = 1;
            mySim.b_cvLambdas_hybrid = 0;
            mySim.b_cv_inParallel = 1;
            
            mySim.locFreeEstimator = LocationFreeEstimator;
%             mySim.locFreeEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLF^2));
%             mySim.locFreeEstimator.regularizationParameter = lambdaLF;
            mySim.locFreeEstimator.enableMissingData = 0;
            mySim.v_lambdas_toTryLF = logspace(-4, -2, 10 );   % 4e-4; % 
            mySim.v_sigmas_toTryLF  = linspace(50, 130, 10);   % 67;   %
            
            mySim.v_lambdas_toTryLB = logspace(-4, -2,  10);   % 4e-4; %   
            mySim.v_sigmas_toTryLB  = linspace(20, 40, 10);    % 25;   %
                
            mySim.locBasedEstimator = LocationBasedEstimator;
            %mySim.locEstimator = WangLocationEstimator;
%             mySim.locBasedEstimator.kernel = @(x, y) ...
%                 exp(-norms(x-y, 2, 1).^2/(kernelSigmaLB^2));
%             mySim.locBasedEstimator.regularizationParameter = lambdaLB;
            %mySim.locBasedEstimator.Xenb = m_source_loc;
                        
            mySim.hybridEstimator = HybridEstimator3; % interpolates the gating function with a monotone kernel
            mySim.hybridEstimator.b_debugPlots = 1;
            mySim.hybridEstimator.max_itrs_alternating = 40;
            mySim.hybridEstimator.b_tryToBalance = 1;
%             mySim.hybridEstimator.h_kernelLF = mySim.locFreeEstimator.kernel;
%             mySim.hybridEstimator.h_kernelLB = mySim.locBasedEstimator.kernel;
%             mySim.hybridEstimator.regularizationParameterLF =lambdaLF;
%             mySim.hybridEstimator.regularizationParameterLB =lambdaLB;
            mySim.hybridEstimator.method_fit_D = 'kernel';
            nu = 5;
            %mySim.hybridEstimator.kernelPar_G = nu: 
            mySim.hybridEstimator.h_kernelG = @(x, y) ...
                (1-tanh((x(1,:)-y(1,:))/nu)).*(1-tanh((x(2,:)-y(2,:))/nu))/4;
            mySim.hybridEstimator.regularizationParameterG = 1e-3;

%             n_allPairs = size(str_dataset.m_pairs, 1);
%             v_trainTestPairs = randperm(n_allPairs, 120)';
%             v_trainPairs = v_trainTestPairs(1:100);
%             v_testPairs  = v_trainTestPairs(101:120);
            
            mySim.b_inParallel = 0;
            v_nTrains = [30 60]; %! 2000:-300:200;
            %train_test_proportion = 4;
            n_test = 100; %! 1000
            for i_nTrain = length(v_nTrains):-1:1
                mySim.n_train = v_nTrains(i_nTrain);
                disp ("starting simulation with nTrain = " + ...
                    string(mySim.n_train))
                %mySim.n_test  = v_nTrains(i_nTrain)*train_test_proportion;
                mySim.n_test = n_test;
                str_NMSE(i_nTrain) = mySim.simulateMonteCarlo(str_dataset);
                save results_2110_partial
            end
            
            tb_NMSE = struct2table(str_NMSE);
            varNames = tb_NMSE.Properties.VariableNames;

            save results_4010

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = varNames;
            F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            F.figureNumber = 4010;
        end

    end
end

            
            
            


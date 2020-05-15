classdef ChannelGainExperiments < ExperimentFunctionSet
    
    properties
        
        b_loadData = 0 % set this property to 1 if you want the experiment
        % methods to load data from the savedResults folder.
        
        % You may define here constants that you will use in the
        % experiments of this file
        
    end
    
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
            
            filename = 'dataset_ChannelGain_1010';
            
            %%%
            save (['datasets' filesep filename], '-struct', 'str_dataset');
            
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
            
            filename = 'dataset_ChannelGain_1030';
            
            %%%
            save (['datasets' filesep filename], '-struct', 'str_dataset');
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
            
            filename = 'dataset_ChannelGain_1060';
            
            %%%
            save (['datasets' filesep filename], '-struct', 'str_dataset');
%             
%             figure(999); 
%             my_datasetGen.generator.plot_environment;
%             F = GFigure.captureCurrentFigure;

        end
        
        function F = experiment_1080(obj, niter)
            % This experiment function creates a dataset to be used in
            % experiment_2060, where we want to average over a number of
            % Monte Carlo realizations
                        
            rng(1)
            
            s_fileName_environment = 'modelFiles/env2.mat';           
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
            
            filename = 'dataset_ChannelGain_1080';
            
            %%%
            save (['datasets' filesep filename], '-struct', 'str_dataset');
%             
%             figure(999); 
%             my_datasetGen.generator.plot_environment;
%             F = GFigure.captureCurrentFigure;

        end

        %% These methods take the datasets in the previous section
        % and run simulations based on them
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
             
            save 
            
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
             
            save results_2040
            
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
             
            save results_2050
            
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
             
            save results_2055
            
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

        % Monte Carlo simulatinos
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
            varNames = t_NMSE.Properties.VariableNames;

            save results_2060

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = varNames;
            F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            
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
            
            mySim.b_inParallel = 0; %!
            v_nTrains = [50 80 120]; %![10 20 40 70 100 200];
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

            save results_2070

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = varNames;
            F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            
            keyboard
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

            save results_2080

            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = varNames;
            F.c_styles = {'-v', '-s', '-h', '--v', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            
            keyboard
        end

    end
end

            
            
            

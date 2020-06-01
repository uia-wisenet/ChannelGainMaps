classdef GlobecomExperiments < ExperimentFunctionSet
    
    methods
        %% Prepare both datasets: 30 independent realizations, and animation
        function F = experiment_1010(obj, niter)
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
            [my_datasetGen.generator, ~] = ...
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

            % define grid:
            x1 = my_datasetGen.generator.boundary(1,1);
            x2 = my_datasetGen.generator.boundary(1,2);
            y1 = my_datasetGen.generator.boundary(2,1);
            y2 = my_datasetGen.generator.boundary(2,2);
            [m_grid_x, m_grid_y] = ndgrid(...
                linspace(x1, x2, gridSize(1)), ...
                linspace(y1, y2, gridSize(2)));
            
            m_locations = cat(2, m_grid_x(:), m_grid_y(:));
            n_locations = numel(m_grid_x);
            m_allPairs  = combnk(1:n_locations, 2);
            n_allPairs = size(m_allPairs, 1);
            
            %% Generate the training set:
            n_trainSamples = 3000;
            m_pairsTrain = m_allPairs(...
                randperm(n_allPairs, n_trainSamples),:);
            
            my_datasetGen.n_realizations = 36;
            my_datasetGen.b_inParallel = 1;
            
            % the Mambo line (30 MC realizations):
            str_dataset = my_datasetGen.generate(m_locations, m_pairsTrain);

            str_dataset.datasetGen = my_datasetGen;
            str_dataset.m_grid_x   = m_grid_x;
            str_dataset.m_grid_y   = m_grid_y;
            
            save (['datasetsGlobecom' filesep 'dataset_30realizations'], ...
                '-struct', 'str_dataset');
            
            %% Set the datasetGen up for the animation:
            
            my_animator = Animator;
            my_animator.t_grid_xy = cat(3, m_grid_x, m_grid_y);
            v_rows_frame = 1:2:gridSize(1);
            my_animator.v_indicesTxPosition = sub2ind(size(m_grid_x), ...
                v_rows_frame, ...
                floor(gridSize(2)/2)*ones(1, length(v_rows_frame)));
            my_animator.generator = my_datasetGen.generator;
            m_pairsMap = my_animator.m_pairs();

            
            m_pairsGen = [m_pairsMap; m_pairsTrain];
            v_indicesMap = 1:size(m_pairsMap,1);
            v_indicesTrain = (size(m_pairsMap,1)+1):size(m_pairsGen,1);
                       
            my_datasetGen.n_realizations = 1;
            my_datasetGen.b_inParallel = 1;
            
            % the Mambo line (slices for animation/storyboard):
            str_dataAnimation = my_datasetGen.generate(m_locations, m_pairsGen);
            
            % if we want to visualize:
            % M = an.create(repmat(str_dataset.v_channelGains, [1 3]));           
            
            str_dataAnimation.datasetGen     = my_datasetGen;
            str_dataAnimation.m_grid_x       = m_grid_x;
            str_dataAnimation.m_grid_y       = m_grid_y;
            str_dataAnimation.v_indicesMap   = v_indicesMap;
            str_dataAnimation.v_indicesTrain = v_indicesTrain;
            str_dataAnimation.animator       = my_animator;
            
            save (['datasetsGlobecom' filesep 'dataset_animation'], ...
                '-struct', 'str_dataAnimation');
            
            F = [];
        end
        
        function F = experiment_2010(obj, niter) % produces comparison plot
            % figure 2
            str_dataset = load(['datasetsGlobecom' filesep ...
                'dataset_30realizations']);
            
            mySim = obj.globecomSimulator();
            
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
            
            tb_NMSE_toDraw = removevars(tb_NMSE, 'locBased_fromHybrid');
            
            F = GFigure;
            F.m_X = v_nTrains;
            F.m_Y = tb_NMSE_toDraw.Variables';
            F.ch_interpreter = 'none';
            F.c_legend = {'locF (step 1)', 'locB (step 2)', 'MoE (step 8)', ...
                'MoE.locF ($\mathbf{f}_p$)', 'MoE.locB ($\mathbf{f}_l$)'};
            F.c_styles = {'-v', '-s', '-h', '--v', '--x', '--s'};
            F.ch_xlabel = 'Number of training samples';
            F.ch_ylabel = 'NMSE';
            F.ch_legendPosition = 'NorthWest';
            F.ch_interpreter = 'Latex';
            F.figureNumber = str2num(ch_expNum);

        end
        
        function F = experiment_2020(obj, niter) % produces storyboard
            % figure 1
            str_dataset = load(['datasetsGlobecom' filesep ...
                'dataset_animation']);
            
            mySim = obj.globecomSimulator();
            
            mySim.n_monteCarloRuns = 1;
            nTrain = 2000;
            
            v_trainPairs = vec(str_dataset.v_indicesTrain(1:nTrain));
            v_mapPairs   = str_dataset.v_indicesMap(:);
            
            [str_NMSE, str_mapEstimates, trueGains] = mySim.simulateOne(...
                str_dataset, v_trainPairs, v_mapPairs);
            
            ch_expNum = whichExp;
            save (['savedResults' filesep 'results_' ch_expNum]);
            
            str_dataset.animator.ch_colormap = 'jet';
            str_dataset.animator.c_plotSource_options{1} = 'pk';
            str_dataset.animator.c_plotSource_options{3} = 'r';
            str_dataset.animator.c_plotSource_options{5} = 5;

            figure(999);
            vs_titles = ["locF", "locB", "MoE.locF", "MoE.locB", ...
                 "MoE.locB(Corrected)", "MoE", "true"];
            h_colorBar = str_dataset.animator.storyBoard([...
                 str_mapEstimates.locFree'...
                 str_mapEstimates.locBased' ...
                 str_mapEstimates.locFree_fromHybrid' ...
                 str_mapEstimates.locBased_fromHybrid' ...
                 str_mapEstimates.locBased_fromHybrid_corrected' ...
                 str_mapEstimates.hybrid(:) ...
                 trueGains(:)], ...
                 vec(vs_titles), [2 5 9]);
            h_colorBar.Position = [0.8333 0.1235 0.0119 0.7840];
            h_colorBar.Limits = [-80 -40];

            F = GFigure.captureCurrentFigure();
        end
    end
    
    methods (Static)
        %% constructs Simulator4 object to be used in exp 2010 and 2020
        function mySim = globecomSimulator() 
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
        end
    end
end

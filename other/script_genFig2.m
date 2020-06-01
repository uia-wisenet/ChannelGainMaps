clear
gsinfo;
load savedResults/results_2520.mat

F = GFigure;
F.m_X = v_nTrains;
F.m_Y = tb_NMSE.Variables';
F.ch_interpreter = 'none';
F.c_legend = {'locF (step 1)', 'locB (step 2)', 'MoE (step 8)', ...
'MoE.locF ($\mathbf{f}_p$)', 'MoE.locB ($\mathbf{f}_l$)', 'MoE.locB CORRECTED'};
F.c_styles = {'-v', '-s', '-h', '--v', '--x' '--s'};
F.ch_xlabel = 'Number of training samples';
F.ch_ylabel = 'NMSE';
F.ch_legendPosition = 'NorthWest';
F.ch_interpreter = 'Latex';
F.figureNumber = str2num(ch_expNum);

F.plot
ff = gcf;
ff.Children(2).YScale = 'log';
ylim([0.25, 1]);
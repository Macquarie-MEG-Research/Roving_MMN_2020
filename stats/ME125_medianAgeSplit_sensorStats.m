%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ME125 (Roving MMF) SENSOR-LEVEL ANALYSIS
% NON-PARAMETRIC CLUSTER-BASED RANDOM PERMUTATION TESTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Set up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('lay.mat');
load ('neighbours_125.mat')


% = PLEASE SPECIFY =

% (1) Folder locations 
% where to read MEG data from:
data_path = '..\\..\\ME125_roving_Phase1_data_37kids\\';
%data_path = '..\\..\\ME125_roving_adult_data\\';

% where to store results:
output_path = 'D:\\Judy\\RA_2020\\ARC_Roving_MMN\\Phase1_Results_young-vs-old\\'; % full path required on Windows, due to back-slash issues
%output_path = 'D:\\Judy\\RA_2020\\ARC_Roving_MMN\\Phase1_Results_adult\\'; 

% (2) The group(s) of participants to analyse: 
% (if you specify two groups, e.g. 'younger', 'older', then these two 
% groups will also be compared with each other at the end)

group_list = {'younger', 'older'};
%group_list = {'adult'};

% The following lists are set up for ME125 (roving) Phase 1 - change if analysing other studies
group.older   = {'2913' '2787' '2697' '2702' '2786' '2716' '2698' '2712' '2872' '2703' '2888' '2811' '2696' '2713' '2904' '2854' '2699' '2858'}; % 18 kids, >=5yo
group.younger = {'2724' '2642' '2866' '2785' '2793' '2738' '2766' '2687' '2629' '2897' '2683' '2695' '2739' '2810' '2632' '2667' '2875' '2912' '2681'}; % 19 kids, <5yo
folders = dir([data_path '2*']);
group.adult = vertcat({folders(:).name});

% (3) Perform baseline correction? if so, specify the baseline interval
DO_BASELINE = false;
ERF_BASELINE = [-0.1 0];

% (4) Other settings
alpha_thresh = 0.05;  % threshold for stats
x_lims       = [0 0.4];
save_to_file = 'yes'; % save figures to file?


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Global field power (separately for young & old)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(data_path)
orig = cd;

folders = dir('2*');

for i=1:length(group_list)

    idx = find(ismember({folders.name}, group.(group_list{i})));

    avg_deviant_planar_all     = [];
    avg_standard_planar_all    = [];
    avg_mmf_planar_all         = [];

    for sub=1:length(idx)
        if strcmp(group_list{i}, 'adult')
            cd([folders(idx(sub)).name '\\']);
        else
            cd([folders(idx(sub)).name '\\ReTHM\\']);
        end

        if ~DO_BASELINE % if no baseline correction needed, just load the saved ERFs (ie. ave)
            load('deviant_ave.mat')
            load('predeviant_ave.mat')
            
        else % otherwise, load the epoched data, perform baseline correction then compute ERF
            load('deviant.mat')
            load('predeviant.mat')

            % do baseline correction on individual epochs
            cfg = [];
            cfg.baseline = ERF_BASELINE;
            deviant = ft_timelockbaseline(cfg, deviant);

            cfg = [];
            cfg.baseline = ERF_BASELINE;
            predeviant = ft_timelockbaseline(cfg, predeviant);

            % Compute ERFs
            deviant_ave                  = ft_timelockanalysis([], deviant);
            %save('deviant_ave_baseline_-100_to_0ms.mat', 'deviant_ave');
            predeviant_ave               = ft_timelockanalysis([], predeviant);
            %save('predeviant_ave_baseline_-100_to_0ms.mat', 'predeviant_ave');
        end
        
        standard_ave = predeviant_ave; % rename

        
        cfg                         = [];
        cfg.method                  = 'triangulation';
        cfg.neighbours              = neighbours; %To get neigbours of the pak for plotting

        cfg.planarmethod            = 'sincos';
        deviant_ave_planar          = ft_megplanar(cfg, deviant_ave);
        standard_ave_planar         = ft_megplanar(cfg, standard_ave);

        % Combine the horizontal and vertical components of the planar gradient
        cfg = [];
        avg_deviant_planar_comb     = ft_combineplanar(cfg, deviant_ave_planar);

        cfg = [];
        avg_standard_planar_comb  = ft_combineplanar(cfg, standard_ave_planar);

        % %  computes the planar gradient magnitude over both directions
        % %  combining the two gradients at each sensor to a single positive-valued number. This
        % %  can be done for single-trial/averaged planar gradient ERFs or single-trial/averaged TFRs.
        avg_deviant_planar_all{sub} = avg_deviant_planar_comb;
        avg_standard_planar_all{sub} = avg_standard_planar_comb;

        cfg=[];
        cfg.operation='subtract';
        cfg.parameter='avg';
        avg_mmf_planar_all{sub}=ft_math(cfg,avg_deviant_planar_all{sub},avg_standard_planar_all{sub});

        cd(orig)
    end


    cfg = [];
    cfg.channel   = 'all';
    cfg.latency   = 'all';
    cfg.parameter = 'avg';
    %cfg.keepindividual = 'yes';
    GA_planar_standard    = ft_timelockgrandaverage(cfg,avg_standard_planar_all{:});
    GA_planar_Dev       = ft_timelockgrandaverage(cfg,avg_deviant_planar_all{:});
    GA_planar_MMF       = ft_timelockgrandaverage(cfg,avg_mmf_planar_all{:});
    % "{:}" means to use data from all elements of the variable

    save([output_path, 'GA_', group_list{i}], 'GA_planar_standard', 'GA_planar_Dev', 'GA_planar_MMF');

    cfg = [];
    cfg.method    = 'power';
    cfg.channel   = 'all';
    GA_planar_standard_GFP = ft_globalmeanfield(cfg,GA_planar_standard);
    GA_planar_Dev_GFP = ft_globalmeanfield(cfg,GA_planar_Dev);
    GA_planar_MMF_GFP = ft_globalmeanfield(cfg,GA_planar_MMF);

    cfg = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'avg';
    GA_difference_planar_GFP = ft_math(cfg,GA_planar_Dev_GFP,GA_planar_standard_GFP);
    GA_difference_planar_GFP.avg = abs(GA_difference_planar_GFP.avg);


    cfg                          = [];
    cfg.xlim                     = [-0.1 0.5];
    %cfg.ylim                     = [0 4e-28];
    cfg.title                    = 'Global Field Power';
    cfg.graphcolor               = 'brk';
    cfg.linestyle                = '-k';
    %cfg.graphcolor              = 'k';
    cfg.linewidth                = 6;
    %cfg.zlim                    = [-4.9 4.9];
    cfg.showlabels               = 'yes';
    cfg.fontsize                 = 6;

    figure;
    % ft_singleplotER(cfg,standard_1_GFP,standard_2_GFP,standard_3_GFP,standard_4_GFP,standard_5_GFP,deviant_GFP)
    ft_singleplotER(cfg,GA_planar_standard_GFP ,GA_planar_Dev_GFP,GA_planar_MMF_GFP); hold on;
    legend('Predeviant','Deviant','MMF')
    title('Global Field Power'); drawnow;
    xlabel('Time (sec)');
    ylabel('Global Mean Field Amplitude (Tesla/cm^{2})');
    set(gca,'fontsize', 40);

    x0=10;
    y0=10;
    width=1400;
    height=1000;
    set(gcf,'position',[x0,y0,width,height])

    print([output_path 'GFP_planar_' group_list{i}],'-dpng');

    % note to self: This looked wrong at first glance (stand_ave.grad.type?)
    % but it's actually fine. PSF: This is not a problem as the grad structure is static. I.e. the geometry of the sensors is not changing between conditions.

    % save individual-subjects ERFs
    for k = 1:length(avg_deviant_planar_all)
        avg_standard_planar_all{1,k}.grad.type = deviant_ave.grad.type;
        avg_deviant_planar_all{1,k}.grad.type = deviant_ave.grad.type;
    end

    save ([output_path,'ERF_planar_',group_list{i}], 'avg_deviant_planar_all','avg_standard_planar_all','avg_mmf_planar_all')

    
    %% Settings for NPCBRPT
    clear stat_wholeEpoch
    %     clear (['stat_wholeEpoch_',group_list{i}])
    cfg = [];
    %cfg.grad = deviant_ave.grad;
    cfg.channel     = 'all';
    cfg.neighbours  = neighbours; % defined as above
    cfg.latency     = x_lims;  % timewindow for the stats. Epoched earlier, but doesn't make sense to do the stats before 0, since an effect before the stimulus (i.e., before 0) would be meaningless/not interpretable.
    cfg.avgovertime = 'no'; %
    cfg.parameter   = 'avg';
    cfg.method      = 'montecarlo';
    cfg.statistic   = 'ft_statfun_depsamplesT';
    cfg.alpha       = alpha_thresh; % threshold for the significant clusters when doing the Mote Carlo p-value comparison step
    cfg.clusteralpha = 0.05; % threshold for initial clustering before correction
    cfg.correctm    = 'cluster';
    %cfg.correcttail = 'prob';
    cfg.correcttail = 0;
    cfg.numrandomization = 1000;
    cfg.minnbchan        = 2; % minimal neighbouring channels

    % DESIGN
    Nsub = sub;
    cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
    cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
    cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
    cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

    stat_wholeEpoch = ft_timelockstatistics(cfg,avg_deviant_planar_all{:},avg_standard_planar_all{:});

    %save (['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/variables/','stat_wholeEpoch_',group_list{i}], 'stat_wholeEpoch')
    save ([output_path,'stat_wholeEpoch_',group_list{i}], 'stat_wholeEpoch')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Plot the clusters (plot sensor-level statistics)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(group_list)

    % Load data generated from code above
    load ([output_path,'stat_wholeEpoch_',group_list{i}]);
    load ([output_path,'GA_',group_list{i}]);

    % Author: Robert Seymour (robert.seymour@mq.edu.au)
    %
    %%%%%%%%%%%
    % Inputs:
    %%%%%%%%%%%
    %
    % - stats       = output from ft_timelockstatistics
    % - alpha       = alpha-level
    % - lay         = layout generated from ft_preparelayout
    % - x_lims      = x-limit of your graph (i.e. which times do you want to
    %               plot?)
    % save_to_file  = 'yes' will result in .png file being saved to file


    % % if group_list{1}
    % %     stat         = 'stat_wholeEpoch_younger.mat';
    % % else if group_list{2}
    % %         stat         = 'stat_wholeEpoch_older.mat';
    % %     end
    % % end

    stat   = stat_wholeEpoch;

    ft_warning('This function only supports plotting POSITIVE clusters at present');

    % FIX HERE - only creates plots for positive clusters; but it does still compute\
    % variables for pos & neg clusters, just doesn't plot them all
    
    % Note: At the moment, there are no neg MMF clusters anyway.
    %       But if need to plot neg clusters, just copy the plotting code 
    %       from the "GROUP COMPARISON" section below

    % Find positive clusters under the specified alpha
    if isfield(stat,'posclusters')
        pos_cluster_pvals = [stat.posclusters(:).prob];

        % If this array is empty, return error
        if isempty(pos_cluster_pvals)
            error('NO POSITIVE CLUSTERS FOUND');
        else
            % Find the p-values of each cluster
            pos_signif_clust = find(pos_cluster_pvals < alpha_thresh);

            % Give the user some feedback in Command Window
            fprintf('Positive Clusters below %.3f alpha level: %d\n',...
                alpha_thresh,length(pos_signif_clust));

            for t = 1:length(pos_signif_clust)
                fprintf('Positive Cluster #%d: %.3f\n',pos_signif_clust(t),...
                    pos_cluster_pvals(pos_signif_clust(t)));
            end

            %% For each positive cluster...
            for t = 1:length(pos_signif_clust)

                % Get the significant channels
                pos = ismember(stat.posclusterslabelmat, pos_signif_clust(t));
                highlight_chan = any(pos(:,:)');

                %pos_cluster_mask(t,:)=highlight_chan;

                % Get the significant times
                index = (any(stat.posclusterslabelmat == pos_signif_clust(t)));
                time_for_topo = stat.time(index');

                stat.index = pos;

                % Find the time of the peak
                cfg = [];
                cfg.latency = [time_for_topo(1) time_for_topo(end)];
                cfg.channel = stat.label(highlight_chan');
                data_for_peak = ft_selectdata(cfg,stat);

                avg_chan  = mean(data_for_peak.stat(:,:));
                time_of_peak = data_for_peak.time(find(max(avg_chan)==avg_chan));

                %% Singleplot: t-values of sig sensors
                cfg                 = [];
                cfg.channel         = stat.label(highlight_chan');
                cfg.maskparameter   = 'index'; % automatically add shaded region
                cfg.xlim            = [0 0.4];
                cfg.linestyle       = '-k';
                cfg.graphcolor      = 'k';
                cfg.linewidth       = 6;
                %cfg.zlim            = [-4.9 4.9];
                cfg.showlabels      = 'yes';
                cfg.fontsize        = 6;
                cfg.layout          = lay;
                cfg.parameter       = 'stat';
                figure;
                ft_singleplotER(cfg,stat); hold on;
                scatter(time_of_peak,max(avg_chan),40,'filled','r');

                % Give the title
                title(sprintf('Cluster: #%d\n Time:  %.3fs to %.3fs\nPeak: %.3fs' ...
                    ,t,time_for_topo(1),...
                    time_for_topo(end),time_of_peak));

                set(gca,'fontsize', 20);
                set(gca, 'Layer','top');

                xlabel('Time (sec)','FontSize',24);
                ylabel('t-value','FontSize',24);

                % Save as png
                if strcmp(save_to_file,'yes')
                    disp('Saving figure to .png file');
                    print(sprintf([output_path, group_list{i}, '_tvalues_pos_cluster_%d'],t),'-dpng','-r200');
                else
                    disp('Not saving figure to file');
                end

                %% Singleplot: average amplitude of sig sensors               
                figure('Name', 'Average of significant channels');
                
                cfg       = [];
                cfg.title = ' '; % hide the display of channel names at the top
                cfg.channel                  = stat.label(highlight_chan');
                cfg.baseline                 = ERF_BASELINE; % makes no diff if we've already done baseline correction earlier
                cfg.xlim                     = [-0.1 0.4];
                cfg.graphcolor               = 'brk';
                cfg.linestyle                = '-k';
                cfg.linewidth                = 3;
                cfg.showlabels               = 'yes';
                cfg.fontsize                 = 6;
                ft_singleplotER(cfg, GA_planar_standard, GA_planar_Dev, GA_planar_MMF); hold on;
                legend('Predeviant','Deviant','MMF')
    
                % Give the title
                title(sprintf('Cluster: #%d\n Time:  %.3fs to %.3fs\nPeak: %.3fs' ...
                    ,t,time_for_topo(1),...
                    time_for_topo(end),time_of_peak));
                
                x0=10;
                y0=10;
                width=800;
                height=550;
                set(gcf,'position',[x0,y0,width,height]) % specify the size of the figure (on screen)

                set(gca,'fontsize', 20);
                set(gca, 'Layer','top');

                xticks([-0.1:0.1:0.4]); % force it to show all tick values
                xlabel('Time (sec)','FontSize',24);
                ylabel('Tesla','FontSize',24);
                box on; % draw a border around the figure

                % create shaded region indicating effect duration
                start_time = time_for_topo(1);
                end_time = time_for_topo(end);
                ylimits = ylim; ylow = ylimits(1); yhigh = ylimits(2);
                x = [start_time end_time end_time start_time]; % specify x,y coordinates of the 4 corners
                y = [ylow ylow yhigh yhigh];
                patch(x,y,'black', 'FaceAlpha',0.3, 'HandleVisibility','off') % draw the shade 
                    % (turn off HandleVisibility so it won't show up in the legends)
                ylim(ylimits); % ensure ylim doesn't get expanded

                % Save as png
                if strcmp(save_to_file,'yes')
                    disp('Saving figure to .png file');
                    print(sprintf([output_path, group_list{i}, '_avg_amplitude_pos_cluster_%d'],t),'-dpng','-r200');
                else
                    disp('Not saving figure to file');
                end
                
                
                %% Topoplot
                cfg                  = [];
                cfg.interpolation    = 'v4';
                cfg.marker           = 'off';
                cfg.highlight        = 'on';
                cfg.highlightchannel = stat.label(highlight_chan);
                cfg.highlightsymbol  = '.';
                cfg.highlightsize   = 20;
                cfg.xlim            = [time_for_topo(1) time_for_topo(end)];
                %   cfg.zlim            = 'maxabs';
                cfg.zlim            = [-5 5];
                cfg.comment         = 'no';
                cfg.fontsize        = 6;
                cfg.layout          = lay;
                cfg.parameter       = 'stat';
                figure;colorbar
                ft_topoplotER(cfg,stat); hold on;
                ft_hastoolbox('brewermap', 1);
                colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

                % Give the title
                title(sprintf('Cluster: #%d\n Time: %.3fs to %.3fs\nPeak: %.3fs' ,t,time_for_topo(1),...
                    time_for_topo(end),time_of_peak));

                set(gca,'fontsize', 20);

                % Save as png
                if strcmp(save_to_file,'yes')
                    disp('Saving figure to .png file');
                    %print(sprintf(['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/png/',...
                    %    group_list{i},'_topoplot_pos_cluster_%d'],t),'-dpng','-r200');
                    print(sprintf([output_path,group_list{i},'_topoplot_pos_cluster_%d'],t),'-dpng','-r200');
                else
                    disp('Not saving figure to file');
                end

            end
        end
    end
end


%% printing stats to the command window

% check the stats: finds all the clusters from initial phase, and checks sig. when compared to
% the null dist.
% run these lines separately

%{
load([output_path 'stat_wholeEpoch_older.mat'])
old_stat=stat_wholeEpoch;
old_stat.posclusters.prob % where pos would = D > S
old_stat.negclusters.prob % where pos would = S > D

load([output_path 'stat_wholeEpoch_younger.mat'])
young_stat=stat_wholeEpoch;
young_stat.posclusters.prob % where pos would = D > S. Needs to be < 0.05 to be sig.
young_stat.negclusters.prob % where pos would = S > D

% if error, probably = no cluster
%}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. GROUP COMPARISON - MMF in younger vs older kids (cluster-based t-test)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ONLY RUN THIS SECTION IF there are two groups specified at the top
if length(group_list) ~= 2
    error('Note: Only 1 group of participants were specified for analysis. Not performing any group comparison.'); % this will terminate the script
end


cd(output_path)

load('ERF_planar_older.mat')
old_mmf = avg_mmf_planar_all;

load('ERF_planar_younger.mat')
young_mmf = avg_mmf_planar_all;


cfg                  = [];
cfg.channel          = 'all';
cfg.neighbours       = neighbours; % defined as above
cfg.latency          = x_lims;
cfg.avgovertime      = 'no'; %
cfg.parameter        = 'avg';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05; % threshold for initial clustering before correction
cfg.method           = 'montecarlo';
cfg.alpha            = alpha_thresh; % threshold for the significant clusters when doing the Mote Carlo p-value comparison step
cfg.correcttail      = 0;
cfg.numrandomization = 1000;
cfg.minnbchan        = 2; % minimal neighbouring channels

% Design Matrix
design = [ones(1,length(group.older)) 2*ones(1,length(group.younger))];
cfg.design      = design; 
cfg.ivar        = 1; % row of design matrix that contains independent variable (the conditions)
%%% EXPLANATION: 'ones' is used to create a matrix. e.g.1., 'ones(4)' creates
%%% a 4x4 matrix of 1s. e.g.2., 'ones(1,4)' creates one row of four 1s. For
%%% b/w-subjects, comparing a row of group A/old (1s) with group B/young (2s).
%%% To create the 2s = 2*ones = 2s. Call this 'deisgn' into the cfg settings.
%%% the ivar is 'group' (with two levels, young vs old


stat_MMFbyGroup           = ft_timelockstatistics(cfg,old_mmf{:},young_mmf{:});
save ([output_path,'stat_MMFbyGroup'], 'stat_MMFbyGroup');


% command window will print # pos and neg clusters BEFORE the correction
% to find out how many remain after correction (run these lines
% individually):
%{
stat_MMFbyGroup.posclusters.prob
stat_MMFbyGroup.negclusters.prob % error probably cuz no neg clusters
%}

%% plot stats
%{
stat=stat_MMFbyGroup;

% Find the clusters under the specified alpha
if isfield(stat,'posclusters')
    pos_cluster_pvals = [stat.posclusters(:).prob];

    % If this array is empty, return error
    if isempty(pos_cluster_pvals)
        error('NO POSITIVE CLUSTERS FOUND');
    else
       % Find the p-values of each cluster
        pos_signif_clust = find(pos_cluster_pvals < alpha);

        % Give the user some feedback in Command Window
        fprintf('Positive Clusters below %.3f alpha level: %d\n',...
            alpha,length(pos_signif_clust));

        for t = 1:length(pos_signif_clust)
            fprintf('Positive Cluster #%d: %.3f\n',pos_signif_clust(t),...
                pos_cluster_pvals(pos_signif_clust(t)));
        end

        %% For each positive cluster...
        for t = 1:length(pos_signif_clust)

            % Get the significant channels
            pos = ismember(stat.posclusterslabelmat, pos_signif_clust(t));
            highlight_chan = any(pos(:,:)');

            %pos_cluster_mask(t,:)=highlight_chan;

            % Get the significant times
            index = (any(stat.posclusterslabelmat == pos_signif_clust(t)));
            time_for_topo = stat.time(index');

            stat.index = pos;

            % Find the time of the peak
            cfg = [];
            cfg.latency = [time_for_topo(1) time_for_topo(end)];
            cfg.channel = stat.label(highlight_chan');
            data_for_peak = ft_selectdata(cfg,stat);

            avg_chan  = mean(data_for_peak.stat(:,:));
            time_of_peak = data_for_peak.time(find(max(avg_chan)==avg_chan));

            %% Singleplot
            cfg                 = [];
            cfg.channel         = stat.label(highlight_chan');
            cfg.maskparameter   = 'index';
            cfg.xlim            = x_lims;
            cfg.linestyle       = '-k';
            cfg.graphcolor      = 'k';
            cfg.linewidth       = 6;
            cfg.zlim            = [-3 3];
            cfg.showlabels      = 'yes';
            cfg.fontsize        = 6;
            cfg.layout          = lay;
            cfg.parameter       = 'stat';
            figure;
            ft_singleplotER(cfg,stat); hold on;
            scatter(time_of_peak,max(avg_chan),40,'filled','r');

            % Give the title
            title(sprintf('Cluster: #%d\n Time:  %.3fs to %.3fs\nPeak: %.3fs' ...
                ,t,time_for_topo(1),...
                time_for_topo(end),time_of_peak));

            set(gca,'fontsize', 20);
            set(gca, 'Layer','top');


            xlabel('Time (sec)','FontSize',24);
            ylabel('t-value','FontSize',24);

            % Save as png
            if strcmp(save_to_file,'yes')
                disp('Saving figure to .png file');
                %print(sprintf(['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/png/',...
                %    group_list{i},'_singleplot_pos_cluster_%d'],t),'-dpng','-r200');
                print(sprintf([output_path,group_list{i},'_singleplot_pos_cluster_%d'],t),'-dpng','-r200');

            else
                disp('Not saving figure to file');
            end


            %% Topoplot
            cfg                  = [];
            cfg.interpolation    = 'v4';
            cfg.marker           = 'off';
            cfg.highlight        = 'on';
            cfg.highlightchannel = stat.label(highlight_chan);
            cfg.highlightsymbol  = '.';
            cfg.highlightsize   = 20;
            cfg.xlim            = [time_for_topo(1) time_for_topo(end)];
            %   cfg.zlim            = 'maxabs';
            cfg.zlim            = [-2 2];
            cfg.comment         = 'no';
            cfg.fontsize        = 6;
            cfg.layout          = lay;
            cfg.parameter       = 'stat';
            figure;colorbar
            ft_topoplotER(cfg,stat); hold on;
            ft_hastoolbox('brewermap', 1);
            colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

            % Give the title
            title(sprintf('Cluster: #%d\n Time: %.3fs to %.3fs\nPeak: %.3fs' ,t,time_for_topo(1),...
                time_for_topo(end),time_of_peak));

            set(gca,'fontsize', 20);

%             % Save as png
%             if strcmp(save_to_file,'yes')
%                 disp('Saving figure to .png file');
%                 %print(sprintf(['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/png/',...
%                     %group_list{i},'_topoplot_pos_cluster_%d'],t),'-dpng','-r200');
%                 print(sprintf([output_path,group_list{i},'_topoplot_pos_cluster_%d'],t),'-dpng','-r200');
%
%             else
%                 disp('Not saving figure to file');
%             end

% NB. The topoplot here doesn't show where the effect is, rather the
% difference in the sensors between young v old MMF
% need to look at old & young topos separately

        end
    end
end
%}

%% Plot MMF old v young (with sig mask)

% NB. Not plotting the GFP, only plotting the channels (averaged) where the effect was

% grand mean/a for old MMF & young MMF
GA_old_mmf   = ft_timelockgrandaverage([],old_mmf{:});
GA_young_mmf = ft_timelockgrandaverage([],young_mmf{:});

% mask_param = stat.mask;


stat=stat_MMFbyGroup;

% find neg clusters that are sig
if isfield(stat,'negclusters')
    neg_cluster_pvals = [stat.negclusters(:).prob];

    % If this array is empty, return error
    if isempty(neg_cluster_pvals)
        error('NO NEGATIVE CLUSTERS FOUND');
    else
       % Find the p-values of each cluster
        neg_signif_clust = find(neg_cluster_pvals < alpha_thresh);
        
        % Give the user some feedback in Command Window
        fprintf('There are %d negative Clusters below %.3f alpha level.\n',...
            alpha_thresh,length(neg_signif_clust));

        for t = 1:length(neg_signif_clust)
            fprintf('Negative Cluster #%d: p = %.3f\n', neg_signif_clust(t),...
                neg_cluster_pvals(neg_signif_clust(t)));
        end

        % For each cluster...
        for t = 1:length(neg_signif_clust)
            % Get the significant channels
            pos = ismember(stat.negclusterslabelmat, neg_signif_clust(t));
            highlight_chan = any(pos(:,:)');
            
            % Get the significant times
            index = (any(stat.negclusterslabelmat == neg_signif_clust(t)));
            time_for_topo = stat.time(index');
                
            % Plot
            [x,y] = find(pos);
            figure; 
            plot(GA_old_mmf.time,mean(GA_old_mmf.avg(sort(unique(x)),:)),'g','LineWidth',5); %average over channels within the cluster
            hold on;
            plot(GA_young_mmf.time,mean(GA_young_mmf.avg(sort(unique(x)),:)),'m','LineWidth',5);
            ylim([-1e-15 1.5e-15]);
            patch([min(y)/1000 min(y)/1000 max(y)/1000 max(y)/1000],[min(ylim) max(ylim) max(ylim) min(ylim)],'k','FaceAlpha',0.1) %shade between time limits of cluster
            %xlim([-0.1 0.5]) %zoom in
            title(sprintf('Cluster Time:  %.3fs to %.3fs\n(p = %.3f)', ...
                time_for_topo(1),time_for_topo(end), stat.negclusters(t).prob));
            xlabel('Time (sec)');
            ylabel('Amplitude (Tesla/cm^{2})')
            legend('Older MMF','Younger MMF', 'Location','northwest')
            set(gca,'fontsize', 40);
            set(gcf,'position',[10,10,1400,1000])

            % Save as png
            if strcmp(save_to_file,'yes')
                disp('Saving figure to .png file');
                print(sprintf([output_path 'MMF_young-vs-old_neg_cluster_%d'], t),'-dpng');
            else
                disp('Not saving figure to file');
            end
                
            % Topoplot
            cfg                  = [];
            cfg.interpolation    = 'v4';
            cfg.marker           = 'off';
            cfg.highlight        = 'on';
            cfg.highlightchannel = stat.label(highlight_chan);
            cfg.highlightsymbol  = '.';
            cfg.highlightsize   = 20;
            cfg.xlim            = [time_for_topo(1) time_for_topo(end)];
            %   cfg.zlim            = 'maxabs';
            %cfg.zlim            = [-5 5];
            cfg.comment         = 'no';
            cfg.fontsize        = 6;
            cfg.layout          = lay;
            cfg.parameter       = 'stat';
            figure;colorbar
            ft_topoplotER(cfg,stat); hold on;
            ft_hastoolbox('brewermap', 1);
            colormap(flipud(brewermap(64,'RdBu'))) % change the colormap

            % Give the title
            title(sprintf('Cluster Time:  %.3fs to %.3fs\n(p = %.3f)', ...
                time_for_topo(1),time_for_topo(end), stat.negclusters(t).prob));

            set(gca,'fontsize', 20);

            % Save as png
            if strcmp(save_to_file,'yes')
                disp('Saving figure to .png file');
                print(sprintf([output_path 'MMF_young-vs-old_topoplot_neg_cluster_%d'], t), '-dpng','-r200');
            else
                disp('Not saving figure to file');
            end

        end
    end
end

% see https://au.mathworks.com/help/matlab/ref/plot.html#btzitot-LineSpec
% for plotting aesthetics

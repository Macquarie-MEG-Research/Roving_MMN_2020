%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ME125 (Roving MMF) SENSOR-LEVEL ANALYSIS
% NON-PARAMETRIC CLUSTER-BASED RANDOM PERMUTATION TESTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% 1. Add Fieldtrip and MQ_MEG_Scripts to your MATLAB path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

addpath /Users/42450500/Documents/MATLAB/fieldTrip/fieldtrip-20200409/; % change path if necessary
ft_defaults
disp('Adding Fieldtrip to your MATLAB path');

addpath /Users/42450500/OneDrive - Macquarie University/phd/other/projects/ME175/analysis/fieldTrip_scripts;
addpath(genpath('/Users/42450500/Documents/MATLAB/fieldTrip/MQ_MEG_Scripts-master'))
disp('Adding MQ_MEG_Scripts + to your MATLAB path');

global ft_default
ft_default.spmversion = 'spm12'; % Force SPM12, SPM8 doesn't go well with mac + 2017b
ft_defaults % This loads the rest of the defaults
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Set up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('lay.mat');
load ('neighbours_125.mat')

% PLEASE SPECIFY where to read MEG data from:
data_path = '..\\..\\ME125_roving_Phase1_data_37kids\\';

% PLEASE SPECIFY where to output results:
output_path = ['D:\\Judy\\RA_2020\\ARC_Roving_MMN\\Phase1_Results_young-vs-old\\']; % full path required on Windows, due to back-slash issues


% PLEASE SPECIFY the groups of participants (for comparison):
group.older   = {'2913' '2787' '2697' '2702' '2786' '2716' '2698' '2712' '2872' '2703' '2888' '2811' '2696' '2713' '2904' '2854' '2699' '2858'}; % 18 kids, >=5yo
group.younger = {'2724' '2642' '2866' '2785' '2793' '2738' '2766' '2687' '2629' '2897' '2683' '2695' '2739' '2810' '2632' '2667' '2875' '2912' '2681'}; % 19 kids, <5yo

group_list = {'younger', 'older'};


% Perform baseline correction before computing ERF?
DO_BASELINE = true;


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
        cd([folders(idx(sub)).name '\\ReTHM\\']);

        if ~DO_BASELINE % if no baseline correction needed, just load the saved ERFs (ie. ave)
            load('deviant_ave.mat')
            load('predeviant_ave.mat')
            
        else % otherwise, load the epoched data, perform baseline correction then compute ERF
            load('deviant.mat')
            load('predeviant.mat')

            % do baseline correction on individual epochs
            cfg = [];
            cfg.baseline = [-0.1 0];
            deviant = ft_timelockbaseline(cfg, deviant);

            cfg = [];
            cfg.baseline = [-0.1 0];
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

    %print(['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/png/', group_list{i},'_Child_MMN_GFP_planar'],'-dpng');
    print([output_path group_list{i} '_Child_GFP_planar'],'-dpng');


    % note to self: This looked wrong at first glance (stand_ave.grad.type?)
    % but it's actually fine. PSF: This is not a problem as the grad structure is static. I.e. the geometry of the sensors is not changing between conditions.


    for k = 1:length(avg_deviant_planar_all)
        avg_standard_planar_all{1,k}.grad.type = deviant_ave.grad.type;
        avg_deviant_planar_all{1,k}.grad.type = deviant_ave.grad.type;
    end

    %save (['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/variables/','avg_deviant_planar_all_',group_list{i}],'avg_deviant_planar_all')
    %save (['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/variables/','avg_standard_planar_all_',group_list{i}],'avg_standard_planar_all')
    %save (['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/variables/','avg_mmf_planar_all_',group_list{i}],'avg_mmf_planar_all')
    save ([output_path,'avg_deviant_planar_all_',group_list{i}],'avg_deviant_planar_all')
    save ([output_path,'avg_standard_planar_all_',group_list{i}],'avg_standard_planar_all')
    save ([output_path,'avg_mmf_planar_all_',group_list{i}],'avg_mmf_planar_all')


    %% Settings for NPCBRPT
    clear stat_wholeEpoch
    %     clear (['stat_wholeEpoch_',group_list{i}])
    cfg = [];
    %cfg.grad = deviant_ave.grad;
    cfg.channel     = 'all';
    cfg.neighbours  = neighbours; % defined as above
    cfg.latency     = [0 0.5];  % timewindow for the stats. Epoched earlier, but doesn't make sense to do the stats before 0, since an effect before the stimulus (i.e., before 0) would be meaningless/not interpretable.
    cfg.avgovertime = 'no'; %
    cfg.parameter   = 'avg';
    cfg.method      = 'montecarlo';
    cfg.statistic   = 'ft_statfun_depsamplesT';
    cfg.alpha       = 0.05; % threshold for the significant clusters when doing the Mote Carlo p-value comparison step
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

%TODO% 
% * Some errors to fix
% * plot average of significant sensors (currently just plotting t-values)


for i=1:length(group_list)

    % Load data generated from code above
    %addpath '/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/variables/';
    load ([output_path,'avg_standard_planar_all_',group_list{i}]);
    load ([output_path,'avg_deviant_planar_all_',group_list{i}]);
    load ([output_path,'stat_wholeEpoch_',group_list{i}]);


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

    alpha        = 0.05;
    x_lims       = [0 0.4];
    save_to_file = 'no';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ft_warning('This function only supports plotting POSITIVE clusters at present');

    % FIX HERE - only creates plots for positive clusters; but it does still compute\
    % variables for pos & neg clusters, just doesn't plot them all
    cd(orig);

    % Find the clusters under the specified alpha
    if isfield(stat,'posclusters')
        pos_cluster_pvals = [stat.posclusters(:).prob];


        % If this array is empty, return error
        if isempty(pos_cluster_pvals)
            error('NO POSITIVE CLUSTERS BELOW ALPHA LEVEL');
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

                pos_cluster_mask(t,:)=highlight_chan;

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
                    %print(sprintf(['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/png/',...
                    %    group_list{i},'_singleplot_pos_cluster_%d'],t),'-dpng','-r200');
                    print(sprintf([output_path, group_list{i}, '_singleplot_pos_cluster_%d'],t),'-dpng','-r200');

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
    
    %cd ('/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/subjects')
    cd(orig)
end


%% printing stats to the command window

% check the stats: finds all the clusters from initial phase, and checks sig. when compared to
% the null dist.
% run these lines separately

%load('/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/variables/stat_wholeEpoch_older.mat')
load([output_path 'stat_wholeEpoch_older.mat'])
old_stat=stat_wholeEpoch;
old_stat.posclusters.prob % where pos would = D > S
old_stat.negclusters.prob % where pos would = S > D

%load('/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/variables/stat_wholeEpoch_younger.mat')
load([output_path 'stat_wholeEpoch_younger.mat'])
young_stat=stat_wholeEpoch;
young_stat.posclusters.prob % where pos would = D > S. Needs to be < 0.05 to be sig.
young_stat.negclusters.prob % where pos would = S > D

% if error, probably = no cluster


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         GROUP COMPARISON - MMF in younger vs older kids
%                      (cluster-based t-test)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%cd '/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/variables/';
cd(output_path)

load('avg_mmf_planar_all_older.mat')
old_mmf = avg_mmf_planar_all;

load('avg_mmf_planar_all_younger.mat')
young_mmf = avg_mmf_planar_all;


cfg                  = [];
cfg.channel          = 'all';
cfg.neighbours       = neighbours; % defined as above
cfg.latency          = [0 0.5];
cfg.avgovertime      = 'no'; %
cfg.parameter        = 'avg';
cfg.statistic        = 'ft_statfun_indepsamplesT';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.05; % threshold for initial clustering before correction
cfg.method           = 'montecarlo';
cfg.alpha            = 0.05; % threshold for the significant clusters when doing the Mote Carlo p-value comparison step
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


stat_MMFbyGrp           = ft_timelockstatistics(cfg,old_mmf{:},young_mmf{:});

%save (['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/variables/','stat_MMFbyGrp'], 'stat_MMFbyGrp');
save ([output_path,'stat_MMFbyGrp'], 'stat_MMFbyGrp');


% command window will print # pos and neg clusters BEFORE the correction
% to find out how many remain after correction (run these lines
% individually):

stat_MMFbyGrp.posclusters.prob
stat_MMFbyGrp.negclusters.prob % error probably cuz no neg clusters


%% plot stats

stat=stat_MMFbyGrp;

% Find the clusters under the specified alpha
if isfield(stat,'posclusters')
    pos_cluster_pvals = [stat.posclusters(:).prob];

    % If this array is empty, return error
    if isempty(pos_cluster_pvals)
        error('NO POSITIVE CLUSTERS BELOW ALPHA LEVEL');
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

            pos_cluster_mask(t,:)=highlight_chan;

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



%% plot MMM old v young (with sig mask)

% grand mean/a for old MMF & young MMF
GA_old_mmf   = ft_timelockgrandaverage([],old_mmf{:});
GA_young_mmf = ft_timelockgrandaverage([],young_mmf{:});

% mask_param = stat.mask;

%% PLOT
% NB. Not plotting the GFP, only plotting the channels (averaged) where the
% effect was
%close all

% which cluster to plot
cluster_n = 1; % select neg cluster #1, which is marginally sig

[x,y] = find(ismember(stat.negclusterslabelmat, cluster_n)); %find members of neg cluster #1
figure;plot(GA_old_mmf.time,mean(GA_old_mmf.avg(sort(unique(x)),:)),'g','LineWidth',5) %average over channels within the cluster
%ylim([-1e-15 2e-15]);
hold on
plot(GA_young_mmf.time,mean(GA_young_mmf.avg(sort(unique(x)),:)),'m','LineWidth',5)
patch([min(y)/1000 min(y)/1000 max(y)/1000 max(y)/1000],[min(ylim) max(ylim) max(ylim) min(ylim)],'k','FaceAlpha',0.1) %shade between time limits of cluster
%xlim([-0.1 0.5]) %zoom in
title(sprintf('p = %.3f', stat.negclusters(cluster_n).prob));
%title(sprintf('Cluster Time:  %.3fs to %.3fs\nCluster Peak: %.3fs' ...
%    ,time_for_topo(1),...
%    time_for_topo(end),time_of_peak));
xlabel('Time (sec)');
ylabel('Amplitude (Tesla/cm^{2})')
legend('Older MMF','Younger MMF', 'Location','northwest')
set(gca,'fontsize', 40);
set(gcf,'position',[10,10,1400,1000])

%print(['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/png/','MMF_sigCluster_youngVold'],'-dpng');
print([output_path 'MMF_sigCluster_young-vs-old'],'-dpng');

% see https://au.mathworks.com/help/matlab/ref/plot.html#btzitot-LineSpec
% for plotting aesthetics
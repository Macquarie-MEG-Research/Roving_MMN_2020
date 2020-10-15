clear all;
close all;
path=pwd

subject_rethm = {'3372','3383',	'3391',	'3437',	'3494',	'3496',	'3531',	'3634',	'3652',	'3678',	'3679',	'3680'};% subjects with autism

subject_norethm = {}; % 0 subjects don't have ReTHM

subject = [subject_rethm,subject_norethm];% Totally 72 subjects

% exclude bad datasets with
%1.  no audio channel for trigger correction
Subj_triggerProblem ={};% 4 subjects '2713'	'2724'	'2810'

%2. too big head movement
    % 1) the co-registration is greater than 5 at three or more markers; 
    % or 2) the difference between pre and post marker positions are
    %       greater than 10 at three or more markers; 
    % or 3) ReTHM was not recorded properly;
Subj_bigHeadmovement = { '3678' '3372' '3496' '3679' }; %3 subjects 

%3. too many bad channels (cutoff is 20)
Subj_badChannels = {}; %2 subjects

Subj_badRethm = {}; %
Subj_PoorCoreg = {'3391'  '3652' '3680'};%

bad_subject = [Subj_triggerProblem,Subj_bigHeadmovement,Subj_badChannels,Subj_badRethm,Subj_PoorCoreg] % 47 subjects


[logic,ind]=ismember (bad_subject,subject);
act_subject= subject;

if logic==1
    act_subject(ind)=[];
end
num_subj = length(act_subject);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% planar gradients %%%%%%%%
disp('Loading Subject List')

avg_deviant_planar_all      = [];
avg_predeviant_planar_all   = [];

for sub = 1:length(act_subject)
    
    disp(['Processing Subject ' act_subject{sub}]);
    
    if ismember (act_subject{sub},subject_rethm)
    cd([path,'\', act_subject{sub}]); %cd([path, act_subject{i},'\ReTHM']); MEG
    elseif ismember (act_subject{sub},subject_norethm)
    cd([path,'\', act_subject{sub}]);
    end   
    load('deviant.mat');
    load('predeviant.mat');

    if sub == 1
        cfg                 = [];
        cfg.feedback        = 'yes';
        cfg.method          = 'triangulation';
        fff     = ft_prepare_neighbours(cfg, deviant);
    end

avg_deviant = ft_timelockanalysis([],deviant);
avg_predeviant = ft_timelockanalysis([],predeviant);



cfg = [];
cfg.neighbours      = fff;
cfg.planarmethod    = 'sincos';
avg_deviant_planar              = ft_megplanar(cfg, avg_deviant);
avg_predeviant_planar       = ft_megplanar(cfg, avg_predeviant);

cfg                         = [];
cfg.method                  = 'sum';
avg_deviant_planar          = ft_combineplanar(cfg,avg_deviant_planar);
avg_predeviant_planar       = ft_combineplanar(cfg,avg_predeviant_planar);

avg_deviant_planar_all{sub} = avg_deviant_planar;
avg_predeviant_planar_all{sub} = avg_predeviant_planar;


end

cfg = [];
cfg.channel   = 'all';
cfg.latency   = 'all';
cfg.parameter = 'avg';
GA_planar_preDev    = ft_timelockgrandaverage(cfg,avg_predeviant_planar_all{:});  
GA_planar_Dev       = ft_timelockgrandaverage(cfg,avg_deviant_planar_all{:});
% "{:}" means to use data from all elements of the variable


cfg = [];
cfg.method    = 'power'
cfg.channel   = 'all';
GA_planar_preDev_GFP = ft_globalmeanfield(cfg,GA_planar_preDev);
GA_planar_Dev_GFP = ft_globalmeanfield(cfg,GA_planar_Dev);

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_difference_planar_GFP = ft_math(cfg,GA_planar_Dev_GFP,GA_planar_preDev_GFP);
GA_difference_planar_GFP.avg = abs(GA_difference_planar_GFP.avg);


cfg                          = [];
cfg.xlim                     = [-0.1 0.4];
cfg.title                    = 'Global Field Power';
cfg.graphcolor               = 'brk';
cfg.linestyle       = '-k';
%cfg.graphcolor      = 'k';
cfg.linewidth       = 6;
%cfg.zlim            = [-4.9 4.9];
cfg.showlabels      = 'yes';
cfg.fontsize        = 6;

figure;
% ft_singleplotER(cfg,standard_1_GFP,standard_2_GFP,standard_3_GFP,standard_4_GFP,standard_5_GFP,deviant_GFP)
ft_singleplotER(cfg,GA_planar_preDev_GFP ,GA_planar_Dev_GFP); hold on; 
legend('Pre-deviant','Deviant')
title(['Global Field Power']); drawnow;
xlabel('Time (sec)');
ylabel('Mean of Global Field Power (planar gradients)');
set(gca,'fontsize', 24);
    
x0=10;
y0=10;
width=1400;
height=1000
set(gcf,'position',[x0,y0,width,height])

print('Child_MMN_GFP_planar','-dpng');

figure; plot(GA_difference_planar_GFP.time,abs(GA_difference_planar_GFP.avg))
[pks,locs,widths,proms]=findpeaks(abs(GA_difference_planar_GFP.avg),GA_difference_planar_GFP.time,'WidthReference','halfprom')
[max1,ind1] = max(pks);
pks(ind1)      = -Inf;
[max2, ind2] = max(pks);
pks(ind2)      = -Inf;
[max3,ind3] = max(pks);
pks(ind3)      = -Inf;
[max4, ind4] = max(pks);
pks(ind4)      = -Inf;
[max5, ind5] = max(pks);
pks(ind5)      = -Inf;


pk_M1=locs(ind1)
pk_M2 = locs(ind2)
%get the time window for each interested peak based on the full width at half predominance
tw_M1 = [locs(ind1)-widths(ind1)/2,locs(ind1)+widths(ind1)/2]
tw_M2 = [locs(ind2)-widths(ind2)/2,locs(ind2)+widths(ind2)/2]

load grad_trans_ChildMMN
load lay_ChildMMN

%define or load neighours here
cfg_neighb        = [];
cfg_neighb.channel  = 'all';%cfg.channel;%[17:40,59:100]
cfg_neighb.feedback = 'yes';
cfg_neighb.method = 'triangulation';     %triangulation  distance template
cfg.grad = grad_trans;
neighbours        = ft_prepare_neighbours(cfg_neighb, avg_deviant_planar_all{1,1});

% stats for M1
clear stat_planar_M1
cfg = [];
cfg.channel     = 'all';
cfg.neighbours  = neighbours; % defined as above
cfg.latency     = tw_M1;%[0.095 0.135] for M1;  [0.185 0.25] for MMN
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;

Nsub = num_subj;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat_planar_M1 = ft_timelockstatistics(cfg,avg_predeviant_planar_all{:},avg_deviant_planar_all{:});

save stat_planar_M1 stat_planar_M1;

% make a plot
cfg = [];
cfg.channel   = 'all';
cfg.latency   = tw_M1;
cfg.parameter = 'avg';
GA_planar_Std_M1       = ft_timelockgrandaverage(cfg,avg_predeviant_planar_all{:});  
GA_planar_Dev_M1       = ft_timelockgrandaverage(cfg,avg_deviant_planar_all{:});

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_planar_DevvsStd_M1 = ft_math(cfg,GA_planar_Dev_M1, GA_planar_Std_M1);

cfg = [];
%cfg.style     = 'blank';
cfg.layout    = lay;
cfg.highlight = 'on';
cfg.highlightchannel = find(stat_M1.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GA_planar_DevvsStd_M1)
title('M1: significant with cluster multiple comparison correction')

% stats for M2
clear stat_M2
cfg = [];
cfg.channel     = 'all';
cfg.neighbours  = neighbours; % defined as above
cfg.latency     = tw_M2;%[0.095 0.135] for M1;  [0.185 0.25] for MMN
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;

Nsub = num_subj;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat_planar_M2 = ft_timelockstatistics(cfg,avg_predeviant_planar_all{:},avg_deviant_planar_all{:});

save stat_planar_M2 stat_planar_M2;

% make a plot
cfg = [];
cfg.channel   = 'all';
cfg.latency   = tw_M2;
cfg.parameter = 'avg';
GA_planar_Std_M2       = ft_timelockgrandaverage(cfg,avg_predeviant_planar_all{:});  
GA_planar_Dev_M2       = ft_timelockgrandaverage(cfg,avg_deviant_planar_all{:});

cfg = [];
cfg.operation = 'subtract';
cfg.parameter = 'avg';
GA_planar_DevvsStd_M2 = ft_math(cfg,GA_planar_Dev_M2, GA_planar_Std_M2);

cfg = [];
%cfg.style     = 'blank';
cfg.layout    = lay;
cfg.highlight = 'on';
cfg.highlightchannel = find(stat_M2.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GA_planar_DevvsStd_M2)
title('M2: significant with cluster multiple comparison correction')

%define or load neighours here
cfg_neighb        = [];
cfg_neighb.channel  = 'all';%cfg.channel;%[17:40,59:100]
cfg_neighb.feedback = 'yes';
cfg_neighb.method = 'triangulation';     %triangulation  distance template
cfg.grad = grad_trans;
neighbours        = ft_prepare_neighbours(cfg_neighb, avg_deviant_planar_all{1,1});

% stats for wholeEpoch
clear stat_wholeEpoch
cfg = [];
cfg.channel     = 'all';
cfg.neighbours  = neighbours; % defined as above
cfg.latency     = [-0.1 0.4];
cfg.avgovertime = 'no';
cfg.parameter   = 'avg';
cfg.method      = 'montecarlo';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'cluster';
cfg.correcttail = 'prob';
cfg.numrandomization = 1000;
cfg.minnbchan        = 2; % minimal neighbouring channels

Nsub = num_subj;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number

stat_wholeEpoch = ft_timelockstatistics(cfg,avg_deviant_planar_all{:},avg_predeviant_planar_all{:});

save stat_wholeEpoch stat_wholeEpoch;


% make a plot
cfg = [];
cfg.highlightsymbolseries = ['*','*','.','.','.'];%
cfg.layout = lay;
cfg.contournum = 0;
cfg.markersymbol = '.';
cfg.alpha = 0.05;
cfg.parameter='stat';
%cfg.zlim = [-5 5];
ft_clusterplot(cfg,stat_wholeEpoch);



% Plot the results (planar gradients)
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path

% To plot a squence of topo plots equally spaced between 0.1 and 0.4 second
timestep                    = 0.05; % in 50 ms steps
sampling_rate               = deviant.fsample;
sample_count                = length(GA_planar_preDev.time);

j                           = [0:timestep:0.4];
% m                           = [1:timestep*sampling_rate:sample_count];  % temporal endpoints in MEEG samples
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path

        figure;
        for k = 1: length(j)-1
            cfg = [];
            cfg.comment         = 'no';
        %     cfg.marker          = 'off';
            cfg.layout       = lay;
        %     cfg.colorbar        = 'southoutside';
            cfg.style           = 'straight';

            subplot(2,round(length(j)-1)/2,k)
            cfg.xlim = [j(k) j(k+1)];
            cfg.zlim = 'maxabs';
            ft_topoplotER(cfg, GA_planar_Dev)
            colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
        %     colorbar
            title(['Time [', num2str(j(k)),' ', num2str(j(k+1)),']']); drawnow;
            hold on
        end
        h = suptitle (['Topographic plots of the deviant response for grand average data (planar gradients)']);
        set (h,'FontSize',12,'FontWeight','bold')
        set(gcf, 'Position', [200, 200, 1000, 1000])
        print(['Topo_deviantavg_planar_grand_average'],'-dpng');
        
        figure;
        for k = 1: length(j)-1
            cfg = [];
            cfg.comment         = 'no';
        %     cfg.marker          = 'off';
            cfg.layout       = lay;
        %     cfg.colorbar        = 'southoutside';
            cfg.style           = 'straight';

            subplot(2,round(length(j)-1)/2,k)
            cfg.xlim = [j(k) j(k+1)];
            cfg.zlim = 'maxabs';
            ft_topoplotER(cfg, GA_planar_preDev)
            colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
        %     colorbar
            title(['Time [', num2str(j(k)),' ', num2str(j(k+1)),']']); drawnow;
            hold on
        end
        h = suptitle (['Topographic plots of the standard response for grand average data (planar gradients)']);
        set (h,'FontSize',12,'FontWeight','bold')
        set(gcf, 'Position', [200, 200, 1000, 1000])
        print(['Topo_standardavg_planar_grand_average'],'-dpng');
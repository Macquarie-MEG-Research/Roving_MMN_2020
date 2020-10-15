path=pwd

%subject_phase1 = {'2629',	'2632',	'2681',	'2683',	'2687',	'2695', '2696',	'2697',	'2699',	'2702',	'2713',	'2716',	'2766',	'2785',	'2786',	'2793',	'2810',	'2872',	'2888',	'2913'	}; %

subject_rethm ={'3370',	'3374',	'3378',	'3380',	'3394',	'3402',	'3407',	'3409',	'3411',	'3412',	'3416',	'3418',	'3419',	'3421',	'3422',	'3423',	'3426',	'3427',	'3428',	'3429',	'3430',	'3433',	'3435',	'3438',	'3439',	'3440',	'3441',	'3443',	'3445',	'3448',	'3493',	'3505',	'3508',	'3573',	'3612',	'3664'}%36 subjects have ReTHM

subject_norethm = {}; % 0 subjects don't have ReTHM

subject = [subject_rethm,subject_norethm];% Totally 37 subjects

% exclude bad datasets with
%1.  no audio channel for trigger correction
Subj_triggerProblem ={};% 4 subjects '2713'	'2724'	'2810'

%2. too big head movement
    % 1) the co-registration is greater than 5 at three or more markers; 
    % or 2) the difference between pre and post marker positions are
    %       greater than 10 at three or more markers; 
    % or 3) ReTHM was not recorded properly;
Subj_bigHeadmovement = { '3419' '3428' '3378' '3409' '3411' '3412' '3418' '3427' '3430' '3493' '3573' '3664'}; %3 subjects

%3. too many bad channels (cutoff is 20)
Subj_badChannels = {'3407'}; %2 subjects

Subj_badRethm = {'3429' '3505'}; % no good data at Phase I 
Subj_PoorCoreg = {};


bad_subject = [Subj_triggerProblem,Subj_bigHeadmovement,Subj_badChannels,Subj_badRethm,Subj_PoorCoreg] % 47 subjects


[logic,ind]=ismember (bad_subject,subject);
act_subject= subject;
%act_subject = subject_phase1; 

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


pk_M1=locs(ind2)
pk_M2 = locs(ind1)
%get the time window for each interested peak based on the full width at half predominance
tw_M1 = [locs(ind2)-widths(ind2)/2,locs(ind2)+widths(ind2)/2]
tw_M2 = [locs(ind1)-widths(ind1)/2,locs(ind1)+widths(ind1)/2]

%below for 22 subjects
pk_M1=locs(ind2)
tw_M1 = [locs(ind3)-widths(ind5)/2,locs(ind5)+widths(ind5)/2]

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
cfg.highlightchannel = find(stat_planar_M1.mask);
cfg.comment   = 'no';
figure; ft_topoplotER(cfg, GA_planar_DevvsStd_M1)
title('M1: significant with cluster multiple comparison correction')

% stats for M2
clear stat_planar_M2
cfg = [];
cfg.channel     = 'all';
cfg.neighbours  = neighbours; % defined as above
cfg.latency     = tw_M2; %[0.2370    0.2570]%[0.095 0.135] for M1;  [0.185 0.25] for MMN  [0.2310    0.2610]
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
cfg.latency   = tw_M2; %[0.2370    0.2570]% 0.2360   0.2560
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
cfg.highlightchannel = find(stat_planar_M2.mask);
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
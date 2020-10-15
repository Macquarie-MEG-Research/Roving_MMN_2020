%% prepare data
close all

% if isempty(dir('*rethm.con'))
%     [file] =  dir('*denoise.con');
% else 
% [file] =  dir('*rethm.con');%('*rethm.con')  denoise.con;
% end

if isempty ( dir('*B1_denoise_rethm.con'))
    [file] =  dir('*B1_denoise.con'); % 
else
    [file] =  dir('*B1_denoise_rethm.con');
end

filename = file.name;
pathname = file.folder;

fname = filename;

% filename = file.name;
% pathname = file.folder;
% fname    = dir('*.hsp');
% fname    = fname.name(1:strfind(fname.name,'.')-1);

hdr                     = ft_read_header(filename,'dataformat','yokogawa_con'); %hdr to get Fs etc.

load ICA_rc.mat;     % load 'data' from Step 4  

%% interpolate bad channels

load grad_trans.mat
load lay

%define or load neighours here
cfg_neighb        = [];
cfg_neighb.channel  = 'all';%cfg.channel;%[17:40,59:100]
cfg_neighb.feedback = 'yes';
cfg_neighb.method = 'triangulation';     %triangulation  distance template
cfg.grad = grad_trans;
neighbours        = ft_prepare_neighbours(cfg_neighb, data);

badchannels            = setdiff(lay.label(1:125), data.label);

cfg = [];
cfg.layout              = lay;
cfg.method         = 'weighted'; %'weighted', 'average', 'spline', 'slap' or 'nan' (default = 'weighted')
cfg.missingchannel      = badchannels; %cell-array, see FT_CHANNELSELECTION for details
cfg.neighbours     = neighbours; %neighbourhood structure, see also FT_PREPARE_NEIGHBOURS
cfg.trials         = 'all';% or a selection given as a 1xN vector (default = 'all')
%cfg.lambda         = regularisation parameter (default = 1e-5, not for method 'distance')
%cfg.order          = order of the polynomial interpolation (default = 4, not for method 'distance')
cfg.grad          = grad_trans; %structure with gradiometer definition, see FT_DATATYPE_SENS
cfg.senstype            = 'meg';
data = ft_channelrepair(cfg, data)


%% Trigger-based trial selection
cfg                         = [];
cfg.dataset                 = filename;
cfg.path                    = pathname;
cfg.trialdef.prestim        = 0.1;         % pre-stimulus interval
cfg.trialdef.poststim       = 0.4;        % post-stimulus interval
cfg.trialdef.eventtype      = 'trigger';
cfg.trialdef.eventvalue     = [1 8];%Trigger numbers
cfg.Fs                  = hdr.Fs;
cfg.First_Channel       = 146; % First trigger channel  146
cfg.Last_Channel        = 158; % Last trigger channel   158
% cfg.Audio_Channel       = 135; % 0 for none  135
% cfg.First_Channel       = 194; % First trigger channel
% cfg.Last_Channel        = 200; % Last trigger channel
% cfg.Audio_Channel       = 167; % -1 for none
if ismember (filename(1:4),Subj_correctedTrig) 
    cfg.Audio_Channel       = -1; % 0 for none  135
    cfg.fixed_offset        = 42; %defualt = []; 42ms for those without audio channel. 
else   
    cfg.Audio_Channel       = 151; % 0 for none  135
    cfg.fixed_offset        = []; 
end
cfg.trialfun            = 'FindTriggers_AudioCh'; %
cfg                     = ft_definetrial(cfg);

alldata                     = ft_redefinetrial(cfg, data);

event                       = cfg.event;
trl                         = cfg.trl;

save event event
save trl trl

 
% Find trials with NaNs
reject_ind = [];
count = 1;

for i = 1: length(alldata.trial)
    aa = find(cell2mat(cellfun(@isnan, alldata.trial(i),'UniformOutput',false)));
    if ~isempty(aa)
        reject_ind(count,1) = i;
        count = count+1;
    end
end


% Get sequence of tones
event = alldata.cfg.event;
ggg = [];
for i = 1:length(event)
   ggg(i,1) = event(i).value;
end

% Find the deviant trials
deviant_trials = find(ggg == 1);

% Remove first trial
deviant_trials(1) = [];

% Now find deviant trials in which previous sequence length is greater than 2
fff = ggg(deviant_trials - 1);
deviant_trials2 = deviant_trials(find(fff > 2));

% Get trial before (predeviant)
predeviant_trials = deviant_trials2 - 1;

% Now remove trials with NaNs
deviant_trials2 = deviant_trials2(~ismember(deviant_trials2,reject_ind));
fprintf('Found %d deviant trials\n', length(deviant_trials2));

predeviant_trials = predeviant_trials(~ismember(predeviant_trials,reject_ind));
fprintf('Found %d predeviant trials\n', length(predeviant_trials));

% Select data

cfg = [];
cfg.trials = deviant_trials2;
deviant = ft_selectdata(cfg,alldata);

cfg.trials = predeviant_trials;
predeviant = ft_selectdata(cfg,alldata);

% Save 

disp('Saving data');
save deviant deviant
save predeviant predeviant

%{
reject_ind = [];
for i = 1: length(alldata.trial)
aa = find(cell2mat(cellfun(@isnan, alldata.trial(i),'UniformOutput',false)));
    if ~isempty(aa)
        reject_ind(i) = 0;
    elseif isempty(aa)
        reject_ind(i) = i;
    end
end

clean_trials = nonzeros (reject_ind)';
clean_event = event(clean_trials);
cfg = [];
cfg.trials = clean_trials;
clean_data = ft_redefinetrial(cfg,alldata);


%% Epoching 

        
        for j=1:length(clean_event)
            if clean_event(j).value == 1
                clean_event(j).type = 'deviant';
            else
                clean_event(j).type = 'standard';
            end
        end
        
        
        for j=1:length(clean_event)-1
            if strcmp(clean_event(j).type,'standard') & strcmp(clean_event(j+1).type,'deviant')
                clean_event(j).type   = 'predeviant';
                clean_event(j+1).type = clean_event(j+1).type;
            else
            end
        end

        trials_idx = [];
        for j=1:length(clean_event)
            if strcmp(clean_event(j).type,'deviant') | strcmp(clean_event(j).type,'predeviant')
            trials_idx(j)=j;
            else
                trials_idx(j)=0;
            end
        end
        
        
selected_trials = nonzeros (trials_idx)';
event_select = clean_event(selected_trials);

cfg = [];
cfg.trials = selected_trials;
selected_data = ft_redefinetrial(cfg,clean_data);

data_clean = selected_data;

%Visual artefact rejection is optional
cfg                         = [];
cfg.method                  = 'summary';
cfg.keepchannel             = 'yes';
cfg.keeptrial               = 'nan';
data_clean                  = ft_rejectvisual(cfg,selected_data);

save data_clean data_clean;
good_trials_idx             = find(~isnan(cell2mat(cellfun(@(isgood)isgood(1),data_clean.trial,'uni',0)))); %just need to evaluate the first element as all samples in bad trial are NaN

event_clean=event_select(good_trials_idx); %remove bad trials events
trl_clean=trl(good_trials_idx,:);

trigger_types  = {event_clean.type};
trigger_values = cell2mat({event_clean.value});


% % deviants
% deviant_trials               = find(ismember(trigger_types,'deviant')); %trigger number 1 is first deviant
% 
% cfg                          = [];
% cfg.trials                   = intersect(deviant_trials,good_trials_idx);
% deviant                      = ft_redefinetrial(cfg,data_clean);
% save deviant deviant
% 
% % pre-deviants/standards
% predeviants_trials               = find(ismember(trigger_types,'predeviant')); %trigger number 1 is first deviant
% 
% cfg                          = [];
% cfg.trials                   = intersect(predeviants_trials ,good_trials_idx);
% predeviant                  = ft_redefinetrial(cfg,data_clean);
% save predeviant predeviant
%}


%%Averaging
cfg                          = [];
% cfg.keeptrials               = 'yes';
deviant_ave                  = ft_timelockanalysis(cfg,deviant);
save deviant_ave deviant_ave

cfg                          = [];
% cfg.keeptrials               = 'yes';
predeviant_ave                  = ft_timelockanalysis(cfg,predeviant);
save predeviant_ave predeviant_ave


%% Some plotting and other steps - not requried for group analysis

% Calculate the difference 
cfg                          = [];
cfg.operation                = 'subtract';
cfg.parameter                = 'avg';
difference_ave               = ft_math(cfg, deviant_ave, predeviant_ave);

%%Planar gradient transfer 

% Calculate the planar gradient of the averaged data: 
cfg                         = [];
cfg.method                  = 'triangulation';
cfg.neighbours              = ft_prepare_neighbours(cfg, grad_trans); %To get neigbours of the pak for plotting

cfg.planarmethod            = 'sincos';
deviant_ave_planar          = ft_megplanar(cfg, deviant_ave);
predeviant_ave_planar       = ft_megplanar(cfg, predeviant_ave);

% Combine the horizontal and vertical components of the planar gradient 
cfg = [];
deviant_ave_planar_comb     = ft_combineplanar(cfg,deviant_ave_planar);
cfg = [];
predeviant_ave_planar_comb  = ft_combineplanar(cfg,predeviant_ave_planar);

% GFP for planar gradient
cfg                          = [];
cfg.method                   = 'power';
deviant_planar_GFP  = ft_globalmeanfield(cfg, deviant_ave_planar_comb);
predeviant_planar_GFP   = ft_globalmeanfield(cfg, predeviant_ave_planar_comb);

cfg                          = [];
cfg.operation                = 'subtract';
cfg.parameter                = 'avg';
difference_planar_GFP               = ft_math(cfg, deviant_planar_GFP , predeviant_planar_GFP );
difference_planar_GFP.avg = abs(difference_planar_GFP .avg);

cfg                          = [];
cfg.xlim                     = [-0.1 0.4];
cfg.title                    = 'Global Field Power';
cfg.graphcolor               = 'brk';
figure;
% ft_singleplotER(cfg,standard_1_GFP,standard_2_GFP,standard_3_GFP,standard_4_GFP,standard_5_GFP,deviant_GFP)
ft_singleplotER(cfg,predeviant_planar_GFP ,deviant_planar_GFP , difference_planar_GFP)
legend('pre-deviant','deviant', 'difference')
title(['GPF of the difference waveform in subject ', num2str(filename(1:4)),' (planar gradients)']); drawnow;
print('GFP_StdvsDev','-dpng');


% Plot the results (planar gradients)
ft_hastoolbox('brewermap', 1);         % ensure this toolbox is on the path

% To plot a squence of topo plots equally spaced between 0.1 and 0.4 second
timestep                    = 0.05; % in 50 ms steps
sampling_rate               = alldata.fsample;
sample_count                = length(deviant_ave_planar_comb.time);

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
            ft_topoplotER(cfg, deviant_ave_planar_comb)
            colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
        %     colorbar
            title(['Time [', num2str(j(k)),' ', num2str(j(k+1)),']']); drawnow;
            hold on
        end
        h = suptitle (['Topographic plots of the deviant response in subject ', num2str(fname(1:4)),' (planar gradients)']);
        set (h,'FontSize',12,'FontWeight','bold')
        set(gcf, 'Position', [200, 200, 1000, 1000])
        print(['Topo_deviantavg_planar'],'-dpng');
        
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
            ft_topoplotER(cfg, predeviant_ave_planar_comb)
            colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
        %     colorbar
            title(['Time [', num2str(j(k)),' ', num2str(j(k+1)),']']); drawnow;
            hold on
        end
        h = suptitle (['Topographic plots of the standard response in subject ', num2str(fname(1:4)),' (planar gradients)']);
        set (h,'FontSize',12,'FontWeight','bold')
        set(gcf, 'Position', [200, 200, 1000, 1000])
        print(['Topo_standardavg_planar'],'-dpng');

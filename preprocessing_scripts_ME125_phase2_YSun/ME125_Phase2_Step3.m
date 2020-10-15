%% Prepare Data for ICA

if isempty ( dir('*merged_B1_denoise_rethm.con'))
    [file] =  dir('*denoise.con'); % 
else
    [file] =  dir('*merged_B1_denoise_rethm.con');
end

filename = file.name;
pathname = file.folder;


hdr                     = ft_read_header(filename,'dataformat','yokogawa_con'); %hdr to get Fs etc.

% high-pass filter with 1 Hz
cfg                         = [];
cfg.trialfun                = 'ft_trialfun_general';  
cfg.channel                 = hdr.grad.label; 
cfg.continuous = 'yes';
cfg.hpfilter         = 'yes';
cfg.hpfilttype       = 'firws';
cfg.hpfreq     = 1;
cfg.hpfiltdf   = 1.5;
cfg.hpfiltwintype    = 'blackman';
cfg.hpfiltdir  = 'onepass-zerophase';
cfg.dftfreq                 = 50; % removal line noise
cfg.headerfile              = filename;
cfg.datafile                = filename;
data_ini                    = ft_preprocessing(cfg);

%lowpass filter
cfg                         = [];
cfg.lpfilter   = 'yes';
cfg.lpfilttype = 'firws';
cfg.lpfreq     = 40;
cfg.lpfiltdf   = 10;
cfg.lpfiltwintype    = 'blackman';
cfg.lpfiltdir  = 'onepass-zerophase';
data_ini                    = ft_preprocessing(cfg,data_ini);

%reject the artifacts and channels that have been marked in script of Step2 
%cd([pwd,'\TSPCA']);
load arft.mat 
load selChLabel.mat
load sat.mat

arft.artfctdef.reject          = 'nan';
data4ICA=ft_rejectartifact(arft, data_ini);

cfg = [];
cfg.channel = selChLabel;
data4ICA = ft_selectdata(cfg, data4ICA);

if size(sat)>0 
satLabel=sat.label;
satTime=sat.time;
        for i=1:length(satLabel)
        [lia,locb]  = ismember(satTime{i},data4ICA.time{1,1});
        data4ICA.trial{1,1}(:,locb)=NaN;
        end     
end

        
%% run ICA
cfg = [];
cfg.method  = 'runica';
cfg.channel = 'all'; 
datacomp = ft_componentanalysis(cfg, data4ICA);

if size(sat)==0
save(['ICA_noSat.mat'], 'datacomp');
elseif size(sat)>0
save(['ICA.mat'], 'datacomp');
end

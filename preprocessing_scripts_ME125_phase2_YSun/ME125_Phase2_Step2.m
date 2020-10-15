close all
%cd([pwd,'\TSPCA']);
load data_0.1hz.mat;
load sat; 
%load arft;
%load selChLabel;
%% reject artifacts and bad channels

% reject artifacts
% Mark the obvious artifacts manually in ft_databrowser.

if size(sat)>0 
satLabel=sat.label;
satTime=sat.time;
        for i=1:length(satLabel)
        [lia,locb]  = ismember(satTime{i},data.time{1,1});
        data.trial{1,1}(:,locb)=NaN;
        end 
end

display(['processing ', data.cfg.previous.headerfile(1:4)]);

%data = ft_selectdata(data, 'channel', selChLabel); % notice the function's input now with key-value pairs, rather than a cfg.


cfg=[];
cfg.viewmode='vertical';
arft=ft_databrowser(cfg,data);

save(['arft.mat'], 'arft');

arft.artfctdef.reject          = 'nan';
data=ft_rejectartifact(arft, data);

% % remove bad channels
cfg                         = [];
cfg.method                  = 'channel';
cfg.alim     = 1e-10;
cfg.keepchannel             = 'no';
cfg.keeptrial               = 'nan';
data = ft_rejectvisual (cfg,data);

selChLabel = data.label;
save ('selChLabel.mat','selChLabel');


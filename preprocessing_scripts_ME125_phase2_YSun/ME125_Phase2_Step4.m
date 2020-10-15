close all
clear data datacomp

%fname    = dir('*.hsp');
%fname    = fname.name(1:strfind(fname.name,'.')-1);


if isempty ( dir('*denoise_rethm.con'))
    [file] =  dir('*denoise.con'); % 
else
    [file] =  dir('*denoise_rethm.con');
end

filename = file.name;
pathname = file.folder;

fname = filename;


%cd([pwd,'\TSPCA']);
load data_0.1hz.mat
load selChLabel.mat
load arft.mat
load lay
load sat.mat

if size(sat)> 0
load ICA.mat
elseif size(sat)== 0
load ICA_noSat
end
%% %%% reject artifacts

cfg = [];
cfg.channel = selChLabel;
data = ft_selectdata(cfg, data);

arft.artfctdef.reject          = 'nan';
data=ft_rejectartifact(arft, data)

if size(sat)>0 
satLabel=sat.label;
satTime=sat.time;
        for i=1:length(satLabel)
            %if ismember(satLabel(i),data.label)
        [lia,locb]  = ismember(satTime{i},data.time{1,1});
        data.trial{1,1}(:,locb)=NaN;
            %end
        end
end


%%%%% reject IC of EOG 
%        plot the components to detect the artifacts
figure
        cfg=[];
        cfg.component = [1:49];       
        cfg.layout = ['lay.mat'];
        cfg.marker             = 'off';
        cfg.comment = 'no';
        ft_topoplotIC (cfg, datacomp);
        
        h = suplabel(['ICA of ', fname(1:4)],'t',[.1 .1 .84 .84]);
        set (h,'FontSize',16,'FontWeight','bold')
        set(gcf, 'Position', [200, 200, 1500, 800])
        screen2png('ICA')

        cfg          = [];
        cfg.channel  = [1:10]; % components to be plotted
        cfg.viewmode = 'component';
        cfg.layout   = ['lay.mat']; % specify the layout file that should be used for plotting
        ft_databrowser(cfg, datacomp);
        
        cfg          = [];
        cfg.channel  = 'all'; % components to be plotted
        cfg.viewmode = 'vertical';
        cfg.layout   = ['lay.mat']; % specify the layout file that should be used for plotting
        ft_databrowser(cfg, data);
        
        IC = input('input the components you want to remove in a []');
        %remove the bad components and backproject the data
        cfg = [];
        cfg.component = [IC];% to be removed component(s)
        cfg.demean = 'no';
        data_removeIC = ft_rejectcomponent(cfg, datacomp, data);
        
        cfg          = [];
        cfg.channel  = 'all'; % components to be plotted
        cfg.viewmode = 'vertical';
        cfg.layout   = ['lay.mat']; % specify the layout file that should be used for plotting
        ft_databrowser(cfg, data_removeIC);
        
        display('Compare the datasets with and without IC removment. If it is OK, press any key to continue.')
        pause
        
        
        data = data_removeIC;
        save ('ICA_rc.mat', 'data')
% cd('E:\YananSun\MMN\MMN_analysed_Wei_YSun\MMN_ICA_Stats\37_subj_planar_Final');
% Set the alpha
alpha = 0.05;

load("lay_ChildMMN");
 
% Find the clusters under this alpha
pos_cluster_pvals = [stat_wholeEpoch.posclusters(:).prob];
pos_signif_clust = find(pos_cluster_pvals < alpha);

for t = 1:length(pos_signif_clust)
    
    % Get the significant channels
    pos = ismember(stat_wholeEpoch.posclusterslabelmat, pos_signif_clust(t));
    highlight_chan = any(pos(:,:)');
    
    % Get the significant times
    index = (any(stat_wholeEpoch.posclusterslabelmat == pos_signif_clust(t)));
    time_for_topo = stat_wholeEpoch.time(index');
    
    pos = ismember(stat_wholeEpoch.posclusterslabelmat, pos_signif_clust(t));
    stat_wholeEpoch.index = pos;
    
    % Find the time of the peak
    cfg = [];
    cfg.latency = [time_for_topo(1) time_for_topo(end)];
    cfg.channel = stat_wholeEpoch.label(highlight_chan');
    data_for_peak = ft_selectdata(cfg,stat_wholeEpoch);
    
    avg_chan  = mean(data_for_peak.stat(:,:));
    time_of_peak = data_for_peak.time(find(max(avg_chan)==avg_chan));
    
    % Plot 
    cfg                 = [];
    cfg.channel         = stat_wholeEpoch.label(highlight_chan');
    cfg.maskparameter   = 'index';
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
    ft_singleplotER(cfg,stat_wholeEpoch); hold on; 
    scatter(time_of_peak,max(avg_chan),40,'filled','r');

    % Give the title
    title(sprintf('Cluster #%d\n Time = %.3fs to %.3fs\nPeak = %.3fs' ...
        ,t,time_for_topo(1),...
        time_for_topo(end),time_of_peak)); 
    
    xlabel('Time (sec)');
    ylabel('t-value');
    set(gca,'fontsize', 24);
    
    % Save as png
    %print(sprintf('Adult_MMN_cluster_erf_%d_omi',t),'-dpng','-r200');
    print(sprintf('Child_MMN_cluster_erf_%d_omi',t),'-dpng','-r200');
end

%% Plot topoplot

for t = 1:length(pos_signif_clust)
    
    % Get the significant channels
    pos = ismember(stat_wholeEpoch.posclusterslabelmat, pos_signif_clust(t));
    highlight_chan = any(pos(:,:)');
    stat_wholeEpoch.highlight_chan = highlight_chan;
    
    % Get the significant times
    index = (any(stat_wholeEpoch.posclusterslabelmat == pos_signif_clust(t)));
    time_for_topo = stat_wholeEpoch.time(index');
    
    pos = ismember(stat_wholeEpoch.posclusterslabelmat, pos_signif_clust(t));
    stat_wholeEpoch.index = pos;
    
    % Find the time of the peak
    cfg = [];
    cfg.latency = [time_for_topo(1) time_for_topo(end)];
    cfg.channel = stat_wholeEpoch.label(highlight_chan');
    data_for_peak = ft_selectdata(cfg,stat_wholeEpoch);
    
    avg_chan  = mean(data_for_peak.stat(:,:));
    time_of_peak = data_for_peak.time(find(max(avg_chan)==avg_chan));
    
    % Plot 
    cfg                  = [];
    cfg.interpolation    = 'v4';
    cfg.marker           = 'off';
    cfg.highlight        = 'yes';
    cfg.highlightchannel = stat_wholeEpoch.label(highlight_chan);
    cfg.highlightsymbol  = '.';
    cfg.highlightsize   = 20;
    cfg.xlim            = [time_for_topo(1) time_for_topo(end)];
    cfg.zlim            = 'maxabs';
    cfg.showlabels      = 'no';
    cfg.comment         = 'no';
    cfg.fontsize        = 6;
    cfg.layout          = lay;
    cfg.parameter       = 'stat';
    figure;
    ft_topoplotER(cfg,stat_wholeEpoch); hold on;
    ft_hastoolbox('brewermap', 1);
    colormap(flipud(brewermap(64,'RdBu'))) % change the colormap
    
    % Give the title
    title(sprintf('Cluster #%d\n Time = %.3fs to %.3fs\nPeak = %.3fs' ,t,time_for_topo(1),...
        time_for_topo(end),time_of_peak)); 
    
    set(gca,'fontsize', 24);
    
    % Save as png
    %print(sprintf('Adult_MMN_cluster_topo_%d_omi',t),'-dpng','-r200');
    print(sprintf('Child_MMN_cluster_topo_%d_omi',t),'-dpng','-r200');

end

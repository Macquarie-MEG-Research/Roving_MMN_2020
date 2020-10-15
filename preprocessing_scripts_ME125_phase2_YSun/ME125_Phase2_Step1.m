%% House-keeping
close all;
if isempty ( dir('*merged_B1_denoise_rethm.con'))
    [file] =  dir('*denoise.con'); % 
else
    [file] =  dir('*merged_B1_denoise_rethm.con');
end
    pathname = file.folder;
    filename = file.name;



hdr                     = ft_read_header(filename,'dataformat','yokogawa_con'); %hdr to get Fs etc.

%% PART 1: MEMES:  Co-registration between MEG, Head Surface, and Brain Sruface in the MRI-library
%         fname    = dir('*.hsp');
%         fname    = fname.name(1:strfind(fname.name,'.')-1);
%         pathname = [pwd,'/'];
%         elpfile  = [fname,'.elp'];
%         hspfile  = [fname,'.hsp'];
%         confile  = filename; % !!! if NOT ReTHM data change the con file name here !!!
%         mrkfile  = [fname,'_INI.mrk'];
%         bad_coil = {};% define bad_coil if applicable !!!
%         path_to_MRI_library         = 'E:/YananSun/MEG110/database_for_MEMES\'; % !!! mannually check !!!
%         
        %child_MEMES_old(pathname,elpfile,hspfile,confile,mrkfile,path_to_MRI_library,bad_coil,'yes');
        % leave MEMES to the source-level analysis.
%% filtering
%highpass filter
cfg                         = [];
cfg.trialfun                = 'ft_trialfun_general';  
cfg.channel                 = hdr.grad.label; 
cfg.continuous = 'yes';
cfg.hpfilter         = 'yes';
cfg.hpfilttype       = 'firws';
cfg.hpfreq     = 0.1; %FIXME should try something lower to appease reviewers 0.1
cfg.hpfiltdf   = 0.15;
cfg.hpfiltwintype    = 'blackman';
cfg.hpfiltdir  = 'onepass-zerophase';
cfg.dftfreq                 = 50; % removal line noise
cfg.headerfile              = filename;
cfg.datafile                = filename;
data                    = ft_preprocessing(cfg);

%lowpass filter
cfg                         = [];
cfg.lpfilter   = 'yes';
cfg.lpfilttype = 'firws'; %'but'  'firws'
cfg.lpfreq     = 40;
cfg.lpfiltdf   = 10;
cfg.lpfiltwintype    = 'blackman';
cfg.lpfiltdir  = 'onepass-zerophase';
data                    = ft_preprocessing(cfg,data);

%% Create layout file for later + save
cfg                         = [];
cfg.grad                    = data.grad; % struct containing gradiometer definition
lay                         = ft_prepare_layout(cfg, data); % creates a 2-D layout of the channel locations
ft_layoutplot(cfg);

%if ~exist('TSPCA', 'dir')
%       mkdir('TSPCA')
%end
    
savepath = [pathname];

save([savepath,'\lay.mat'], 'lay');
save([savepath,'\data_0.1hz.mat'], 'data');


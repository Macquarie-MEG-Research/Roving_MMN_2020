%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         ME125 (Roving MMF) SOURCE RECONSTRUCTION
%          FOR ALL PARCELLATIONS IN THE AAL ATLAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Set up
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% = PLEASE SPECIFY =

% (1) Add necessary paths
addpath(genpath([pwd '/../'])); % get access to all scripts in current repo
addpath(genpath('D:/Judy/GitHub/')); % path to MQ_MEG_Scripts & MEMES
addpath(genpath('C:/Users/43606024/Documents/MATLAB/fieldtrip-20190702/template/')); % path to FT templates

% (2) Run adult or child data?
thisrun = 'child';
% THIS SETTING IS IMPORTANT for a number of steps (MEMES etc), 
% not just for setting the appropriate paths below

% Which group(s) of participants to analyse? (see lists below)
%group_list = {'younger', 'older'}; % split kids into two groups by age
group_list = {'child'}; % all kids in one group
%group_list = {'adult'}; % all adults in one group
    
% Note: If you specify two groups here, e.g. {'younger', 'older'}, 
% then these two groups will also be compared with each other at the end.

% (3) Specify relevant paths below

% Where to read MEG data from:
data_path_child = 'D:/Judy/RA_2021/ARC_Roving_MMN/Paul-HDD/data/';
data_path_adult = 'D:/Judy/RA_2021/ARC_Roving_MMN/ME125_roving_adult_data/';

% Where to store results:
output_path_child = 'D:/Judy/RA_2021/ARC_Roving_MMN/SourceReconstruction_forML_child/';
output_path_adult = 'D:/Judy/RA_2021/ARC_Roving_MMN/SourceReconstruction_forML_adult/'; 

% Location of MRI database (needed for MEMES):
MRI_library_child = 'D:/Judy/MRI_databases/database_for_MEMES_child/';
MRI_library_adult = 'D:/Judy/MRI_databases/new_HCP_library_for_MEMES/';

% Where are the data located inside each subject folder?
MEG_folder_child = '/'; %'/ReTHM/';
MEG_folder_adult = '/';
    
% (4) The following lists are set up for ME125 (roving) Phase 1 -> change if analysing other studies

% Lists of child subjects
%{
group.older   = {'2913' '2787' '2697' '2702' '2786' '2716' '2698' '2712' '2872' '2703' '2888' '2811' '2696' '2713' '2904' '2854' '2699' '2858'}; % 18 kids, >=5yo
group.younger = {'2724' '2642' '2866' '2785' '2793' '2738' '2766' '2687' '2629' '2897' '2683' '2695' '2739' '2810' '2632' '2667' '2875' '2912' '2681'}; % 19 kids, <5yo
group.child = [group.older group.younger];
%}
folders_phase1 = dir([data_path_child '2*']);
folders_phase2 = dir([data_path_child '3*']);
group.child = vertcat({folders_phase1(:).name, folders_phase2(:).name});

% List of adult subjects
folders = dir([data_path_adult '2*']);
group.adult = vertcat({folders(:).name});
    

% = ADJUST THESE SETTINGS AS NECESSARY =

% Coreg settings
%{
if strcmp(thisrun, 'child')
    coreg_version = 'MEMES_5mm'; % select which version of MEMES results to use
else
    coreg_version = 'MEMES'; % select which version of MEMES results to use
end
%}
coreg_quality_check = false; % if 'true', will produce plots of headmodel/mesh/sensors/etc


% Use beamformer or mne for source analysis? Acceptable options: 'mne', 'lcmv'
source_method = 'lcmv';
fixedori = 'no'; % 'no' == free orientation (not supported in the mne method, will use default 'yes')
if strcmp(source_method, 'mne')
    fixedori = 'yes';
end

% Stats settings
alpha_thresh = 0.05;  % threshold for stats
x_lims       = [0 0.4];
%save_to_file = 'yes'; % save stats figures to file?

% Location of template brain & MRI library
mri = 'single_subj_T1.nii'; % standard brain from the MNI database
aal_atlas = 'ROI_MNI_V4.nii';
%aal_atlas = 'C:/Users/43606024/Documents/MATLAB/fieldtrip-20190702/template/atlas/aal/ROI_MNI_V4.nii';
    
% Create subfolder in the output location & set appropriate suffix for VE filename (used below)
ori_suffix = ''; % no suffix if using fixed orientation
if strcmp(fixedori, 'no')
    ori_suffix = '_freeori';
end
if strcmp(source_method, 'lcmv')
    suffix = ori_suffix;
else
    suffix = ['_' source_method ori_suffix];
end


% = Getting ready to start =
if strcmp(thisrun, 'child')
    data_path = data_path_child;
    output_path = output_path_child;
    path_to_MRI_library = MRI_library_child;
    MEG_folder = MEG_folder_child;  
elseif strcmp(thisrun, 'adult')
    data_path = data_path_adult;
    output_path = output_path_adult;
    path_to_MRI_library = MRI_library_adult;
    MEG_folder = MEG_folder_adult;
else
    error(sprintf('Please specify a valid run option: adult, child.\nScript terminated.\n'));
end

output_path = [output_path source_method ori_suffix '/'];

cd(data_path)
orig = cd;
%folders = dir('2*');

folders = group.child; % SPECIFY HERE


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Start the subject loop - run MEMES & create VE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Run MEMES for each subject
%{
for j=1:length(folders)
    SubjectID = folders(j).name;
    cd([orig, '/', SubjectID, MEG_folder])
    
    coreg_output = [pwd '/' coreg_version '/']; % where to store the output from MEMES

    % if headmodel etc haven't been generated, do this now
    if ~exist([coreg_output 'headmodel.mat'], 'file')
        
        % find digitisation files
        %{
        elp_file  = dir('*.elp'); % find the .elp file
        filename_base = elp_file.name(1:strfind(elp_file.name,'.')-1); % get the base filename (ie. remove suffix)
        elpfile = [filename_base, '.elp'];
        hspfile = [filename_base, '.hsp'];
        mrkfile = [filename_base, '_ini.mrk']; % choose which marker file to use
        confile = [filename_base, '_B1_denoise_rethm.con'];        
        %}
        elp_file  = dir('*.elp');
        elpfile = elp_file.name;
        hsp_file  = dir('*.hsp');
        hspfile = hsp_file.name;
        mrk_file  = dir('*INI.mrk'); % choose which marker file to use
        mrkfile = mrk_file.name;
        con_file  = dir('*.con'); 
        confile = con_file.name;

        % MANUALLY SPECIFY the bad marker coils (max 2) for each subject
        % Enter as: {'LPAred','RPAyel','PFblue','LPFwh','RPFblack'}
        bad_coil = ''; 
        if strcmp(SubjectID, '2818') || strcmp(SubjectID, '2853')
            bad_coil = {'LPAred', 'RPFblack'};
        elseif strcmp(SubjectID, '2673') || strcmp(SubjectID, '2766') || strcmp(SubjectID, '2786')
            bad_coil = {'LPAred'};
        elseif strcmp(SubjectID, '2874')
            bad_coil = {'LPFwh', 'RPFblack'};
        end

        % call MEMES3 (new version Oct 2019)
        if length(bad_coil) >= 2     % we use 'rot3dfit' by default (faster), however this sometimes has issues (creates upside-down coreg)
            realign_method = 'icp';  % when there are 2 or more bad coils, so in that case we use 'icp' instead
        else
            realign_method = 'rot3dfit';
        end
        
        if strcmp(thisrun, 'child')
            [headshape_downsampled] = downsample_headshape_child(hspfile); 
            [grad_trans] = mq_realign_sens(pwd,elpfile,hspfile,confile,mrkfile,bad_coil,realign_method);
            child_MEMES(pwd, grad_trans, headshape_downsampled, path_to_MRI_library, 3);
        elseif strcmp(thisrun, 'adult')
            cfg = [];
            %cfg.facial_info = 'no'; % did not collect facial points during digitisation
            [headshape_downsampled] = downsample_headshape_new(cfg,hspfile); 
            [grad_trans] = mq_realign_sens(pwd,elpfile,hspfile,confile,mrkfile,bad_coil,realign_method);
            MEMES3(pwd, grad_trans, headshape_downsampled, path_to_MRI_library, 'best', [0.95:0.01:1.05], 5, 3);
        end        
        % call MEMES3 (old version 2018)
        %MEMES3_old_2018(pwd, elpfile, hspfile, confile, mrkfile, MRI_folder, bad_coil, 'best', [0.99:0.01:1.01], 5, 'yes')

        % close the figures MEMES created (each subject creates 5
        % figures - becomes too many when running in batch)
        close all;

        % move the MEMES output into the coreg_output folder.
        if ~exist(coreg_output, 'dir')
            mkdir(coreg_output);
        end                
        movefile('*trans*', coreg_output);
        movefile('*shape*', coreg_output);
        movefile('*model*', coreg_output);
        movefile('*quality*', coreg_output);
        if strcmp(thisrun, 'child')
            movefile('*error_age*', coreg_output);
            movefile('MEMES_output.mat', coreg_output);
        elseif strcmp(thisrun, 'adult')
            movefile('*example*', coreg_output);
            movefile('*scaling*', coreg_output);
            movefile('*realigned*', coreg_output);
        end
    end
end
%}

%% Create VEs for each subject
for j=1:length(folders)
    %SubjectID = folders(j).name;
    SubjectID = folders{j};

    disp('Loading relevant data');
    
    cd([orig,'/',SubjectID,MEG_folder])
    load('deviant.mat');
    load('predeviant.mat');
    
    %cd (coreg_version);
    load('headmodel.mat');
    load('sourcemodel3d.mat');
    load('grad_trans.mat');
    
    if strcmp(thisrun, 'child')
        load('MEMES_output.mat');
    
        % Here we are loading the MRI chosen during MEMES coreg (ages 2-5 to 7-5)
        load([path_to_MRI_library '/' MEMES_output.MRI_winner...
            '/mri_realigned.mat']);

        % Transform this MRI based on the two matrices computed during MEMES coreg
        mri_realigned = ft_transform_geometry(MEMES_output.fid_matrix,...
            mri_realigned);
        mri_realigned = ft_transform_geometry(MEMES_output.trans_matrix,...
            mri_realigned);       
    elseif strcmp(thisrun, 'adult')
        load('mri_realigned_MEMES.mat');
        mri_realigned = mri_realigned_MEMES; 
        % No need to perform transformation here, as MEMES3 has already done
        % that for you
    end
    
    % Make a figure to check you've marked LPA and RPA the right way round(!)
    if coreg_quality_check
        ft_determine_coordsys(mri_realigned, 'interactive', 'no');
        hold on; % add the subsequent objects to the figure
        drawnow; % workaround to prevent some MATLAB versions (2012b and 2014b) from crashing
        ft_plot_vol(headmodel);
        ft_plot_sens(grad_trans);
    end
    
    temp                 = load('standard_sourcemodel3d5mm.mat'); 
    template_sourcemodel = temp.sourcemodel;
    template_sourcemodel = ft_convert_units(template_sourcemodel, 'mm');
    
    cfg                = [];
    cfg.grid.warpmni   = 'yes';
    cfg.grid.template  = template_sourcemodel; % standard sourcemodel
    cfg.grid.nonlinear = 'yes';
    cfg.mri            = mri_realigned; % individual mri
    sourcemodel        = ft_prepare_sourcemodel(cfg); % creates individual sourcemodel
    % (the grid points map 1-to-1 onto the template grid points, with the .pos field
    % specifying the actual coordinates of these grid points in subject space)
    sourcemodel        = ft_convert_units(sourcemodel,'mm');
    
    if coreg_quality_check
        figure;
        ft_plot_sens(grad_trans, 'style', '*b'); % plot the MEG sensor locations
        ft_plot_vol(headmodel, 'edgecolor', 'cortex'); alpha 0.4; % plot the single shell (i.e. brain shape)
        ft_plot_mesh(sourcemodel.pos(sourcemodel.inside,:)); % plot all vertices (ie. grid points) that are inside the brain
    end
    
    % Prepare Leadfield
    cfg            = [];
    cfg.grad       = grad_trans;
    cfg.headmodel  = headmodel; % individual headmodel (from coreg)
    cfg.reducerank = 2; % Should check this is appropriate - also check the rank of the data as we project out mouth artifacts earlier
    cfg.channel    = deviant.label; % use the actual channels present in our data (i.e. ensure that rejected sensors are also removed here)
    cfg.grid       = sourcemodel; % individual sourcemodel (warped from template grid)
    grid           = ft_prepare_leadfield(cfg); % sourcemodel + leadfield
    %lf = grid;  % computes the forward model for many dipole locations on a regular sourcemodel and stores it for efficient inverse modelling
    
    % make a figure of the single subject{i} headmodel, and grid positions
    if coreg_quality_check
        figure; hold on;
        ft_plot_vol(headmodel,  'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
        ft_plot_mesh(grid.pos(grid.inside,:),'vertexsize',20);
        ft_plot_sens(grad_trans, 'style', 'r*','edgealpha',0.3); view([90,90]);
        %print('lf_headmodel_sens','-dpng','-r100');
    end
    
    %% Compute covariance matrix
    cfg                  = [];
    cfg.covariance       = 'yes';
    cfg.vartrllength     = 2;
    if strcmp(source_method, 'mne')
        cfg.covariancewindow = [-0.1 0]; % calculate covariance matrix
                                         % on the timepoints before zero 
    else
        cfg.covariancewindow = [0 0.5];
    end
    avg_deviant          = ft_timelockanalysis(cfg, deviant);
    avg_standard         = ft_timelockanalysis(cfg, predeviant);
    
    % Make a dummy variable with covariance matrices averaged
    avg_combined     = avg_deviant;
    avg_combined.cov = (avg_deviant.cov + avg_standard.cov) ./ 2;
    
    %% Source reconstruction
    % perform source reconstruction
    cfg                   = [];
    cfg.channel           = deviant.label;
    cfg.grad              = grad_trans;
    cfg.method            = source_method;
    cfg.grid              = grid;
    cfg.headmodel         = headmodel;
    if strcmp(source_method, 'mne')
        cfg.mne.keepfilter   = 'yes';
        cfg.mne.fixedori     = fixedori;
        %cfg.mne.fixedori     = 'no'; % this setting is not supported
        cfg.mne.prewhiten = 'yes';
        cfg.mne.lambda    = 3;
        cfg.mne.scalesourcecov = 'yes'; 
    else
        cfg.lcmv.keepfilter   = 'yes';
        cfg.lcmv.fixedori     = fixedori;
        cfg.lcmv.projectnoise = 'yes';
        %cfg.lcmv.weightnorm    = 'nai';
        cfg.lcmv.lambda       = '5%';
    end
    sourceall             = ft_sourceanalysis(cfg, avg_combined);

    % source localisation for deviant and standard trials using the common
    % spatial filter (not required for ROI analysis)
    % cfg.lcmv.filter        = sourceall.avg.filter;
    % source_deviant         = ft_sourceanalysis(cfg, avg_deviant);
    % source_standard      = ft_sourceanalysis(cfg, avg_standard);
    
    
    % Load Atlas (contains parcellation of brain into regions/tissues/parcels)
    atlas = ft_read_atlas(aal_atlas);
    atlas = ft_convert_units(atlas, 'mm');% ensure that atlas and template_sourcemodel are expressed in the same units
    
    % Interpolate the atlas onto template sourcemodel (10mm grid),
    % because the atlas may not be at the same resolution as your grid
    % (e.g. you created a grid with 6000 vertices, but atlas may only have 2000 vertices)
    cfg              = [];
    cfg.interpmethod = 'nearest';
    cfg.parameter    = 'tissue';
    atlas_interpo    = ft_sourceinterpolate(cfg, atlas, template_sourcemodel);
    
    
    % Define our ROIs (can combine multiple parcels together to form one ROI)
    %ROIs = {{'Frontal_Inf_Oper_L';'Frontal_Inf_Tri_L'},{'Frontal_Inf_Oper_R';'Frontal_Inf_Tri_R'},...
    %    {'Temporal_Sup_L'},{'Temporal_Sup_R'},{'Heschl_L'},{'Heschl_R'}};
    ROIs = atlas.tissuelabel; % use all parcellations in the atlas
    
    % Frontal merged for IFG
    % Temporal_Sup_L = STG
    % Heschl = A1
    
    %ROIs_label = {'LIFG','RIFG','LSTG','RSTG','LA1','RA1'}; %Labels for the groupings
    ROIs_label = ROIs;
    
    %% centroid
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%create_virtual_sensor_Centroid appropriate
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%labels - can do dynamically
    for k = 1:length(ROIs)
        ROI_name = ROIs_label{k};
        disp([num2str(k) ': ' ROI_name])
        
        % skip small parcels which don't have enough points for the centroid method
        if strcmp(ROI_name, 'Vermis_1_2')
            continue;
        end
        
        % for this ROI, find a list of vertices that belong to it, and
        % extract the spatial filter for each vertex in cue window & target window
        vertices_all = []; % will hold a single list of all vertices (from all parcels belonging to this ROI)
        %for i = 1:length(ROIs{k})
            %indx         = find(ismember(atlas_interpo.tissuelabel, ROIs{k}{i})); % find index of the required tissue label
            indx         = find(ismember(atlas_interpo.tissuelabel, ROIs{k})); % find index of the required tissue label
            vertices     = find(atlas_interpo.tissue == indx); % find vertices that belong to this tissue label
            % add vertices from the current parcel to the overall list
            vertices_all = [vertices_all; vertices];
        %end
        
        % get the spatial filter (i.e. set of weights) for each vertex
        if strcmp(fixedori, 'no') % for free orientation, there are 3 sets of weights (1 on each axis) for each vertex
            vertices_filters = cat(3, sourceall.avg.filter{vertices_all}); 
        else % for fixed orientation, there is 1 set of weights for each vertex
            vertices_filters = cat(1, sourceall.avg.filter{vertices_all});
        end

        % create virtual sensor for this ROI
        if strcmp(fixedori, 'no') % free dipole orientation
            VE_S = create_virtual_sensor_freeori(ROI_name, vertices_all, vertices_filters, avg_combined, avg_standard, 1, 'centroid', headmodel, sourcemodel); 
            VE_D = create_virtual_sensor_freeori(ROI_name, vertices_all, vertices_filters, avg_combined, avg_deviant, 1, 'centroid', headmodel, sourcemodel); 
        else % fixed dipole orientation
            VE_S = create_virtual_sensor_Centroid(ROI_name, vertices_all, vertices_filters, avg_combined, avg_standard, 1, headmodel, sourcemodel);
            VE_D = create_virtual_sensor_Centroid(ROI_name, vertices_all, vertices_filters, avg_combined, avg_deviant, 1, headmodel, sourcemodel);
        end
        
        if ~isempty(VE_S) % successful
            ROI_activity_standard.(ROI_name) = VE_S;
        else
            fprintf(['No solution for ', ROI_name, ' in cue window.']);
        end
        
        if ~isempty(VE_D) % successful
            ROI_activity_deviant.(ROI_name) = VE_D;
        else
            fprintf(['No solution for ', ROI_name, ' in cue window.']);
        end
        
    end
    
    % VE = [];
    % VE.label = labels;
    % try
    %     VE.sampleinfo = data_clean.sampleinfo;
    % catch
    %     disp('No sampleinfo field');
    % end
    % VE.time  = data_clean.time;
    %
    %
    % % For every VE...
    % for i = 1:length(labels)
    %     fprintf('ROI: %10s done\n',labels{i});
    %     % Create VE using the corresponding filter
    %     for trial=1:(length(data_clean.trial))
    %         % Multiply the filter with the data for each trial
    %         VE.trial{trial}(i,:) = sourceall.avg.filter{i,1}(:,:)...
    %             *data_clean.trial{trial}(:,:);
    %     end
    % end
    %save(ROI_output_file, 'ROI_activity');
    
    
    %%%%%%%%%%%%%%%%%%BROKEN ABOVE
    %%%%%%
    
    %     %% Now create VEs using the computed filters
    %     labels = {'L_A1','R_A1','L_STG','R_STG','L_IFG','R_IFG'};
    %
    %     [VE_deviant]    = mq_create_VE(deviants,source_deviant,labels);
    %     [VE_standard] = mq_create_VE(standards,source_standard,labels);
    %
    %     disp('Saving data');
    %     save VE_deviant VE_deviant
    %     save VE_standard VE_standard
    
    
    %% Now produce VE plot for 6 ROIs
    
    % % Timelock analysis
    % cfg             = [];
    % VE_deviant_ERF  = ft_timelockanalysis(cfg,VE_deviant);
    % VE_standard_ERF = ft_timelockanalysis(cfg,VE_standard);
    
    % Get max value for all ROIs for consistent plotting
    % ylimmm = max([max(max(abs(VE_deviant_ERF.avg))) ...
    %     max(max(abs(VE_standard_ERF.avg)))]).*1.1;
    
%    labels = fields(ROI_activity_deviant);
    labels = {'Frontal_Inf_Oper_L','Frontal_Inf_Oper_R','Temporal_Sup_L','Temporal_Sup_R','Heschl_L','Heschl_R'}'; % select 6 out of 116 parcels
    for i = 1:length(labels)
        
        minyd(i) = min(min(min(ROI_activity_deviant.(labels{i}).avg)));
        maxyd(i) = max(max(max(ROI_activity_deviant.(labels{i}).avg)));
        
        minys(i) = min(min(min(ROI_activity_standard.(labels{i}).avg)));
        maxys(i) = max(max(max(ROI_activity_standard.(labels{i}).avg)));
        
    end
    
    maxy = max([maxyd;maxys])+0.3*max([maxyd;maxys]);
    miny = min([minyd;minys])-abs(0.3*min([minyd;minys]));
    
    % Plot the Figure
    figure;
    set(gcf, 'Position',  [100, 100, 800, 1600])
    % For each ROI
    limit_idx = [sort(repmat([1:2:length(labels)]',2,1)) sort(repmat([2:2:length(labels)]',2,1))];
    
    
    for i = 1:length(labels)       
        cfg            = [];
        cfg.channel    = labels{i};
        cfg.linewidth  = 6;
        %     cfg.ylim       = [miny(limit_idx(i,1)) maxy(limit_idx(i,2))];
        cfg.xlim       = [-0.1 0.4];
        cfg.showlabels = 'yes';
        %cfg.graphcolor      = [0 0 0;1/255*[190,190,190]];
        cfg.fontsize   = 6;
        cfg.parameter  = 'avg';
        subplot(3,2,i);
        ft_singleplotER(cfg,ROI_activity_deviant.(labels{i}),ROI_activity_standard.(labels{i})); % blue = deviant
        
        xlabel('Time (sec)');
        set(gca,'fontsize', 14);
        legend('deviant','predeviant')
        
        % Give the title
        title(sprintf('%s',labels{i}),'Interpreter', 'none','FontSize',18);
    end
    
    % Save this png
    print([output_path 'individual_subjects_VE_plots/' SubjectID '_VE'],'-dpng','-r200');
    
    close all
    
    % save the ROI activities
    VE_standard = [];
    VE_deviant  = [];
    
    labels = fields(ROI_activity_deviant);
    
    VE_standard.time   = ROI_activity_standard.(labels{1}).time;
    VE_standard.dimord = ROI_activity_standard.(labels{1}).dimord;
    
    VE_deviant.time   = ROI_activity_deviant.(labels{1}).time;
    VE_deviant.dimord = ROI_activity_deviant.(labels{1}).dimord;
    
    for i=1:length(labels)        
        VE_standard.label{i} = labels{i};
        VE_standard.avg(i,:) = ROI_activity_standard.(labels{i}).avg;
        
        VE_deviant.label{i} = labels{i};
        VE_deviant.avg(i,:) = ROI_activity_deviant.(labels{i}).avg;
    end
    
    VE_predeviant = VE_standard; % rename
    save([output_path SubjectID '_VE.mat'], 'VE_predeviant', 'VE_deviant');
    %save(['VE_predeviant' suffix '.mat'], 'VE_predeviant');
    %save(['VE_deviant' suffix '.mat'], 'VE_deviant');
    
    cd(orig)
end

%% THE END


%{

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Grand average VE (separately for young & old)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for q=1:length(group_list)
    
    disp('Now performing group analysis by age group');
    
    idx   = find(ismember({folders.name},group.(group_list{q})));
    
    % remove subjects with bad coreg
    if strcmp(thisrun, 'adult')
        idx([3,4,5,6,8,10]) = [];
    end
    
    % Load the data for all subjects into two arrays
    % Put all of the VEs into a structure: Virtual electrode for deviant, and standard,
    % (avg across all subjects)
    
    VE_deviant_all  = [];
    VE_standard_all = [];
    
    for sub=1:length(idx)        
        SubjectID = folders(idx(sub)).name;
        disp(SubjectID);
        cd([orig,'/',SubjectID,MEG_folder,coreg_version])

        load(['VE_deviant' suffix '.mat']); % load each individual's VE D & S, and add to the '..._all' variable
        load(['VE_standard' suffix '.mat']);
        
        VE_deviant_all{sub}  = VE_deviant;
        VE_standard_all{sub} = VE_standard;
        
        clear VE_deviant VE_standard % not sure this is necessary, it will be overwritten by next load (?)
        cd(orig)        
    end
    
    % Grandaverage BUT keep the individuals (so we can plot 95% confidence
    % intervals
    cfg                  = [];
    cfg.parameter        = 'avg';
    cfg.keepindividual   = 'yes';
    VE_deviant_grandavg  = ft_timelockgrandaverage(cfg,VE_deviant_all{:});
    VE_standard_grandavg = ft_timelockgrandaverage(cfg,VE_standard_all{:});
    
        
    %% Plot S, D and difference waveform
    
    cfg           = [];
    cfg.operation = 'subtract';
    cfg.parameter = 'individual';
    MMF_all       = ft_math(cfg,VE_deviant_grandavg,VE_standard_grandavg);
    
    save ([output_path,'MMF_all_',group_list{q}],'MMF_all')
    
%{    
% % % %     save (['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/variables/','MMF_all_',group_list{q}],'MMF_all')
% % % %     
% % % %     mask_param_all = [];
% % % %     
% % % %     cmap_1 = [0.3294 0.6706 0.3569];
% % % %     
% % % %     figure;
% % % %     set(gcf, 'Position',  [100, 100, 800, 1600]);
% % % %     
% % % %     subplots = [1 2 5 6 9 10];
% % % %     
% % % %     for i = 1:length(labels)
% % % %         % Calculate mean and 95% confidence intervals for deviants
% % % %         [mean_MMF, CI_MMF] = mq_get_confidence_cousineau(squeeze(...
% % % %             MMF_all.individual(:,i,:)));
% % % %         cfg                  = [];
% % % %         cfg.channel          = VE_deviant_all{1,1}.label{i};
% % % %         cfg.latency          = [0 0.5];
% % % %         %cfg.dim             = VE_deviant_all{1,i}.dim;
% % % %         cfg.method           = 'montecarlo';
% % % %         cfg.statistic        = 'ft_statfun_depsamplesT';
% % % %         cfg.parameter        = 'avg';
% % % %         cfg.correctm         = 'cluster';
% % % %         cfg.computecritval   = 'yes';
% % % %         cfg.numrandomization = 2000;  % NB. Only did 1000 for the sensor level. More computationally difficult to do >1000 for sensor-lvl since more channels to process. Can do way more for source.
% % % %         %cfg.clusteralpha    = 0.001;
% % % %         cfg.tail             = 0;    % Two sided testing
% % % %         cfg.alpha            = 0.05;
% % % %         % Design Matrix
% % % %         nsubj           = length(VE_deviant_all);
% % % %         cfg.design(1,:) = [1:nsubj 1:nsubj];
% % % %         cfg.design(2,:) = [ones(1,nsubj) ones(1,nsubj)*2];
% % % %         cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
% % % %         cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)
% % % %         %
% % % %         stat            = ft_timelockstatistics(cfg,VE_deviant_all{:},...
% % % %             VE_standard_all{:});
% % % %         %     cfg = [];
% % % %         %     cfg.parameter = 'stat';
% % % %         %     cfg.maskparameter = 'mask';
% % % %         %     cfg.linewidth = 4;
% % % %         %     subplot(3,2,i);ft_singleplotER(cfg,stat)
% % % %         %     ylabel('T-value');
% % % %         %     xlabel('Time (s)');
% % % %         %     title(sprintf('%s',labels{i}),'Interpreter', 'none','FontSize',18);
% % % %         %     set(gca,'FontSize',14);
% % % %         % This creates a line with the significant times for plotting (plotted at
% % % %         % -1.9e-13)
% % % %         mask_param                = double(stat.mask);
% % % %         mask_param(mask_param==0) = NaN;
% % % %         %         mask_param(mask_param==1) = -1.4e-13;
% % % %         mask_param(mask_param==1) = -.9e-13;
% % % %         mask_param_all(i,:) = mask_param;
% % % %         % Plot using boundedline
% % % %         subplot(6,2,subplots(i)+2);boundedline(x,mean_MMF,CI_MMF(2,:),...
% % % %             'alpha','transparency',0.3,'cmap',cmap_1);
% % % %         % Label and adjust lims
% % % %         if i>4
% % % %             xlabel('Time (sec)');
% % % %         end
% % % %         
% % % %         if logical(mod(i,2))
% % % %             ylabel('Amplitude (Tesla/cm^{2})');
% % % %         end
% % % %         
% % % %         %     ylim([-1.5e-13 1.5e-13]);
% % % %         % %     ylim([-1e-13 1e-13]);
% % % %         xlim([-0.01 0.51]);
% % % %         %yticks([-1e-13 0 1e-13])
% % % %         % Adjust FontSize
% % % %         set(gca,'fontsize', 14);
% % % %         % Give the subplot a title (ROI)
% % % %         %title(sprintf('%s',labels{i}),'Interpreter', 'none','FontSize',18);
% % % %         % Plot a line indicating the significant times computed earlier
% % % %         hold on; drawnow;
% % % %         plot([0:.001:0.5],mask_param,'-k','LineWidth',3);
% % % %         ax = gca;
% % % %         ax.YRuler.TickLabelFormat = '%.0f';
% % % %         % end
% % % %         %
% % % %         % print('MMN_VE_all2','-dpng','-r300');
% % % %         %% Plot Deviant and Standard with output from Statistics
% % % %         % colors
% % % %         cmap_2 = [0.9020 0.0824 0.2706; 0.3098 0.1686 1.0000];
% % % %         labels = {'Left_IFG','Right_IFG','Left_STG','Right_STG','Left_A1','Right_A1'};
% % % %         % mask_param_all2 = mask_param_all.*(-2.4e-13 / -1.4000e-13);
% % % %         %mask_param_all2 = mask_param_all.*(-2.2e-13 / -2.2000e-13);
% % % %         %figure;
% % % %         set(gcf, 'Position',  [100, 100, 800, 1600]);
% % % %         % For every ROI
% % % %         %for i = 1:length(labels)
% % % %         % Calculate mean and 95% confidence intervals for deviants
% % % %         [mean_deviant, CI_deviant] = mq_get_confidence_cousineau(squeeze(...
% % % %             VE_deviant_grandavg.individual(:,i,:)));
% % % %         
% % % %         % Calculate mean and 95% confidence intervals for standards
% % % %         [mean_standard, CI_standard] = mq_get_confidence_cousineau(squeeze(...
% % % %             VE_standard_grandavg.individual(:,i,:)));
% % % %         
% % % %         %     %FIX ME HERE &&&&&&
% % % %         %     % Plot using boundedline
% % % %         %     subplot(6,2,subplots(i));boundedline(x,mean_deviant,CI_deviant(2,:),...
% % % %         %         x,mean_standard,...
% % % %         %         CI_standard(2,:),'alpha','transparency',0.3,'cmap',cmap_2);
% % % %         
% % % %         subplot(6,2,subplots(i));plot(x,mean_standard,x,mean_deviant);
% % % %         % Label and adjust lims
% % % %         %xlabel('Time (sec)');
% % % %         if logical(mod(i,2))
% % % %             ylabel('Amplitude (Tesla/cm^{2})');
% % % %         end
% % % %         
% % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %         %%%%%%%%%%%%%%% Y LIMITS FOR EACH ROI %%%%%%%%%%%%%%%%%%%%
% % % %         %         ylim([-10e-5 5e-5]);
% % % %         limits = {([-9e-6 8e-6]),([-9e-6 8e-6]),([-3e-5 2e-5]),([-3e-5 2e-5]),([-10e-5 6e-5]),([-10e-5 6e-5])};
% % % %         ylim(limits{i});
% % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % %         
% % % %         xlim([0 0.51]);
% % % %         % Adjust FontSize
% % % %         set(gca,'fontsize', 14);
% % % %         % Give the subplot a title (ROI)
% % % %         title(sprintf('%s',labels{i}),'Interpreter', 'none','FontSize',18);
% % % %         % Plot a line indicating the significant times computed earlier
% % % %         hold on; drawnow;
% % % %         plot(0:.001:0.5,mask_param_all(i,:),'-k','LineWidth',3);
% % % %         ax = gca;
% % % %         ax.YRuler.TickLabelFormat = '%.0f';
% % % %     end
% % % %     
% % % %     print(['/Users/42450500/OneDrive - Macquarie University/phd/data/MEG/ME175/inUse/group/png/', group_list{q},'_VE_S-D-MMF'],'-dpng','-r300');
%}    
    %%     Plot S v D

    mask_param_all = [];
    %cmap   = [0.9020 0.0824 0.2706; 0.3098 0.1686 1.0000];
    cmap   = ['r'; 'b'];

    x = VE_deviant_grandavg.time;
    
    % colors
    labels_full = {'Left IFG','Right IFG','Left STG','Right STG','Left A1','Right A1'};
    
    figure;
    set(gcf, 'Position',  [100, 100, 800, 1600]);
    
    % For every ROI
    for i = 1:length(labels_full)

        % Calculate mean and 95% confidence intervals for deviants
        [mean_deviant, CI_deviant] = mq_get_confidence_cousineau(squeeze(...
            VE_deviant_grandavg.individual(:,i,:)));
        
        % Calculate mean and 95% confidence intervals for standards
        [mean_standard, CI_standard] = mq_get_confidence_cousineau(squeeze(...
            VE_standard_grandavg.individual(:,i,:)));
        
        %%% ADDED top
        cfg                  = [];
        cfg.channel          = VE_deviant_all{1,1}.label{i};
        cfg.latency          = x_lims;
        %cfg.dim             = VE_deviant_all{1,i}.dim;
        cfg.method           = 'montecarlo';
        cfg.statistic        = 'ft_statfun_depsamplesT';
        cfg.parameter        = 'avg';
        cfg.correctm         = 'cluster';
        cfg.computecritval   = 'yes';
        cfg.numrandomization = 2000;  % NB. Only did 1000 for the sensor level. More computationally difficult to do >1000 for sensor-lvl since more channels to process. Can do way more for source.
        %cfg.clusteralpha    = 0.001;
        cfg.tail             = 0;    % Two sided testing
        cfg.alpha            = alpha_thresh;
        % Design Matrix
        nsubj           = length(VE_deviant_all);
        cfg.design(1,:) = [1:nsubj 1:nsubj];
        cfg.design(2,:) = [ones(1,nsubj) ones(1,nsubj)*2];
        cfg.uvar        = 1; % row of design matrix that contains unit variable (in this case: subjects)
        cfg.ivar        = 2; % row of design matrix that contains independent variable (the conditions)
        %
        stat            = ft_timelockstatistics(cfg,VE_deviant_all{:},VE_standard_all{:});
        mask_param                = double(stat.mask);
        mask_param(mask_param==0) = NaN;
        %         mask_param(mask_param==1) = -1.4e-13;
        mask_param(mask_param==1) = -.9e-13;
        mask_param_all(i,:) = mask_param;   
        %%% ADDED bottom

        
        % Plot using boundedline
        %         subplot(3,2,i);plot(x,mean_standard,x,mean_deviant);
        % OR
        %     % Plot using boundedline (CIS
        h = subplot(3,2,i);
        boundedline(x,mean_deviant,CI_deviant(2,:),...
            x,mean_standard,...
            CI_standard(2,:),'alpha','transparency',0.3,'cmap',cmap);
        
        
        % Label and adjust lims
        xlabel('Time (sec)');
        ylabel('Amplitude (Tesla/cm^{2})')
        %ylim([-4.2e-13 3.5e-13]);
        %ylim([-2.3e-13 2.3e-13]);
        
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%% Y LIMITS FOR EACH ROI %%%%%%%%%%%%%%%%%%%%
        %         ylim([-10e-5 5e-5]);
        %limits = {([-13e-6 12e-6]),([-13e-6 12e-6]),([-4e-5 3e-5]),([-4e-5 3e-5]),([-1.5e-4 1.5e-4]),([-1.5e-4 1.5e-4])};
        %ylim(limits{i});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        xlim([-0.1 0.41]);
        
        % Adjust FontSize
        set(gca,'fontsize', 14);
        xticks([0:0.1:0.4]);
        
        %%% ADDED below       
        hold on; drawnow;
        plot([0:.001:0.4], mask_param, '-k', 'LineWidth',3);
        %ax = gca;
        %ax.YRuler.TickLabelFormat = '%.0f';
        %%% ADDED above
    
        % Give the subplot a title (ROI)
        title(sprintf('%s',labels_full{i}),'Interpreter', 'none','FontSize',18);
        
        % Line width & legend
        lines = findall(h, 'Type','line');
        set(lines(2:3), 'Linewidth',2); % line thickness
        lgnd = legend(flip(lines(2:3)), 'Deviant','Predeviant', 'Location','southwest'); 
        set(lgnd,'color','none'); % make the legend background transparent
        legend('boxoff') % remove the box around legend
        
    end
    
    print([output_path, 'VE_SvsD_', group_list{q}],'-dpng','-r300');
    
    close all
    %clear mean_deviant mean_standard
    
    cd(orig);

end % end of Step 3


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. young vs old statistical comparison
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ONLY RUN THIS SECTION IF there are two groups specified at the top
if length(group_list) ~= 2
    warning('Only one group of participants were specified for analysis. Not performing any group comparison.'); % this will terminate the script
    return;
end

load([output_path 'MMF_all_older.mat'])
old_mmf = MMF_all;

load([output_path 'MMF_all_younger.mat'])
young_mmf = MMF_all;

for i = 1:length(ROIs_label)
    cfg                  = [];
    cfg.channel          = young_mmf.label{i};
    cfg.latency          = x_lims;
    %cfg.dim             = VE_deviant_all{1,i}.dim;
    cfg.method           = 'montecarlo';
    cfg.statistic        = 'ft_statfun_indepsamplesT';
    cfg.parameter        = 'individual';
    cfg.correctm         = 'cluster';
    cfg.computecritval   = 'yes';
    cfg.numrandomization = 2000;  % NB. Only did 1000 for the sensor level. More computationally difficult to do >1000 for sensor-lvl since more channels to process. Can do way more for source.
    %cfg.clusteralpha    = 0.001;
    cfg.tail             = 0;    % Two sided testing
    cfg.alpha            = alpha_thresh;
    
    % Design Matrix
    design = [ones(1,size(old_mmf.individual,1)) 2*ones(1,size(young_mmf.individual,1))];
    cfg.design      = design;
    cfg.ivar        = 1; % row of design matrix that contains independent variable (the conditions)
    %%% EXPLANATION: 'ones' is used to create a matrix. e.g.1., 'ones(4)' creates
    %%% a 4x4 matrix of 1s. e.g.2., 'ones(1,4)' creates one row of four 1s. For
    %%% b/w-subjects, comparing a row of group A/old (1s) with group B/young (2s).
    %%% To create the 2s = 2*ones = 2s. Call this 'deisgn' into the cfg settings.
    %%% the ivar is 'group' (with two levels, young vs old

    stat_MMFbyGrpROI.(ROIs_label{i}) = ft_timelockstatistics(cfg,old_mmf,young_mmf);
    
end

save ([output_path 'stat_MMFbyGroup_ROI'], 'stat_MMFbyGrpROI');

% check for sig effects
for i = 1:length(ROIs_label)
    length(find(stat_MMFbyGrpROI.(ROIs_label{i}).mask))
end


%% plot MMF: young vs old

% PLEASE SPECIFY: plot adult along with kids?
% (note - the stats were done on younger vs older kids, these were not statistically compared with the adult data)
PLOT_ADULT = false;

if PLOT_ADULT
    load([orig '/../Phase1_Source_Results_adult/' source_method suffix '/MMF_all_adult_goodcoreg19.mat']);
    adult_mmf = MMF_all;
end

figure; set(gcf, 'Position',  [100, 100, 800, 1600]);

cmap   = [1 0 1; 0 1 0; 0.9 0.4 0.15]; % magenta, green, orange

% For every ROI
for i = 1:length(labels_full)

    x = young_mmf.time;

    % Calculate mean and 95% confidence intervals for deviants
    [mean_young, CI_young] = mq_get_confidence_cousineau(squeeze(...
        young_mmf.individual(:,i,:)));
    
    [mean_old, CI_old] = mq_get_confidence_cousineau(squeeze(...
        old_mmf.individual(:,i,:)));

    if PLOT_ADULT
        [mean_adult, CI_adult] = mq_get_confidence_cousineau(squeeze(...
            adult_mmf.individual(:,i,:)));
    end
 
    % to create the bar indicating sig time intervals 
    % (the sig effects are only for the younger vs older kids comparison!!)
    mask_param                = double(stat_MMFbyGrpROI.(ROIs_label{i}).mask);
    mask_param(mask_param==0) = NaN;
    %         mask_param(mask_param==1) = -1.4e-13;
    mask_param(mask_param==1) = -.9e-13;
    mask_param_all(i,:) = mask_param;   


    % Plot using boundedline
    %         subplot(3,2,i);plot(x,mean_standard,x,mean_deviant);
    % OR
    %     % Plot using boundedline (CIS
    h = subplot(3,2,i);
    if PLOT_ADULT
        boundedline(x,mean_young,CI_young(2,:),...
                    x,mean_old,CI_old(2,:),...
                    x,mean_adult,CI_adult(2,:),...
                    'alpha','transparency',0.3,'cmap',cmap); 
    else
        boundedline(x,mean_young,CI_young(2,:),...
                    x,mean_old,CI_old(2,:),...
                    'alpha','transparency',0.3,'cmap',cmap);    
    end
    
    % Label and adjust lims
    xlabel('Time (sec)');
    ylabel('Amplitude (Tesla/cm^{2})');        

    xlim([-0.1 0.4]);
    xticks([0:0.1:0.4]);

    % Adjust FontSize
    set(gca,'fontsize', 14);
    
    hold on; drawnow;
    plot([0:.001:0.4], mask_param, '-k', 'LineWidth',3);

    % Give the subplot a title (ROI)
    title(sprintf('%s',labels_full{i}),'Interpreter', 'none','FontSize',18);
    
    % Line width & legend
    lines = findall(h, 'Type','line');
    set(lines(2:end), 'Linewidth',2); % line thickness
    if PLOT_ADULT
        lgnd = legend(flip(lines(2:end)), 'Younger MMF','Older MMF','Adult MMF', 'Location','northwest');    
    else
        lgnd = legend(flip(lines(2:end)), 'Younger MMF','Older MMF', 'Location','northwest');    
    end
    set(lgnd,'color','none'); % make the legend background transparent
    legend('boxoff') % remove the box around legend
end

if PLOT_ADULT
    print([output_path, 'MMF_VE_young-old-adult'],'-dpng','-r300');
else
    print([output_path, 'MMF_VE_young-vs-old'],'-dpng','-r300');
end

%}
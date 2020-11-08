% Load coordinates of centroids from file, then warp into subject space
% (in preparation for calling create_virtual_sensor_Centroid)
%
% @param mri_realigned: correctly realigned MRI file
% 
% @output pos_grid: centres of mass in individual space
%
function pos_grid = load_centroids(templates_dir, mri_realigned)
    % Load centre of mass information (in mm)
    % These coordinates were generated using:
    % https://github.com/mingruixia/BrainNet-Viewer/blob/master/BrainNet_GenCoord.m
    centre_of_mass = load([templates_dir 'Node_AAL116.txt']);

    % convert mri to mm for consistency
    mri_realigned = ft_convert_units(mri_realigned, 'mm');

    % Align centres of mass into subject space 
    % I.e. we are performing inverse non-linear warping from MNI-->individual
    %cfg = [];
    %cfg.template = fullfile(templates_dir, 'single_subj_T1.nii'); % template brain in MNI space (matches AAL atlas)
    %cfg.nonlinear   = 'yes';
    norm = ft_volumenormalise([], mri_realigned); % transformation matrix from MNI <--> individual
    posback = ft_warp_apply(norm.params, centre_of_mass, 'sn2individual');
    pos_grid = ft_warp_apply(pinv(norm.initial), posback); % xyz-coordinates of each parcel's centre of mass, in individual space

    % convert back to cm
    pos_grid = pos_grid./10; 
end
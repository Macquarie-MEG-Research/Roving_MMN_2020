% create a virtual sensor for the given ROI
% use this function if you chose cfg.fixedori='no' in ft_sourceanalysis
% https://mne.tools/stable/auto_tutorials/source-modeling/plot_dipole_orientations.html
%
% step (1): compute the reconstructed source activity for each vertex 
% along all 3 dipole directions (i.e. obtaining 3 timecourses for this vertex)
% step (2): do a vector combination at every time sample, in order to 
% combine into 1 timecourse (only keep the vector magnitude, discard the orientation)
% step (3): collapse all vertices into one VE by taking a plain average
% over all vertices (the magnitude of activity is positive at all vertices)
%
% @param vertices_filters: 3 x channel x vertex (i.e. for each vertex
%                          in this ROI, 3 sets of weights are provided)
% @param erf_combined:     average erf across all 4 conditions
% @param erf:              separate erf for each condition
%
function VE = create_virtual_sensor_freeori(ROI_name, vertices_filters, erf_combined, erf, conds)
    VE = [];

    % each cycle handles one cond
    for i = conds
        all_vertices_timecourses = [];

        % each cycle handles one vertex
        for v = 1:size(vertices_filters, 3)
            % grab the 3 sets of weights for this vertex
            F = vertices_filters(:,:,v);

            % compute source timecourse in all 3 dipole orientations
            timecourses_xyz = F * erf.avg(:,:); % estimated source activity = filter * erf (i.e. s = w * X)
            timecourse_combined = vecnorm(timecourses_xyz); % vector combination at every time sample, 
                                                         % to obtain the length of the vector (i.e. absolute magnitude 
                                                         % of brain activity, regardless of orientation)
            % add the timecourse for this vertex to the list
            all_vertices_timecourses = [all_vertices_timecourses; timecourse_combined];
        end

        % take a plain average over all vertices
        % can also try PCA
        VE_timecourse = mean(all_vertices_timecourses);

        % put it into a timelock structure for later calling ft_timelockstatistics (in stats_ROI.m)
        VE.time = erf.time;
        VE.avg = VE_timecourse;
        VE.label = {ROI_name};
        VE.dimord = 'chan_time';
    end        
end 
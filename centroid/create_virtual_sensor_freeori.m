% create a virtual sensor for the given ROI
% use this function if you chose cfg.fixedori='no' in ft_sourceanalysis
% https://mne.tools/stable/auto_tutorials/source-modeling/plot_dipole_orientations.html
%
% step (1): compute the reconstructed source activity for each vertex 
% along all 3 dipole directions (i.e. obtaining 3 timecourses for this vertex)
% step (2): do a vector combination at every time sample, in order to 
% combine into 1 timecourse (only keep the vector magnitude, discard the orientation)
% step (3): collapse all vertices into one VE by taking a plain average
% over all vertices (the magnitude of activity is positive at all
% vertices), or a weighted average (i.e. centroid method)
%
% @param vertices_filters: 3 x channel x vertex (i.e. for each vertex
%                          in this ROI, 3 sets of weights are provided)
% @param erf_combined:     average erf across all 4 conditions
% @param erf:              separate erf for each condition
% @param ROI_method:       'plainavg' or 'centroid'
% @param headmodel,sourcemodel: only required if using 'centroid' method
%                               (otherwise can just supply [])
%
function VE = create_virtual_sensor_freeori(ROI_name, vertices, vertices_filters, erf_combined, erf, conds, ROI_method, headmodel, sourcemodel)
    VE = [];
    
    if strcmp(ROI_method, 'centroid')
        % compute the coordinates of the centroid
        vertices_coords = get_coordinates_for_vertices(sourcemodel, vertices); % coordinates are in cm
        centroid = find_centroid(vertices_coords); % compute the centroid of this ROI 
        %plot_ROI_centre_of_mass(ROI_name, vertices_coords, centroid, headmodel); % for quality check
    end
    
    
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

        if strcmp(ROI_method, 'plainavg')
            % take a plain average over all vertices
            VE_timecourse = mean(all_vertices_timecourses);
        elseif strcmp(ROI_method, 'centroid')
            % weight the timecourse for each vertex by its distance to centre
            timecourses_weighted = [];
            for vertex = 1:size(vertices_coords,1) % each cycle processes one vertex
                distance = pdist([vertices_coords(vertex,:).*10 ; centroid.*10]); % calculate distance to centre in mm
                timecourses_weighted(vertex,:) = exp((-distance.^2)./400) .* all_vertices_timecourses(vertex,:); % weight
            end
            % take the mean of all the weighted timecourses (i.e. collapse into one timecourse)
            VE_timecourse = mean(timecourses_weighted, 1); 
        else
            error('Error: invalid ROI method. Supported options: plainavg, centroid.\n');
        end
        % can also try PCA method
        
        % put it into a timelock structure for later calling ft_timelockstatistics (in stats_ROI.m)
        VE.time = erf.time;
        VE.avg = VE_timecourse;
        VE.label = {ROI_name};
        VE.dimord = 'chan_time';
    end        
end 
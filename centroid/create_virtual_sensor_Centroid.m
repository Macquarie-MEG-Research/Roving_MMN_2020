% create a virtual sensor for the given ROI (based on the spatial filter
% for each vertex within the ROI), to estimate its activity
%
% uses Centroid method for collapsing all vertices within the ROI into a single VE
%
% @param vertices:      a list of vertices within this ROI (specified by their vertex index)
% @param vertices_filters: a list of vertices within this ROI, each entry contains the spatial filter (set of weights) for one vertex
% @param erf_combined:  average erf across all 4 conditions
% @param erf:           separate erf for each condition
%
function VE = create_virtual_sensor_Centroid(ROI_name, vertices, vertices_filters, erf_combined, erf, conds, headmodel, sourcemodel)

    % compute the coordinates of the centroid
    vertices_coords = get_coordinates_for_vertices(sourcemodel, vertices); % coordinates are in cm
    centroid = find_centroid(vertices_coords); % compute the centroid of this ROI 
    %plot_ROI_centre_of_mass(ROI_name, vertices_coords, centroid, headmodel); % for quality check


    % Create VE by collapsing the activities at all vertices within ROI into one timecourse,
    % weighting the activity at each vertex based on its distance from centre of mass
    %
    %"To generate a single regional timecourse, individual voxel signals
    % were weighted according to their distance from the centre of mass"
    % Brookes et al., (2016)

    VE = [];

    % if there are any vertices in this ROI, create virtual sensor to represent this ROI
    if ~isempty(vertices_filters)
        for i = conds
            % put it into a timelock structure for later calling ft_timelockstatistics (in stats_ROI.m)
            VE.time = erf.time;
            VE.label = {ROI_name};
            VE.dimord = 'chan_time';

            % a list of reconstructed source activities, one for each vertex
            timecourses = vertices_filters(:,:) * erf.avg(:,:); % estimated source activity = filter * erf (i.e. s = w * X) 

            % weight the timecourse for each vertex by its distance to centre
            timecourses_weighted = [];
            for vertex = 1:size(vertices_coords,1) % each cycle processes one vertex
                distance = pdist([vertices_coords(vertex,:).*10 ; centroid.*10]); % calculate distance to centre in mm
                timecourses_weighted(vertex,:) = exp((-distance.^2)./400) .* timecourses(vertex,:); % weight
            end

            % take the mean of all the weighted timecourses (i.e. collapse into one timecourse)
            VE.avg = mean(timecourses_weighted, 1);
        end

        % Preserve .sampleinfo field to avoid warnings later
        %VE.sampleinfo = data.sampleinfo;
    end    
end    
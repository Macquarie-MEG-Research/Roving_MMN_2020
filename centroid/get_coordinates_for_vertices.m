% retrieve the xyz-coordinates of the given vertices
function coordinates = get_coordinates_for_vertices(sourcemodel, vertices)
    % each cycle processes one vertex
    for i = 1:length(vertices)
        coordinates(i,:) = sourcemodel.pos(vertices(i), :); %take the warped coordinates
    end

    % alternative method:
    %coordinates = sourcemodel.pos(vertices,:);
end
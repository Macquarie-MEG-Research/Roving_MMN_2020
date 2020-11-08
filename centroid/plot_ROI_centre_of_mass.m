% plot vertices & centroid in the given ROI
function plot_ROI_centre_of_mass(ROI_name, vertices_coords, centroid, headmodel_singleshell)

    figure('Name',ROI_name, 'NumberTitle','off'); hold on;
    ft_plot_vol(headmodel_singleshell,  'facecolor', 'cortex', 'edgecolor', 'none');
    alpha 0.5; camlight; hold on;

    try           
        ft_plot_mesh(vertices_coords, 'vertexcolor','blue', 'vertexsize',10); hold on;
        ft_plot_mesh(centroid, 'vertexcolor','red', 'vertexsize',30); hold off;
    catch
        fprintf('In ROI "%s": Cannot plot vertices & centroid.\n', ROI_name);
    end
end
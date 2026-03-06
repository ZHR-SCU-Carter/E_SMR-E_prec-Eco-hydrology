function [SM_check, date_check, time_line, years, len] = extractAndVisualizeSMData(Site_Name, check_element, theta_crisis)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % - This script is designed for the visual inspection of soil moisture data across different sites in the "SiteAll" database. 
    % 
    % - It also provides a visualization of the relationship between soil moisture and the critical threshold θ_cri.
    % 
    % - The user is required to manually input three variables: Site_Name, check_element, and theta_crisis.
    %
    % - Site_Name refers to the name assigned to each site within the SiteAll database.
    %
    % - The SiteAll database is a master structure composed of multiple site-specific structures. 
    %   Each site structure includes geographic coordinates, yearly records, timestamps in datetime format, 
    %   daily averaged data, and corresponding element names.
    %
    % - The check_element corresponds to the element name in the SiteAll database, such as SW0.1 or SW0.2, 
    %   which represent soil moisture at depths of 10 cm and 20 cm, respectively.
    %
    % - The theta_crisis value is determined in advance using the script "1 - theta_cri_determination.m".
    % 
    % - The results returned by this function are subsequently used as input for "3 - analyzeESmrprecEffect.m" 
    %   to compute E_SMR and E_prec through the partitioning of spring soil moisture amplification.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %
    % Load data
    load SiteAll.mat

    % Prepare the full time axis
    date_check = (SiteAll.(Site_Name)(1).timestr(1):SiteAll.(Site_Name)(end).timestr(end))';
    time_line = date_check;
    years = [SiteAll.(Site_Name).year];
    len = numel(years);
    
    % Extract the element data
    SM_check = [];
    element_index = find(strcmp(SiteAll.(Site_Name)(1).elementName, check_element));
    for ii = 1:length(SiteAll.(Site_Name))
        SM_check = [SM_check; SiteAll.(Site_Name)(ii).data(:, element_index)];
    end

    % Visualization
    figure('OuterPosition', [42.6, 622.6, 1836.8, 508.8])
    plot(date_check, SM_check);
    title(strcat(Site_Name, '--', check_element));
    hold on
    yline(theta_crisis, '--');

    % Fill areas above threshold
    above_threshold = SM_check > theta_crisis;
    regions = bwconncomp(above_threshold);
    for k = 1:regions.NumObjects
        region_indices = regions.PixelIdxList{k};
        fill_dates = date_check(region_indices);
        fill_moisture = SM_check(region_indices);
        x_vertices = [fill_dates', flip(fill_dates')];
        y_vertices = [fill_moisture', theta_crisis * ones(size(fill_moisture'))];
        patch('XData', x_vertices, 'YData', y_vertices, 'FaceColor', 'g', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
end

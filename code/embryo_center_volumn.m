close all; clear; clc;

%%
counter = 1;

home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\InSitu_preImpEmbryo\';
cd(home_dir);
addpath(genpath('scripts/'))
species = 'mouse';
addpath(genpath('\\sodium\broad_thechenlab\anusinha\embryo_geometry'));
output_directory = sprintf('%s/embryo_geometry/',home_dir);
mkdir(output_directory);

plate = 203;
well = 1;
%embryos = [10 11 12 15:17 20 24 26 28 29 31 34];%plate 102 well 3
%embryos = [6 15 16 21];%plate 102 well 4
%embryos = [1 2 3 10 13];%plate 201 well 1
%embryos = [2 3 11];%plate 201 well 2
%embryos = [2 4 12];%plate 202 well 1


%blast
%embryos = [10];%plate 102 well 3
%embryos = [10 12 22];%plate 102 well 4
embryos = [2 12 13 14 15 16 8];%plate 203 well 1

expt = '20220726';%'20220501'; 

for embryo=embryos
    % Add a name for the first column
    p_w_e = sprintf('p_%d_w_%d_e_%d', plate, well, embryo);
    
    disp(sprintf('----- Processing embryo ... %s', p_w_e));

                
    embryo_dir = sprintf('%s/projects/%s/%s/transcriptome/plate_%d_well_%d/embryo_%d',home_dir,species,expt,plate,well,embryo);
    embryo_bounds = dlmread(sprintf('%s/processed/bounds.txt',embryo_dir));

    % Parameters to set
    ReReg_dir = sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/ReReg',home_dir,species,expt, plate, well, embryo);

        if exist(ReReg_dir, 'dir')
            time = datestr(datetime(now,'ConvertFrom','datenum'));
            disp(sprintf('%s ----- Loading  segmentation mask from ReReg dir ...', time));
            segmentation_mask_directory = ReReg_dir;
        else
            time = datestr(datetime(now,'ConvertFrom','datenum'));
            disp(sprintf('%s ----- Loading  segmentation mask  ...', time));
            segmentation_mask_directory = sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/',home_dir,species,expt, plate, well, embryo);
        end 

    VOXEL_SIZE_XY = 0.1625;
    VOXEL_SIZE_Z = 0.3986799;

    SURFACE_AREA_CONNECTIVITY = 6;

    % Polygonalization parameters
    EROSION_KERNEL_RADIUS = 0.5;
    EROSION_KERNEL_SIZE = 1 + 2*floor([EROSION_KERNEL_RADIUS/VOXEL_SIZE_XY, EROSION_KERNEL_RADIUS/VOXEL_SIZE_XY, EROSION_KERNEL_RADIUS/VOXEL_SIZE_Z]);
    EROSION_KERNEL = ones(EROSION_KERNEL_SIZE);

    POLYGONALIZATION_NUM_FACES = 100;


    % Functions to compute
    COMPUTE_DIMENSIONS = 1;
    COMPUTE_VOLUMES = 1;
    COMPUTE_CELL_CENTROIDS = 1;
    COMPUTE_SURFACE_AREA = 1;
    COMPUTE_CONTACT_AREA = 1;
    COMPUTE_REDUCED_REPRESENTATION = 1;
    COMPUTE_EMBRYO_CENTROIDS = 1;



    %%

        %disp(['Starting embryo ' num2str(currentMask) ' out of ' num2str(num_embryos) ' at ' char(datetime)]);    
        %disp(['Embryo name: ' segmentation_mask_filenames{currentMask}]);    
        disp(['-Loading data and extracting cell masks...']);

        % expected time: ~20 seconds
        tic

        % Load segmentation mask
        %segmentation_mask = load3DImage_uint8([segmentation_mask_directory segmentation_mask_filenames{currentMask}]);
        segmentation_mask = read_3d_tif(sprintf('%s/p%dw%de%d-membSeg-relabeled.tif',segmentation_mask_directory,plate,well,embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));

        % Extract metadata descriptor from filename
        metadata = p_w_e;

        % Get size of mask
        [xdim, ydim, zdim] = size(segmentation_mask);

        % Get values in label matrix (cell id numbers) as row vector
        cell_ids = unique(segmentation_mask(:));
        cell_ids = cell_ids(cell_ids ~= 0)';

        num_cells = size(cell_ids, 2);
        
    
        % Extract masks for individual cells
        individual_cell_masks = cell(num_cells,1);
        for currentCell = cell_ids
            disp(['--Processing mask/borders for cell ' num2str(currentCell) ' out of ' num2str(num_cells)]);
            individual_cell_masks{currentCell} = (segmentation_mask == currentCell);
        end     

        % CELL VOLUMES
        cell_volumes_voxels = zeros(num_cells, 1);
        cell_volumes_um3 = zeros(num_cells, 1);

        if COMPUTE_VOLUMES
            for currentCell = cell_ids
                disp(['--Computing volumes for cell ' num2str(currentCell) ' out of ' num2str(num_cells)]);
                [volume_voxels, volume_um3] = compute_volume(individual_cell_masks{currentCell}, VOXEL_SIZE_XY, VOXEL_SIZE_Z);
                cell_volumes_voxels(currentCell) = volume_voxels;
                cell_volumes_um3(currentCell) = volume_um3;
                
            end

            T_volumes_table = table(cell_ids', cell_volumes_voxels, cell_volumes_um3);
            T_volumes_table.Properties.VariableNames = {'cell_id', 'cell_volume_voxel', 'cell_volume_um3'};

            writetable(T_volumes_table, [output_directory metadata '_volume.csv'], 'WriteVariableNames', true);        
        end


        % CELL CENTROIDS
        cell_centroids_voxels = zeros(num_cells, 3);
        cell_centroids_um = zeros(num_cells, 3);

        if COMPUTE_CELL_CENTROIDS
            for currentCell = cell_ids
                disp(['--Computing centroids for cell ' num2str(currentCell) ' out of ' num2str(num_cells)]);
                [centroid_voxels, centroid_um] = compute_centroid(individual_cell_masks{currentCell}, VOXEL_SIZE_XY, VOXEL_SIZE_Z);
                cell_centroids_voxels(currentCell, :) = centroid_voxels;
                cell_centroids_um(currentCell, :) = centroid_um;            
            end

            T_centroids_table = table(cell_ids', cell_centroids_voxels(:, 1), cell_centroids_voxels(:, 2), cell_centroids_voxels(:, 3), cell_centroids_um(:,1), cell_centroids_um(:,2), cell_centroids_um(:,3));
            T_centroids_table.Properties.VariableNames = {'cell_id_dup', 'cell_centroid_voxel_x', 'cell_centroid_voxel_y', 'cell_centroid_voxel_z', 'cell_centroid_um_x', 'cell_centroid_um_y', 'cell_centroid_um_z'};

            writetable(T_centroids_table, [output_directory metadata '_centroid.csv'], 'WriteVariableNames', true);
        end

        %T_centroids_table.Properties.VariableNames{1} = strcat(T_centroids_table.Properties.VariableNames{1}, '_dup');

        
         % VOLUME DIMENSIONS
        if COMPUTE_DIMENSIONS
            %disp(['--Writing volume dimensions for embryo ' num2str(currentMask) ' out of ' num2str(num_embryos)]);
            T_dimensions = array2table([xdim, ydim, zdim, num_cells]);
            T_dimensions.Properties.VariableNames = {'embryo_x_voxels', 'embryo_y_voxels', 'embryo_z_voxels', 'embryo_cell_count'};
            writetable(T_dimensions, [output_directory metadata '_dimensions.csv'], 'WriteVariableNames', true);
        end   
        
        T_dimensions_table = repmat(T_dimensions, num_cells, 1);
        % COMPUTE GLOBAL EMBRYO CENTROID

        if COMPUTE_EMBRYO_CENTROIDS                
            tic                

            %disp(['-Saving embryo centroids and radii for embryo ' num2str(currentMask)]);
            [embryo_centroid_voxels, embryo_centroid_um, embryo_radius_um] = compute_embryo_centroid(segmentation_mask, SURFACE_AREA_CONNECTIVITY, VOXEL_SIZE_XY, VOXEL_SIZE_Z);

            T_centroid = [array2table(embryo_centroid_voxels), array2table(embryo_centroid_um)];
            T_centroid.Properties.VariableNames = {'embryo_centroid_voxel_x', 'embryo_centroid_voxel_y', 'embryo_centroid_voxel_z', 'embryo_centroid_um_x', 'embryo_centroid_um_y', 'embryo_centroid_um_z'};

            writetable(T_centroid, [output_directory metadata '_global_centroid.csv'], 'WriteRowNames', false, 'WriteVariableNames', true);

            T_radius = array2table(embryo_radius_um);
            T_radius.Properties.VariableNames = {'radius'};

            writetable(T_radius, [output_directory metadata '_radius.csv'], 'WriteRowNames', false, 'WriteVariableNames', true);

            toc
        end
        
        T_centroid_embryo_table = repmat(T_centroid, num_cells, 1);
        T_radius_table = repmat(T_radius, num_cells, 1);
        
        Table = [T_volumes_table T_centroids_table T_dimensions_table T_centroid_embryo_table T_radius_table];
        
        Table.cell_id = cellstr(strcat(sprintf('%s_',p_w_e), string(Table.cell_id)));


        %disp(['Completed embryo ' num2str(currentMask) ' out of ' num2str(num_embryos) ' at ' char(datetime)]);
        
        if counter == 1
            embryo_geometry_table = Table;
            counter = counter+1;
        else
            embryo_geometry_table = [embryo_geometry_table; Table];
        end
end


distance = table('Size', [0, 1], 'VariableTypes', {'double'}, 'VariableNames', {'distance_to_center'});

numRows = height(embryo_geometry_table);

for i = 1:numRows
    % Define the coordinates of the center of the object
    x1 = table2array(embryo_geometry_table(i,'embryo_centroid_voxel_x'));
    y1 = table2array(embryo_geometry_table(i,'embryo_centroid_voxel_y'));
    z1 = table2array(embryo_geometry_table(i,'embryo_centroid_voxel_z'));

    center = [x1, y1, z1];
    % Define the coordinates of the dot in space
    x2 = table2array(embryo_geometry_table(i,'cell_centroid_voxel_x'));
    y2 = table2array(embryo_geometry_table(i,'cell_centroid_voxel_y'));
    z2 = table2array(embryo_geometry_table(i,'cell_centroid_voxel_z'));    
    dot = [x2, y2, z2];

    % Calculate the Euclidean distance
    d = sqrt(sum((center - dot).^2));
    distance(i, 1) = array2table(d);
end

embryo_geometry_table = [embryo_geometry_table distance];

%writetable(embryo_geometry_table, [output_directory 'embryo_geometry_table_onlyMorula.csv'], 'WriteRowNames', false, 'WriteVariableNames', true);
writetable(embryo_geometry_table, [output_directory 'embryo_geometry_table_blast.csv'], 'WriteRowNames', false, 'WriteVariableNames', true);


close all; clear;




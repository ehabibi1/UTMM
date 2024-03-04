function [cellxgene] = alignment_pipeline_mechanotranscriptome(home_dir, species, expt, exptxlsfile,sampleId, viz, doDisplacement, dropUncertainMeasures, round)
%alignment_pipeline Align the embryo mechanics
%   projectname - string, e.g. 'mechano_transcriptome_project'
%   experimentname - string, e.g. 'joint_measurement_may_2021'
%   exptxlsfile - path to excel spreadsheet (first sheet should be
%   called 'Sheet1'), e.g. 'Embryo Samples Mechanics Metadata.xlsx'
%   sampleId - integer, row of sample in exptxlsfile, starting from first
%   row after header in excel file, (i.e. header is sampleId=0)
%   showRegistration - bool, e.g. true, if set to true will show registered
%   live (mechanical) dapi image overlayed on top of membrane; if false
%   will show fixed morpho image dapi nuclei.
%   doDisplacement - bool, e.g. true, whether or not to apply x,y,z
%   displacement field to live image transformation.
% set(gcf,'Position', [1639 1107 946 933])
% if strcmp(getenv('HOSTNAME'),'pyadolla-centos-vm-2'),
%     wd = fullfile('/mnt/disks/datadisk1/gsbucket/ehsan/', expt);
%     xlspath = fullfile(wd,exptxlsfile);
% end

if strcmp(getenv('COMPUTERNAME'),'GPDC2-CF0'),
    wd = fullfile(sprintf('%s/%s/%s/', home_dir, species, expt));
    xlspath = fullfile(wd,sprintf('mechanics/%s',exptxlsfile));
else,
    [ret, hostname] = system('hostname')
    wd = fullfile(sprintf('%s/%s/%s/', home_dir, species, expt));
    xlspath = fullfile(wd,sprintf('mechanics/%s',exptxlsfile));
end

xldat = readtable(xlspath);

plate = xldat.Plate(sampleId);
well = xldat.Well(sampleId);
embryo = xldat.EmbryoNum(sampleId);


liveimagefname = xldat.LiveFileName{sampleId};
livenucleifname = xldat.LiveProcessedFileName{sampleId};
fixednucleifname = xldat.FixedRegistrationNucleiImagesFileName{sampleId};
fixedmembranefname = xldat.FixedRegistrationMembraneImagesFileName{sampleId};
%nuclei_segmentation_fname = xldat.NucleiSegmentationFileName{sampleId};
membrane_segmentation_fname = xldat.MembraneSegmentationFileName{sampleId};
cellcount_table_fname = xldat.CellCountTableFileName{sampleId};
fixedsamplename = xldat.FixedSampleName{sampleId};
iscontrastadjusted = xldat.Contrast_adjusted{sampleId};
mechxlsfname = xldat.ResultMechanicalMeasurementFile{sampleId};

results_path = fullfile(wd, 'joint_analysis');
registration_path = fullfile(results_path,'live-to-fixed-registration',fixedsamplename);
affine_reg_path = fullfile(registration_path,'Exported_data');
displacement_field_path = fullfile(registration_path,'displacement_registration','Exported_data');

affine_reg_fpath = fullfile(affine_reg_path,'transform_global_img_moving.txt');
displacement_field_fpath_x = fullfile(displacement_field_path, 'transform_global_img_moving.x.tif');
displacement_field_fpath_y = fullfile(displacement_field_path, 'transform_global_img_moving.y.tif');
displacement_field_fpath_z = fullfile(displacement_field_path, 'transform_global_img_moving.z.tif');



%experiments_path = fullfile(wd,'mechanics',experimentname);
%fixednuclei_img_fpath = fullfile(experiments_path, 'fixed_morphology_images','staining',fixednucleifname);
fixednuclei_img_fpath = fullfile(wd, 'joint_analysis','fixed_morphology_images',fixednucleifname);
livenuclei_img_fpath = fullfile(wd, 'joint_analysis','live_morphology_images',livenucleifname);


ReReg_dir = sprintf('%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/ReReg',wd,plate, well, embryo);
   
    if exist(ReReg_dir, 'dir')
        time = datestr(datetime(now,'ConvertFrom','datenum'));
        disp(sprintf('%s ----- Loading  segmentation mask from ReReg dir ...', time));
        fixedmembrane_img_fpath  = sprintf('%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/ReReg/%s', wd, plate, well, embryo,fixedmembranefname);

        membrane_segmentation_fpath = sprintf('%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/ReReg/%s', wd, plate, well, embryo,membrane_segmentation_fname);
    else
        time = datestr(datetime(now,'ConvertFrom','datenum'));
        disp(sprintf('%s ----- Loading  segmentation mask  ...', time));
        fixedmembrane_img_fpath  = sprintf('%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/%s', wd, plate, well, embryo,fixedmembranefname);        
        membrane_segmentation_fpath = sprintf('%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/%s', wd, plate, well, embryo,membrane_segmentation_fname);
    end 
    

%nuclei_segmentation_fpath  = fullfile(experimetns_path,'fixed_morphology_images','segmented_sembryo','tiff','nuclei',nuclei_segmentation_fname);
%mechanical_measurements_xls_fpath = fullfile(wd, 'joint_analysis','mechanical_measurements',mechxlsfname);
mechanical_measurements_xls_fpath = sprintf('%s/joint_analysis/mechanical_measurements/round%s/%s',wd,round,mechxlsfname);

cellxgenexmech_dir = sprintf('%s/joint_analysis/cellxgenexmech',wd);
if ~exist(cellxgenexmech_dir, 'dir') mkdir(cellxgenexmech_dir), end
%[~,samplename,~]=fileparts(membrane_segmentation_fname)
if ~dropUncertainMeasures
    %cellcounts_table_output_fpath = fullfile(results_path,'expression_measurements', ['cell_counts_table_' samplename '_with_msd.txt']);
    cellcounts_table_output_fpath = sprintf('%s/cell_counts_table_plate_%d_well%d_embryo_%d_with_msd_round%s.txt', cellxgenexmech_dir, plate, well, embryo, round);
    
else
    %cellcounts_table_output_fpath = fullfile(results_path,'expression_measurements_drop_uncertain', ['cell_counts_table_' samplename '_with_msd_drop_unc.txt']);
    cellcounts_table_output_fpath = sprintf('%s/cell_counts_table_plate_%d_well%d_embryo_%d_with_msd_round%s_drop_unc.txt', cellxgenexmech_dir, plate, well, embryo, round);
end

cellcounts_table_fpath =  sprintf('%s/transcriptome/plate_%d_well_%d/embryo_%d/figure/cell_counts_table/%s', wd, plate, well, embryo,cellcount_table_fname);

[~,livename,~]=fileparts(liveimagefname);

%round = xldat.ResultMechanicalMeasurementRound(sampleId);
live_img_series_info_fpath =  sprintf('%s/mechannics/experiment/mechanical_measurements/round%s/results/%s.xlsx', wd, round,livename);
%live_img_series_info_fpath = fullfile(results_path,'mechanical_measurements','mouse_embryo_20210505','live_image_series_infos',[livename '.xlsx']);

%missing_cell_mechano_measurement_table_output_fpath = fullfile(experiments_path,'mechanical_measurements','mouse_embryo_20210505','missing_mechanical_measurements',[livename ' xyz.xlsx']);
missing_cell_mechano_measurement_table_output_fpath = fullfile(wd, 'joint_analysis','missing_mechanical_measurements',mechxlsfname);

%javaaddpath '/opt/MATLAB/R2022a/java/jar/mij.jar'
%javaaddpath '/opt/MATLAB/R2022a/java/jar/ij.jar'

%addpath(genpath('/home/pyadolla/Fiji.app/scripts/'));
% addpath(genpath('\ssodium\broad_thechenlab\ehsan\analysis\InSitu_preImpEmbryo\mechanics_transcriptome_alignment_pipeline'));

%ImageJ 2021
MIJ.start %2017
%IJ=ij.IJ(); 


livenuclei_imp = ij.IJ.openImage(livenuclei_img_fpath);
livenuclei_cal = livenuclei_imp.getCalibration();
livenuclei_pxWidth = livenuclei_cal.pixelWidth();
livenuclei_pxHeight = livenuclei_cal.pixelHeight();
livenuclei_pxDepth = livenuclei_cal.pixelDepth();
livenuclei_width = livenuclei_imp.getWidth();
livenuclei_height = livenuclei_imp.getHeight();
livenuclei_depth = livenuclei_imp.getNSlices();
livenuclei_imp.close();
Rlive = imref3d([livenuclei_height livenuclei_width livenuclei_depth], livenuclei_pxWidth, livenuclei_pxHeight, livenuclei_pxDepth);
livenuclei = FastTiff(livenuclei_img_fpath);

fixednuclei_imp = ij.IJ.openImage(fixednuclei_img_fpath);
fixednuclei_cal = fixednuclei_imp.getCalibration();
fixednuclei_pxWidth = fixednuclei_cal.pixelWidth();
fixednuclei_pxHeight = fixednuclei_cal.pixelHeight();
fixednuclei_pxDepth = fixednuclei_cal.pixelDepth();
fixednuclei_width = fixednuclei_imp.getWidth();
fixednuclei_height = fixednuclei_imp.getHeight();
fixednuclei_depth  = fixednuclei_imp.getNSlices();
fixednuclei_imp.close();
%Rfixed = imref3d([fixednuclei_height fixednuclei_width fixednuclei_depth],fixednuclei_pxWidth,fixednuclei_pxHeight,fixednuclei_pxDepth);
%Manually set pxWidth, pxHeight, pxDepth, because python script for
%saving fixed image tiffs did not set calibration information in tif
%files.
Rfixed = imref3d([fixednuclei_height fixednuclei_width fixednuclei_depth],0.1625000,0.1625000,0.3986799);
% Rfixed = imref3d([fixednuclei_height fixednuclei_width fixednuclei_depth],fixednuclei_pxWidth, fixednuclei_pxHeight, fixednuclei_pxDepth);
fixednuclei = FastTiff(fixednuclei_img_fpath);
fixedmembrane = FastTiff(fixedmembrane_img_fpath);


membrane_segmentation = FastTiff(membrane_segmentation_fpath);
membrane_segmentation = imresize3(membrane_segmentation,0.5,'nearest');

transforms = cli_lineartransformread(affine_reg_fpath,'itk');
affine_transform = inv(transforms)';

if (doDisplacement)
    dx = FastTiff(displacement_field_fpath_x);
    dy = FastTiff(displacement_field_fpath_y);
    dz = FastTiff(displacement_field_fpath_z);
end

mechano_table = readtable(mechanical_measurements_xls_fpath);
mechano_xyz = table2array(mechano_table(:,{'x','y','z'}));
mechano_msd = mechano_table.MSD;

mechano_xyz_tr = transformPointsForward(affine3d(affine_transform), mechano_xyz);

cellxgene = readtable(cellcounts_table_fpath);
cellxgene=cellxgene(~isnan(cellxgene.x),:);
cellxgene.cell_index = [1:size(cellxgene,1)]';

[Xf,Yf,Zf]=meshgrid(linspace(Rfixed.XWorldLimits(1),Rfixed.XWorldLimits(2),size(fixednuclei,2)),...
                        linspace(Rfixed.YWorldLimits(1),Rfixed.YWorldLimits(2),size(fixednuclei,1)),...
                        linspace(Rfixed.ZWorldLimits(1),Rfixed.ZWorldLimits(2),size(fixednuclei,3)));

vizshowmemb = bitor(VizOpts.Visualize, VizOpts.Membrane);
if eq(vizshowmemb, (bitand(vizshowmemb, viz)))
    h_fixedmemb = showEmbryo(fixedmembrane,Xf,Yf,Zf,Rfixed,100,150,0.05);
    hold on;

    colors=colormap(jet(16));
end


if exist(missing_cell_mechano_measurement_table_output_fpath)
    missing_img_series_infos = readtable(missing_cell_mechano_measurement_table_output_fpath);
    if isempty(missing_img_series_infos)
        disp('no entries in missing mechano measurement file, so not using any missing measurements.');
    else
        offset_series = mechano_table.Filename{1};
        missing_img_series_infos.Properties.RowNames = missing_img_series_infos.Filename;
        %missing_img_series_infos.Filename=[];

        live_img_series_infos = readtable(live_img_series_info_fpath);
        live_img_series_infos.Properties.RowNames = live_img_series_infos.Filename;
        live_img_series_infos.Filename=[];

        mechano_table.Properties.RowNames = mechano_table.Filename;
        %mechano_table.Filename=[];
        %subtract offset series
        missing_img_series_infos.x = missing_img_series_infos.x - live_img_series_infos.x(offset_series);
        missing_img_series_infos.y = missing_img_series_infos.y - live_img_series_infos.y(offset_series);
        missing_img_series_infos.z = missing_img_series_infos.z - live_img_series_infos.z(offset_series);
        missing_img_series_infos(intersect(missing_img_series_infos.Properties.RowNames,mechano_table.Properties.RowNames),:)=[]; %remove those series in mechano_table
        %add mechano_table values for offset series
        missing_img_series_infos.x = missing_img_series_infos.x + mechano_table.x(offset_series);
        missing_img_series_infos.y = missing_img_series_infos.y + mechano_table.y(offset_series);
        missing_img_series_infos.z = missing_img_series_infos.z + mechano_table.z(offset_series);
        missing_img_series_xyz = table2array(missing_img_series_infos(:,{'x','y','z'}));
        missing_img_series_msd = table2array(missing_img_series_infos(:,{'MSD'}));
        missing_img_series_xyz_tr = transformPointsForward(affine3d(affine_transform), missing_img_series_xyz);
        mechano_xyz = [mechano_xyz; missing_img_series_xyz];
        mechano_xyz_tr = [mechano_xyz_tr; missing_img_series_xyz_tr];
        mechano_msd = [mechano_msd; missing_img_series_msd];
        loc_offset = size(mechano_table,1);
        mechano_table = [mechano_table; missing_img_series_infos];
    end
end
cellj = setdiff(unique(membrane_segmentation(:)),0);
oldindices = cellj;
tmp = zeros(size(membrane_segmentation));
for i=1:length(cellj)
    tmp(membrane_segmentation==cellj(i))=i;
end
membrane_segmentation = tmp;
cellj = setdiff(unique(membrane_segmentation(:)),0);

hold on;
centers=[];
for j = 1:length(cellj)
   labelidx = find(membrane_segmentation == cellj(j));
   temp = zeros(size(membrane_segmentation));
   temp(labelidx) = 1;
   temp = imerode(temp, true(7));

   [nX,nY,nZ] = size(temp);
   [X,Y,Z] = meshgrid(1:nY,1:nX,1:nZ);
   [f v] = isosurface(X,Y,Z,temp,0.5);
   v=v.*repmat([Rfixed.PixelExtentInWorldX*2, ...
                Rfixed.PixelExtentInWorldY*2, ...
                Rfixed.PixelExtentInWorldZ*2],size(v,1),1);
   %v(:,1) = RfixedImageExtentInWorldX - v(:,1);
   %v(:,3) = Rfixed.ImageExtentInWorldZ - v(:,3);
   p = patch('Faces',f,'Vertices',v,'EdgeColor','none','FaceColor','g','FaceAlpha',0.1, 'EdgeAlpha',0); %0.3
   reducepatch(p,40);

   cinds=unique(convhull(v,'simplify',true));
   centers = [centers; centroid(v(cinds,:))];
   [distances, surface_points]=point2trimesh('Faces',f,'Vertices',v,'QueryPoints',mechano_xyz_tr);

   %dt = delaunayTriangulation(v);
   %ti = pointLocation(dt,mechano_xyz_tr);
   %locs=find(~isnan(ti));
   locs = find(distances > 0);
   cell2measurement{j} = locs;
end
for i=1:size(mechano_xyz_tr,1)
    text(mechano_xyz_tr(i,1),mechano_xyz_tr(i,2),mechano_xyz_tr(i,3),sprintf('%d',i));
    scatter3(mechano_xyz_tr(i,1),mechano_xyz_tr(i,2),mechano_xyz_tr(i,3),500, 'filled','MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.3,'MarkerFaceColor','r');%colors(i,:));
end

if ~dropUncertainMeasures
    points_outside_cells = setdiff(1:size(mechano_xyz_tr,1),sort(cat(1,cell2measurement{:})));
    closestcell_inds=zeros(1,length(points_outside_cells));
%     x=cellxgene{:,find(cellxgene.Properties.VariableNames == "x")};
%     y=cellxgene{:,find(cellxgene.Properties.VariableNames == "y")};
%     z=cellxgene{:,find(cellxgene.Properties.VariableNames == "z")};
%     P = [y x z];
%     for i = 1:length(points_outside_cells)
%         loc = points_outside_cells(i);
%         qp = mechano_xyz_tr(loc,:);
%         closestcell_inds(i) = dsearchn(P,qp);
%         cell2measurement{1,closestcell_inds(i)} = [cell2measurement{1,closestcell_inds(i)}; loc];
%     end
    %find assignment based on closest distance of cell membrane to patch centers (points)
    %for i = 1:length(points_outside_cells)
    %    distances=[];
    %    loc = points_outside_cells(i);
    %    for j = 1:length(cellj)
    %        labelidx = find(membrane_segmentation == cellj(j));
    %        temp = zeros(size(membrane_segmentation));
    %        temp(labelidx) = 1;
    %        temp = imerode(temp, true(7));

    %        [nX,nY,nZ] = size(temp);
    %        [X,Y,Z] = meshgrid(1:nY,1:nX,1:nZ);
    %        [f v] = isosurface(X,Y,Z,temp,0.5);
    %        v=v.*repmat([Rfixed.PixelExtentInWorldX*2, ...
    %                     Rfixed.PixelExtentInWorldY*2, ...
    %                     Rfixed.PixelExtentInWorldZ*2],size(v,1),1);
    %        qp = mechano_xyz_tr(loc,:);
    %        [distance, surface_point]=point2trimesh('Faces',f,'Vertices',v,'QueryPoints',qp);
    %        distances = [distances distance];
    %    end
    %    distances = abs(distances);
    %    closestcell_inds(i) = find(distances==min(distances));
    %    cell2measurement{1,closestcell_inds(i)} = [cell2measurement{1,closestcell_inds(i)}; loc];
    %end

    a3d=affine3d(affine_transform);
    wbf=waitbar(0,'Assigning remaining patches...');
    for i = 1:length(points_outside_cells)
      loc = points_outside_cells(i);
      overlaps=[];
      for j = 1:length(cellj)
        labelidx = find(membrane_segmentation == cellj(j));
        temp = zeros(size(membrane_segmentation));
        temp(labelidx) = 1;
        temp = imerode(temp, true(7));

        [nX,nY,nZ] = size(temp);
        [X,Y,Z] = meshgrid(1:nY,1:nX,1:nZ);
        [f v] = isosurface(X,Y,Z,temp,0.5);
        v=v.*repmat([Rfixed.PixelExtentInWorldX*2, ...
                     Rfixed.PixelExtentInWorldY*2, ...
                     Rfixed.PixelExtentInWorldZ*2],size(v,1),1);
        dt = delaunayTriangulation(v);

        waitbar(((i-1)*length(cellj) + (j-1))/(length(cellj)*length(points_outside_cells)),wbf,'Assigning remaining patches...');
        v=[];
        v0=[];
        for xd=[-32 32]
            for yd=[-32 32]
                pt = mechano_xyz(loc,:) + [xd*0.3661187 yd*0.3661187  0];%[xd*Rlive.PixelExtentInWorldX yd*Rlive.PixelExtentInWorldY 0];
                % pt = mechano_xyz(loc,:) + [xd*Rlive.PixelExtentInWorldX yd*Rlive.PixelExtentInWorldY 0];
                v = [v; transformPointsForward(a3d, pt)];
                v0 = [v0; pt];
            end
        end
        [x,y]=meshgrid(linspace(v0(1,1),v0(3,1),25),linspace(v0(1,2),v0(2,2),25));
        z = ones(size(x))*v0(1,3);
        [xt,yt,zt]=transformPointsForward(a3d,x(:),y(:),z(:));
        xt = reshape(xt,[],1);
        yt = reshape(yt,[],1);
        zt = reshape(zt,[],1);
        ti = pointLocation(dt,[xt yt zt]);
        pct_cvr=sum(~isnan(ti))/size(xt,1);
        overlaps =[overlaps; pct_cvr];
      end
      [max_pct_cvr, max_ind] =max(overlaps);
      if max_pct_cvr > 0.25 % >25% overlap of patch with cell
        cell2measurement{cellj(max_ind)} = [cell2measurement{cellj(max_ind)}; loc];
      end
    end
    close(wbf);
    % now add overlapping patches if they exist
end


if dropUncertainMeasures
    cell2msd = nan(1,numel(vertcat(cell2measurement{:})));
    count=1;
    tmpcellxgene = cell2table(cell(0, length(cellxgene.Properties.VariableNames)),'VariableNames',cellxgene.Properties.VariableNames);
    tmpoldindices=[];
    for i = 1:length(cell2measurement)
        locs = cell2measurement{i};
        if ~isempty(locs)
            for j=1:length(locs)
                tmpcellxgene = [tmpcellxgene; cellxgene(i,:)];
                cell2msd(count) = mechano_msd(locs(j));
                tmpoldindices =[tmpoldindices; oldindices(i)];
                count = count+1;
            end
        else
            tmpcellxgene = [tmpcellxgene; cellxgene(i,:)];
            cell2msd(count)=nan;
            tmpoldindices = [tmpoldindices; oldindices(i)];
            count = count+1;
        end
    end
    tmpcellxgene.msd = cell2msd';
    tmpcellxgene.oldindices = tmpoldindices;
    cellxgene = tmpcellxgene;
else
    cell2msd = nan(1,length(cell2measurement));
    cell2msd_sd = nan(1,length(cell2measurement));
    for i=1:length(cell2measurement)
        locs = cell2measurement{i};
        if ~isempty(locs)
            cell2msd(i) = mean(mechano_msd(locs));
            cell2msd_sd(i) = std(mechano_msd(locs));
        end
    end
    cellxgene.msd = cell2msd';
    cellxgene.msd_sd = cell2msd_sd';
    cellxgene.oldindics = oldindices;
end

if eq(VizOpts.Visualize, (bitand(VizOpts.Visualize, viz)))

    if eq(VizOpts.Registered, (bitand(VizOpts.Registered, viz)))
       livenuclei_post_global_t = imwarp(livenuclei,Rlive,affine3d(affine_transform),'OutputView',Rfixed);
       %livenuclei_post_dense_t = imwarp(livenuclei_post_global_t,cat(4,dx/Rfixed.PixelExtentInWorldX,dy/Rfixed.PixelExtentInWorldY,dz/Rfixed.PixelExtentInWorldZ));
       bscale =max(fixedmembrane(:))/max(livenuclei_post_global_t(:));
       vizshowregpatch = bitor(VizOpts.Registered, VizOpts.PatchBounds);
       if eq(vizshowregpatch, bitand(vizshowregpatch, bitxor(viz,VizOpts.PatchBounds)))
           hold on; h_nuclei = showEmbryo(bscale*livenuclei_post_global_t,Xf,Yf,Zf,Rfixed,100,bscale*25,0.2); hold off;% alpha(0.2);
           disp('showing registered nuclei');
       end
    elseif eq(VizOpts.Fixed, (bitand(VizOpts.Fixed, viz)))
       bscale = max(fixedmembrane(:))/max(fixednuclei(:));
       hold on; h_nuclei = showEmbryo(bscale*fixednuclei,Xf,Yf,Zf,Rfixed,100,bscale*25,0.2); hold off;
       disp('showing fixed nuclei');
    else
       disp('not showing nuclei');
    end

    hold on;
    for i=1:size(mechano_xyz_tr,1)
        %text(M_tr(i,1),M_tr(i,2),M_tr(i,3),sprintf('[%.1f, %.1f, %.1f]',mechano_xyz_tr(i,1),mechano_xyz_tr(i,2),mechano_xyz_tr(i,3)));
        scatter3(mechano_xyz_tr(i,1),mechano_xyz_tr(i,2),mechano_xyz_tr(i,3),500, 'filled','MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.3,'MarkerFaceColor','r');%colors(i,:));
    end

    vizassignment = bitor(VizOpts.CHullCenter, VizOpts.Assignment);
    if eq(vizassignment, (bitand(vizassignment, viz)))
        for i=1:size(centers,1)
            for j=1:length(cell2measurement{i})
                mind=cell2measurement{i}(j);
                plot3([mechano_xyz_tr(mind,1); centers(i,1)],...
                      [mechano_xyz_tr(mind,2); centers(i,2)],...
                      [mechano_xyz_tr(mind,3); centers(i,3)], 'r--');
            end
        end
    end
    vizchullcenter = bitor(VizOpts.CHullCenter, VizOpts.SegCenter);
    if eq(vizchullcenter, bitand(vizchullcenter, bitxor(viz,VizOpts.SegCenter)))
        for i=1:size(centers,1)
            scatter3(centers(i,1), centers(i,2), centers(i,3), 500, 'filled','MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.3,'MarkerFaceColor','c');
            text(centers(i,1), centers(i,2), centers(i,3),num2str(i));
        end
    end
    vizsegcenter = bitor(VizOpts.SegCenter, VizOpts.CHullCenter);
    if eq(vizsegcenter, bitand(vizsegcenter, bitxor(viz,VizOpts.CHullCenter)))
         for i=1:size(centers,1)
            x=cellxgene{i,find(cellxgene.Properties.VariableNames == "x")};
            y=cellxgene{i,find(cellxgene.Properties.VariableNames == "y")};
            z=cellxgene{i,find(cellxgene.Properties.VariableNames == "z")};
            scatter3(y,x,z,500, 'filled','MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.3,'MarkerFaceColor','y');
         end
    end
    if eq(VizOpts.PatchBounds, bitand(VizOpts.PatchBounds, viz))
        hold on;
        a3d=affine3d(affine_transform);
        for i=1:size(mechano_xyz,1)
            v=[];
            v0=[];
            for xd=[-32 32]
                for yd=[-32 32]
                    pt = mechano_xyz(i,:) + [xd*0.3661187 yd*0.3661187  0];%[xd*Rlive.PixelExtentInWorldX yd*Rlive.PixelExtentInWorldY 0];
                    % pt = mechano_xyz(i,:) + [xd*Rlive.PixelExtentInWorldX yd*Rlive.PixelExtentInWorldY 0];
                    v = [v; transformPointsForward(a3d, pt)];
                    v0 = [v0; pt];
                end
            end
            if eq(VizOpts.Registered, bitand(VizOpts.Registered, viz))
                [x,y]=meshgrid(linspace(v0(1,1),v0(3,1),50),linspace(v0(1,2),v0(2,2),50));
                z = ones(size(x))*v0(1,3);
                [xt,yt,zt]=transformPointsForward(a3d,x(:),y(:),z(:));
                xt = reshape(xt,size(x));
                yt = reshape(yt,size(y));
                zt = reshape(zt,size(z));
                %disp(norm(v(3,:)-v(1,:)))
                %disp(norm(v(2,:)-v(1,:)))
                ph=slice(Xf,Yf,Zf,bscale*livenuclei_post_global_t,xt,yt,zt);
                set(ph, 'FaceColor', 'interp', 'EdgeColor', 'none','FaceAlpha',0.3, 'EdgeAlpha',0.3);
            else
                patch('Faces',[1 2 4 3], 'Vertices', v, 'FaceColor','cyan','FaceAlpha',0.3, 'EdgeColor','none');
            end
            line(v([1 3],1), v([1 3],2), v([1 3], 3),'Color','red','LineWidth',2);
            line(v([1 2],1), v([1 2],2), v([1 2], 3), 'Color','blue','LineWidth',2);
            text(mechano_xyz_tr(i,1),mechano_xyz_tr(i,2),mechano_xyz_tr(i,3),mechano_table.Filename{i},'FontSize',8,'Color','red','VerticalAlignment','top');
        end
    end
    vizmsd = bitor(VizOpts.MSD, VizOpts.CHullCenter);
    if eq(vizmsd, bitand(vizmsd, bitxor(viz,VizOpts.CHullCenter)))
        colors=colormap(autumn(100));
        msd_col = round(cellxgene.msd*(100/max(cellxgene.msd)));
        %figure
        hold on
        for i=1:length(msd_col)
            if ~isnan(msd_col(i))
                x=cellxgene{i,find(cellxgene.Properties.VariableNames == "x")};
                y=cellxgene{i,find(cellxgene.Properties.VariableNames == "y")};
                z=cellxgene{i,find(cellxgene.Properties.VariableNames == "z")};
                scatter3(y,x,z,500, 'filled','MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.3,'MarkerFaceColor',colors(msd_col(i),:));
                text(y,x,z,num2str(i));
            end
        end
    end
%     if eq(VizOpts.PickPatchesForMissing, bitand(VizOpts.PickPatchesForMissing, viz))
%         live_img_series_infos = readtable(live_img_series_info_fpath);
%         offset_series = mechano_table.Filename{1};
%         live_img_series_infos.Properties.RowNames = live_img_series_infos.Filename;
%         live_img_series_infos.Filename=[];
%         mechano_table.Properties.RowNames = mechano_table.Filename;
%         mechano_table.Filename=[];
%         %subtract offset series
%         live_img_series_infos.x = live_img_series_infos.x - live_img_series_infos.x(offset_series);
%         live_img_series_infos.y = live_img_series_infos.y - live_img_series_infos.y(offset_series);
%         live_img_series_infos.z = live_img_series_infos.z - live_img_series_infos.z(offset_series);
%         live_img_series_infos(mechano_table.Properties.RowNames,:)=[]; %remove those series in mechano_table
%         %add mechano_table values for offset series
%         live_img_series_infos.x = live_img_series_infos.x + mechano_table.x(offset_series);
%         live_img_series_infos.y = live_img_series_infos.y + mechano_table.y(offset_series);
%         live_img_series_infos.z = live_img_series_infos.z + mechano_table.z(offset_series);
%         live_img_series_xyz = table2array(live_img_series_infos(:,{'x','y','z'}));
% 
%         live_img_series_xyz_tr = transformPointsForward(affine3d(affine_transform), live_img_series_xyz);
%         a3d=affine3d(affine_transform);
%         %non-covered cells
%         empty_cells = find(cellfun(@isempty,cell2measurement));
% 
%         wbf=waitbar(0,'Assigning patches...');
%         for j = 1:length(empty_cells)
%             labelidx = find(membrane_segmentation == empty_cells(j));
%             temp = zeros(size(membrane_segmentation));
%             temp(labelidx) = 1;
%             temp = imerode(temp, true(7));
% 
%             [nX,nY,nZ] = size(temp);
%             [X,Y,Z] = meshgrid(1:nY,1:nX,1:nZ);
%             [f v] = isosurface(X,Y,Z,temp,0.5);
%             v=v.*repmat([Rfixed.PixelExtentInWorldX*2, ...
%                          Rfixed.PixelExtentInWorldY*2, ...
%                          Rfixed.PixelExtentInWorldZ*2],size(v,1),1);
%             dt = delaunayTriangulation(v);
% %             [distances, surface_points]=point2trimesh('Faces',f,'Vertices',v,'QueryPoints',live_img_series_xyz_tr);
% %             %dt = delaunayTriangulation(v);
% %             %ti = pointLocation(dt,live_img_series_xyz_tr);
% %             %locs=find(~isnan(ti));
% %             locs = find(distances > 0);
% %             if isempty(locs)
% %                 distances = abs(distances);
% %                 locs = find(distances==min(distances));
% %             end
%             locs=[];
%             for i=1:size(live_img_series_xyz,1)
%                 waitbar(((j-1)*size(live_img_series_xyz,1)+(i-1))/(size(live_img_series_xyz,1)*length(empty_cells)),wbf,'Assigning patches...');
%                 v=[];
%                 v0=[];
%                 for xd=[-32 32]
%                     for yd=[-32 32]
%                         pt = live_img_series_xyz(i,:) + [xd*0.3661187 yd*0.3661187  0];%[xd*Rlive.PixelExtentInWorldX yd*Rlive.PixelExtentInWorldY 0];
%                         v = [v; transformPointsForward(a3d, pt)];
%                         v0 = [v0; pt];
%                     end
%                 end
%                 [x,y]=meshgrid(linspace(v0(1,1),v0(3,1),25),linspace(v0(1,2),v0(2,2),25));
%                 z = ones(size(x))*v0(1,3);
%                 [xt,yt,zt]=transformPointsForward(a3d,x(:),y(:),z(:));
%                 xt = reshape(xt,[],1);
%                 yt = reshape(yt,[],1);
%                 zt = reshape(zt,[],1);
%                 ti = pointLocation(dt,[xt yt zt]);
%                 pct_cvr=sum(~isnan(ti))/size(xt,1);
%                 if pct_cvr > 0.25 % >25% overlap of patch with cell
%                     locs = [locs; i];
%                 end
%             end
%             if isempty(locs)
%                 disp(['Found no patches in lif file to assign to cell ' num2str(empty_cells(j)) '!']);
%             end
%             cell2measurement{empty_cells(j)} = locs;
% 
%             for i=1:length(locs)
%                 scatter3(live_img_series_xyz_tr(locs(i),1),live_img_series_xyz_tr(locs(i),2),live_img_series_xyz_tr(locs(i),3),500, 'filled','MarkerEdgeAlpha',0.3,'MarkerFaceAlpha',0.3,'MarkerFaceColor','g');
%                 text(live_img_series_xyz_tr(locs(i),1),live_img_series_xyz_tr(locs(i),2),live_img_series_xyz_tr(locs(i),3),live_img_series_infos.Properties.RowNames{locs(i)});
%                 if eq(vizassignment, (bitand(vizassignment, viz)))
%                     plot3([live_img_series_xyz_tr(locs(i),1); centers(empty_cells(j),1)],...
%                     [live_img_series_xyz_tr(locs(i),2); centers(empty_cells(j),2)],...
%                     [live_img_series_xyz_tr(locs(i),3); centers(empty_cells(j),3)], 'g--');
%                 end
%                 if eq(VizOpts.PatchBounds, bitand(VizOpts.PatchBounds, viz))
%                     a3d=affine3d(affine_transform);
%                     v=[];
%                     v0=[];
%                     for xd=[-32 32]
%                         for yd=[-32 32]
%                             pt = live_img_series_xyz(locs(i),:) + [xd*0.3661187 yd*0.3661187  0];%[xd*Rlive.PixelExtentInWorldX yd*Rlive.PixelExtentInWorldY 0];
%                             v = [v; transformPointsForward(a3d, pt)];
%                             v0 = [v0; pt];
%                         end
%                     end
%                     if eq(VizOpts.Registered, bitand(VizOpts.Registered, viz))
%                         [x,y]=meshgrid(linspace(v0(1,1),v0(3,1),50),linspace(v0(1,2),v0(2,2),50));
%                         z = ones(size(x))*v0(1,3);
%                         [xt,yt,zt]=transformPointsForward(a3d,x(:),y(:),z(:));
%                         xt = reshape(xt,size(x));
%                         yt = reshape(yt,size(y));
%                         zt = reshape(zt,size(z));
% 
%                         ph=slice(Xf,Yf,Zf,bscale*livenuclei_post_global_t,xt,yt,zt);
%                         set(ph, 'FaceColor', 'interp', 'EdgeColor', 'none','FaceAlpha',0.3, 'EdgeAlpha',0.3);
%                     else
%                         patch('Faces',[1 2 4 3], 'Vertices', v, 'FaceColor','cyan','FaceAlpha',0.3, 'EdgeColor','none');
%                     end
%                     line(v([1 3],1), v([1 3],2), v([1 3], 3),'Color','red','LineWidth',2);
%                     line(v([1 2],1), v([1 2],2), v([1 2], 3), 'Color','blue','LineWidth',2);
%                 end
%             end
%         end
%         close(wbf);
%         n=length(vertcat(cell2measurement{empty_cells}));
%         headers = {'Filename','x','y','z','patchx','patchy','Cell'};
%         T = cell2table(cell(0,7));
%         T.Properties.VariableNames = headers;
%         for j = 1:length(empty_cells)
%             locs=cell2measurement{empty_cells(j)};
%             for i = 1:length(locs)
%                 T=[T; {replace(live_img_series_infos.Properties.RowNames{locs(i)},'Series','S'),[],[],[],[],[],empty_cells(j)}];
%             end
%         end
%         writetable(T, missing_cell_mechano_measurement_table_output_fpath);
%     end
    axis equal;
    axis vis3d;
    hold off;
    set(gcf,'Position', [1639 1107 946 933])
end

% uncomment below to save cellxgene table - Payman
% writetable(cellxgene, cellcounts_table_output_fpath);
%cellxgenexmech_dir = sprintf('%s/joint_analysis/cellxgenexmech',wd);
%if ~exist(cellxgenexmech_dir, 'dir') mkdir(cellxgenexmech_dir), end
%writetable(cellxgene, sprintf('%s/cell_counts_table_plate_%d_well%d_embryo_%d_with_msd_round%d_dropUncertain_withSD.txt', cellxgenexmech_dir, plate, well, embryo, round));


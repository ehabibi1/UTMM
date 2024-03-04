% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Embryo Ex-Seq (non-ExM) registration of FoV
function[] = chech_seg(home_dir, expt, species,plate,well,embryo, bounderies, ReReg)

series = 1;
timepoint = 1;
p_w_e = sprintf('plate_%d_well_%d_embryo_%d', plate, well, embryo);
%num_embryos = size(embryos,2);
embryo_dir = sprintf('%s/projects/%s/%s/transcriptome/plate_%d_well_%d/embryo_%d',home_dir,species,expt,plate,well,embryo);

%embryo_bounds = bounderies(p_w_e).hyb;
embryo_bounds = dlmread(sprintf('%s/processed/bounds.txt',embryo_dir));
%% register Dapi from morphology to dapi from hyb

%% Load morphology staining (membrane and dapi)
reg_seg_dir = sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/',home_dir,species,expt, plate, well, embryo);    

if ReReg==1
    ReReg_seg_dir = sprintf('%s/ReReg',reg_seg_dir);
    hyb_dapi = uint16(read_3d_tif(sprintf('%s/processed/offset/cy01_dapi.tif',embryo_dir),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6)));
    stain_dapi_stack = read_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_dapi.tif', ReReg_seg_dir, plate, well, embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
    membSeg_stack = read_3d_tif(sprintf('%s/p%dw%de%d-membSeg.tif',ReReg_seg_dir, plate, well, embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
else
    hyb_dapi = uint16(read_3d_tif(sprintf('%s/processed/offset/cy01_dapi.tif',embryo_dir),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6)));
    stain_dapi_stack = read_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_dapi.tif', reg_seg_dir, plate, well, embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
    membSeg_stack = read_3d_tif(sprintf('%s/p%dw%de%d-membSeg.tif',reg_seg_dir, plate, well, embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
end

fig = figure;
imshowpair(capImage(max(hyb_dapi,[],3),99,'prc'),capImage(max(stain_dapi_stack,[],3),99,'prc'))
title(sprintf('embryo %d corrMax', embryo))

if ReReg==1
    saveas(fig,sprintf('%s/Dapi_corrMax_p%dw%de%d.png',ReReg_seg_dir, plate, well, embryo));
else
    saveas(fig,sprintf('%s/Dapi_corrMax_p%dw%de%d.png',reg_seg_dir, plate, well, embryo));
end

end

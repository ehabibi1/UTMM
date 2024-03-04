% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Embryo Ex-Seq (non-ExM) registration of FoV
function[] = reRegister_process_reg_hyb_to_morpho_server_finereg(home_dir, expt,species, plate,well,embryo, bounderies)

%% set up environment
tic
%home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\ExSeq_Zack';
%data_dir = '\\sodium\broad_thechenlab\ehsan\DataRepository\SplintR';
%SERVER
%home_dir = '/broad/thechenlab/ehsan/analysis/ExSeq_Zack';
%data_dir = '/broad/thechenlab/ehsan/DataRepository/SplintR';
cd(home_dir);
addpath(genpath('scripts/'))
 
%%
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Registering morphology to hyb --> plate %d well %d embryo %d ...',time, plate, well, embryo));

%% Load hyb DAPI

% all files have 1 series, 1 timepoint for now

series = 1;
timepoint = 1;
p_w_e = sprintf('plate_%d_well_%d_embryo_%d', plate, well, embryo);
%num_embryos = size(embryos,2);
embryo_dir = sprintf('%s/projects/%s/%s/transcriptome/plate_%d_well_%d/embryo_%d',home_dir,species,expt,plate,well,embryo);

%embryo_bounds = bounderies(p_w_e).hyb;
embryo_bounds = dlmread(sprintf('%s/processed/bounds.txt',embryo_dir));
%% register Dapi from morphology to dapi from hyb

%% Load morphology staining (membrane and dapi)

hyb_dapi = uint16(read_3d_tif(sprintf('%s/processed/offset/cy01_dapi.tif',embryo_dir),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6)));

reg_seg_dir = sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/',home_dir,species,expt, plate, well, embryo);    
stain_Ecad_stack = read_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Ecad.tif', reg_seg_dir, plate, well, embryo), embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
stain_Ezr_stack = read_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Ezr.tif', reg_seg_dir, plate, well, embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
stain_dapi_stack = read_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_dapi.tif', reg_seg_dir, plate, well, embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
stain_merged_EzrEcad = read_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_EcadEzr_merged.tif', reg_seg_dir, plate, well, embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
membSeg_stack = read_3d_tif(sprintf('%s/p%dw%de%d-membSeg.tif', reg_seg_dir, plate, well, embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));


ReReg_seg_dir = sprintf('%s/ReReg',reg_seg_dir);
if ~exist(ReReg_seg_dir, 'dir') mkdir(ReReg_seg_dir), end


fig = figure;
imshowpair(capImage(max(hyb_dapi,[],3),99,'prc'),capImage(max(stain_dapi_stack,[],3),99,'prc'))
title(sprintf('embryo %d corrMax', embryo))
saveas(fig,sprintf('%s/Dapi_corrMax_p%dw%de%d.png',ReReg_seg_dir, plate, well, embryo));




%% CROP 

%[stain_dapi_offset_stack2, tform] = corrRegMax_morph(hyb_dapi, stain_dapi_stack);



%offset_stain_Ecad_stack2 = imwarp3D(hyb_dapi, stain_Ecad_stack, tform);
%offset_stain_Ezr_stack2 = imwarp3D(hyb_dapi, stain_Ezr_stack, tform);
%offset_stain_merged_EzrEcad = imwarp3D(hyb_dapi, stain_merged_EzrEcad, tform);
%offset_membSeg_stack = imwarp3D_labelMatrix(hyb_dapi, membSeg_stack, tform);

[optimizer, metric] = imregconfig('monomodal');

optimizer.MaximumStepLength = 0.01;

tform = imregtform(stain_dapi_stack, hyb_dapi, 'affine', optimizer, metric);
stain_dapi_offset_stack2 = imwarp(stain_dapi_stack,tform,'OutputView',imref3d(size(hyb_dapi)));
 
fig = figure;
imshowpair(capImage(max(hyb_dapi,[],3),99,'prc'),capImage(max(stain_dapi_offset_stack2,[],3),99,'prc'))
title(sprintf('embryo %d affine', embryo))
saveas(fig,sprintf('%s/Dapi_affine_p%dw%de%d.png',ReReg_seg_dir, plate, well, embryo));



 offset_stain_Ecad_stack2 = imwarp(stain_Ecad_stack,tform,'OutputView',imref3d(size(hyb_dapi)));
 offset_stain_Ezr_stack2 = imwarp(stain_Ezr_stack,tform,'OutputView',imref3d(size(hyb_dapi)));
 offset_stain_merged_EzrEcad = imwarp(stain_merged_EzrEcad,tform,'OutputView',imref3d(size(hyb_dapi)));
 offset_membSeg_stack = imwarp(membSeg_stack,tform,'nearest','OutputView',imref3d(size(hyb_dapi))); 

time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Writing registered files ...',time));


write_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_dapi.tif',ReReg_seg_dir, plate,well,embryo),stain_dapi_offset_stack2);
write_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Ecad.tif',ReReg_seg_dir, plate,well,embryo),offset_stain_Ecad_stack2);
write_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Ezr.tif',ReReg_seg_dir, plate,well,embryo),offset_stain_Ezr_stack2);
write_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_merged_EzrEcad.tif',ReReg_seg_dir, plate,well,embryo),offset_stain_merged_EzrEcad);
write_3d_tif(sprintf('%s/p%dw%de%d-membSeg.tif',ReReg_seg_dir, plate,well,embryo),offset_membSeg_stack);


filename = sprintf('%s/regtform_Rereg_hybdapi_to_membSeg.txt',ReReg_seg_dir);

dlmwrite(filename,tform.T)

disp(sprintf('%s: saved stain stacks',sec2time(toc)))

end

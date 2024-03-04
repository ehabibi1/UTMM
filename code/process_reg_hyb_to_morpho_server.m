% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Embryo Ex-Seq (non-ExM) registration of FoV
function[] = process_reg_hyb_to_morpho_server(data_dir, home_dir, expt,plate,well,embryos_in_fov,embryo, bounderies, Ecad_channel, Ezr_channel, dapi_channel)


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
embryo_dir = sprintf('%s/plate_%d_well_%d/embryo_%d',expt,plate,well,embryo);

%embryo_bounds = bounderies(p_w_e).hyb;
embryo_bounds = dlmread(sprintf('%s/processed/bounds.txt',embryo_dir));
%% register Dapi from morphology to dapi from hyb

%% Load morphology staining (membrane and dapi)

stain_reader = bfGetReader(sprintf('%s/%s/ImageData/Plate %d, Well %d/P%d_W%d_Embryo_%s/1A -- Morphology/*_Well_%d_Embryo_%s_*.tif',data_dir,expt,plate,well,plate,well,embryos_in_fov,well,embryos_in_fov));
stain_xlen = stain_reader.getSizeX; stain_ylen = stain_reader.getSizeY; stain_zlen = stain_reader.getSizeZ; 

stain_Ecad_stack = zeros(stain_ylen,stain_xlen,stain_zlen,'uint16');
%stain_Zo1_stack = zeros(stain_ylen,stain_xlen,stain_zlen,'uint16');
stain_Ezr_stack = zeros(stain_ylen,stain_xlen,stain_zlen,'uint16');
stain_dapi_stack = zeros(stain_ylen,stain_xlen,stain_zlen,'uint16');

for z=1:stain_zlen
    stain_Ecad_stack(:,:,z) = readPlane(stain_reader,series,z,Ecad_channel,timepoint);
    %stain_Zo1_stack(:,:,z) = readPlane(stain_reader,series,z,2,timepoint);
    stain_Ezr_stack(:,:,z) = readPlane(stain_reader,series,z,Ezr_channel,timepoint);
    stain_dapi_stack(:,:,z) = readPlane(stain_reader,series,z,dapi_channel,timepoint);
end

disp(sprintf('%s: loaded stain stack + DAPI',sec2time(toc)))

%% CROP 

%if length(split(embryos_in_fov,"+"))>1
    time = datestr(datetime(now,'ConvertFrom','datenum'));
    disp(sprintf('%s ----- Loading files ...',time));

    %embryo_bounds = dlmread(sprintf('%s/plate_%d_well_%d/embryo_%d/processed/bounds.txt',expt,plate,well,embryo));
    hyb_dapi = uint16(read_3d_tif(sprintf('%s/processed/offset/hyb_dapi.tif',embryo_dir),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6)));

    time = datestr(datetime(now,'ConvertFrom','datenum'));
    disp(sprintf('%s ----- Cropping morphology files ...',time));
    stain_Ecad_stack = imcrop_xyz(stain_Ecad_stack,bounderies(p_w_e).moropho);
    %stain_Zo1_stack = imcrop_xyz(stain_Zo1_stack,bounderies(p_w_e).moropho);
    stain_Ezr_stack = imcrop_xyz(stain_Ezr_stack,bounderies(p_w_e).moropho);
    stain_dapi_stack = imcrop_xyz(stain_dapi_stack,bounderies(p_w_e).moropho);
% else
%     time = datestr(datetime(now,'ConvertFrom','datenum'));
%     disp(sprintf('%s ----- Reading Hyb ... ', time));
%     hyb_reader = bfGetReader(sprintf('%s/%s/ImageData/Plate %d, Well %d/P%d_W%d_Embryo_%s/9 -- RolonyDetection/*_Well_%d_Embryo_%s_*.tif',data_dir,expt,plate,well,plate,well,embryos_in_fov,well,embryos_in_fov));
%     num_hyb_channels = hyb_reader.getSizeC;
%     hyb_dapi_channel = 2;
%     xlen(1) = hyb_reader.getSizeX; ylen(1) = hyb_reader.getSizeY; zlen(1) = hyb_reader.getSizeZ;
% 
%     hyb_dapi = zeros(ylen(1),xlen(1),zlen(1),'uint16');
%     
%     for z=1:zlen(1)
%         hyb_dapi(:,:,z) = readPlane(hyb_reader,series,z,hyb_dapi_channel,timepoint);
%     end
% end 
%% Offset DAPI and apply to stain
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Calculating and applying offsets ...',time));
%[stain_dapi_offset_stack stain_offsets] = get_offset_xyz(dapi_offset_stacks{1},stain_dapi_stack,max_dims);
[stain_dapi_offset_stack stain_offsets] = get_offset_xyz(hyb_dapi,stain_dapi_stack,[embryo_bounds(4), embryo_bounds(5), embryo_bounds(6)]);
offset_stain_Ecad_stack = apply_offset_xyz(stain_Ecad_stack,stain_offsets,[embryo_bounds(4), embryo_bounds(5), embryo_bounds(6)]);
offset_stain_Ezr_stack = apply_offset_xyz(stain_Ezr_stack,stain_offsets,[embryo_bounds(4), embryo_bounds(5), embryo_bounds(6)]);

%offset_stain_Zo1_stack = apply_offset_xyz(stain_Zo1_stack,stain_offsets,[embryo_bounds(4), embryo_bounds(5), embryo_bounds(6)]);
%figure; imshowpair(capImage(max(hyb_dapi,[],3),99,'prc'),capImage(max(offset_stain_Zo1_stack,[],3),99,'prc'))

%############################

time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- affine registration ...',time));
[optimizer, metric] = imregconfig('monomodal');

tform = imregtform(stain_dapi_offset_stack, hyb_dapi, 'affine', optimizer, metric);
stain_dapi_offset_stack2 = imwarp(stain_dapi_offset_stack,tform,'OutputView',imref3d(size(hyb_dapi)));
offset_stain_Ecad_stack2 = imwarp(offset_stain_Ecad_stack,tform,'OutputView',imref3d(size(hyb_dapi)));
offset_stain_Ezr_stack2 = imwarp(offset_stain_Ezr_stack,tform,'OutputView',imref3d(size(hyb_dapi)));



% figure; imshowpair(capImage(max(hyb_dapi,[],3),99,'prc'),capImage(max(stain_dapi_offset_stack,[],3),99,'prc'))
% figure; imshowpair(capImage(max(hyb_dapi,[],3),99,'prc'),capImage(max(stain_dapi_offset_stack2,[],3),99,'prc'))
% 
% figure; imshowpair(capImage(max(hyb_dapi,[],3),99,'prc'),capImage(max(offset_stain_Ecad_stack,[],3),99,'prc'))
% figure; imshowpair(capImage(max(hyb_dapi,[],3),99,'prc'),capImage(max(offset_stain_Ecad_stack2,[],3),99,'prc'))
% 
% figure; imshowpair(capImage(max(hyb_dapi,[],3),99,'prc'),capImage(max(offset_stain_Ezr_stack,[],3),99,'prc'))
% figure; imshowpair(capImage(max(hyb_dapi,[],3),99,'prc'),capImage(max(offset_stain_Ezr_stack2,[],3),99,'prc'))
% 
% 
% figure; imshowpair(capImage(max(offset_stain_Ecad_stack,[],3),99,'prc'),capImage(max(offset_stain_Ezr_stack,[],3),99,'prc'))
% figure; imshowpair(capImage(max(offset_stain_Ecad_stack2,[],3),99,'prc'),capImage(max(offset_stain_Ezr_stack2,[],3),99,'prc'))


disp(sprintf('%s: applied offset to stain stack',sec2time(toc)))

reg_seg_dir = sprintf('%s/%s/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d',home_dir,expt, plate, well, embryo);
if ~exist(reg_seg_dir, 'dir') mkdir(reg_seg_dir), end

time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Writing registered files ...',time));
write_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_dapi.tif',reg_seg_dir, plate,well,embryo),stain_dapi_offset_stack2);
write_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Ecad.tif',reg_seg_dir, plate,well,embryo),offset_stain_Ecad_stack2);
write_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Ezr.tif',reg_seg_dir, plate,well,embryo),offset_stain_Ezr_stack2);
%write_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_dapi.tif',reg_seg_dir, plate,well,embryo),imcrop_xyz(stain_dapi_offset_stack2,embryo_bounds));
%write_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Ecad.tif',reg_seg_dir, plate,well,embryo),imcrop_xyz(offset_stain_Ecad_stack2,embryo_bounds));
%write_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Zo1.tif',reg_seg_dir, plate,well,embryo),imcrop_xyz(offset_stain_Zo1_stack,embryo_bounds));
%write_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Ezr.tif',reg_seg_dir, plate,well,embryo),imcrop_xyz(offset_stain_Ezr_stack2,embryo_bounds));

disp(sprintf('%s: saved stain stacks',sec2time(toc)))

end

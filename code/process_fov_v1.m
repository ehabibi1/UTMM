% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Embryo Ex-Seq (non-ExM) registration of FoV

function[] = process_fov_v1(data_dir, home_dir, raw_data_dir_name, expt,species, plate,well,embryos_in_fov,embryo, bounderies)

%% set up environment
tic
%home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\InSitu_preImpEmbryo';
%data_dir = '\\sodium\broad_thechenlab\ehsan\DataRepository\SplintR';
%SERVER
%home_dir = '/broad/thechenlab/ehsan/analysis/InSitu_preImpEmbryo/';
%data_dir = '/broad/thechenlab/ehsan/DataRepository/SplintR/';
cd(home_dir);
addpath(genpath('scripts/'))
 
%%
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Processing plate %d well %d embryo %d ...',time, plate, well, embryo));

%% Load hyb DAPI

% all files have 1 series, 1 timepoint for now

series = 1;
timepoint = 1;
p_w_e = sprintf('plate_%d_well_%d_embryo_%d', plate, well, embryo);
%num_embryos = size(embryos,2);
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Reading Hyb ... ', time));
hyb_reader = bfGetReader(sprintf('%s/%s/ImageData/Plate %d, Well %d/P%d_W%d_Embryo_%s/9 -- RolonyDetection/*_Well_%d_Embryo_%s_*.tif',data_dir,raw_data_dir_name,plate,well,plate,well,embryos_in_fov,well,embryos_in_fov));
current_file = hyb_reader.getCurrentFile();
disp(sprintf('%s ----- Loading ... %s',time, current_file));
num_hyb_channels = hyb_reader.getSizeC;
hyb_dapi_channel = 2;
xlen(1) = hyb_reader.getSizeX; ylen(1) = hyb_reader.getSizeY; zlen(1) = hyb_reader.getSizeZ;

dapi_stacks{1} = zeros(ylen(1),xlen(1),zlen(1),'uint16');
    
for z=1:zlen(1)
    dapi_stacks{1}(:,:,z) = readPlane(hyb_reader,series,z,hyb_dapi_channel,timepoint);
end

%% set region of interest (per embryo)

figure; imshow(capImage(max(dapi_stacks{1},[],3),99,'prc'),[])


%crop test
embryo_bounds = bounderies(p_w_e).hyb;
crop_test = imcrop_xyz(dapi_stacks{1},embryo_bounds(1,:));
figure; imshow(capImage(max(crop_test,[],3),99,'prc'),[])

%% make embryo directories

%for i=1:num_embryos

embryo_dir = sprintf('%s/projects/%s/%s/transcriptome/plate_%d_well_%d/embryo_%d',home_dir,species,expt,plate,well,embryo);

mkdir(embryo_dir)

time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Making directories ... ', time));
fig_dir = sprintf('%s/figure',embryo_dir);
mkdir(fig_dir)
processed_dir = sprintf('%s/processed',embryo_dir);
mkdir(processed_dir)    

offset_dir = sprintf('%s/processed/offset',embryo_dir);
mkdir(offset_dir)
    
deconv_dir = sprintf('%s/deconv',embryo_dir);
mkdir(deconv_dir)
    
%end


%for i=1:num_embryos
%dlmwrite(sprintf('%s/processed//bounds.txt',embryo_dir),embryo_bounds);
%end



%% save hyb images

hyb_probe_channel = 1;

hyb_probe_stack = zeros(ylen(1),xlen(1),zlen(1),'uint16');
    
for z=1:zlen(1)
    hyb_probe_stack(:,:,z) = readPlane(hyb_reader,series,z,hyb_probe_channel,timepoint);
end

% load hyb probe and stain

%for i=1:num_embryos
    
write_3d_tif(sprintf('%s/processed/offset/hyb_probe.tif',embryo_dir),imcrop_xyz(hyb_probe_stack,embryo_bounds));
write_3d_tif(sprintf('%s/processed/offset/hyb_dapi.tif',embryo_dir),imcrop_xyz(dapi_stacks{1},embryo_bounds));
    
%end

% code to test your crop
%this if condition is needed when you crop all images even if there is one
%embyo in the FOV
%if length(split(embryos_in_fov,"+"))>1
       dapi_stacks{1} = imcrop_xyz(dapi_stacks{1},embryo_bounds);
       hyb_probe_stack = imcrop_xyz(hyb_probe_stack,embryo_bounds);
       embryo_bounds(1,1) = 1;
       embryo_bounds(1,2) = 1;
       xlen(1) = embryo_bounds(1,4);
       ylen(1) = embryo_bounds(1,5);
       zlen(1) = embryo_bounds(1,6);
       dlmwrite(sprintf('%s/processed/bounds.txt',embryo_dir),embryo_bounds);
%end

%dlmwrite(sprintf('%s/processed/bounds.txt',embryo_dir),embryo_bounds);

%for i=1:num_embryos
%dlmwrite(sprintf('%s/processed/bounds.txt',embryo_dir),embryo_bounds);
%end

figure; imshow(capImage(max(dapi_stacks{1},[],3),99,'prc'),[])
figure; imshow(capImage(max(hyb_probe_stack,[],3),99,'prc'),[])

%figure; 
%for z=1:size(crop_test,3)
%    imshow(capImage(crop_test(:,:,z),99,'prc'),[])
%end


%clearvars hyb_probe_stack hyb_stain_stack

disp(sprintf('%s: wrote hyb and assoc. images',sec2time(toc)))


%% load data
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Loading seq data ... ', time));
% set channel and cycle info

dapi_channel = 5;
num_cycles = 7; %Round
%offsets = {};
%dapi_offset_stacks = {};
dapi_offset_stacks{1} = dapi_stacks{1};

figure; imshow(capImage(max(dapi_offset_stacks{1},[],3),99,'prc'),[])

for cycle=1:num_cycles
    
    % load reader
    %cycle_number_dir = cycle+1;
    %seq_reader = bfGetReader(sprintf('%s/MECHANICS (May 2021)/ImageData/Plate %d, Well %d/P%d_W%d_Embryo_%s/2 -- SeqN-1_Ligation_1/*_Well_%d_Embryo_%d_*.tif',data_dir,plate,well,plate,well,embryos_in_fov,well,embryos(1)));
    
    if cycle==1
        seq_reader = bfGetReader(sprintf('%s/%s/ImageData/Plate %d, Well %d/P%d_W%d_Embryo_%s/2 -- SeqN-1_Ligation_1/*_Well_%d_Embryo_%s_*.tif',data_dir,raw_data_dir_name,plate,well,plate,well,embryos_in_fov,well,embryos_in_fov));
    elseif cycle==2
        seq_reader = bfGetReader(sprintf('%s/%s/ImageData/Plate %d, Well %d/P%d_W%d_Embryo_%s/3 -- SeqN-1_Ligation_2/*_Well_%d_Embryo_%s_*.tif',data_dir,raw_data_dir_name,plate,well,plate,well,embryos_in_fov,well,embryos_in_fov));
    elseif cycle==3
        seq_reader = bfGetReader(sprintf('%s/%s/ImageData/Plate %d, Well %d/P%d_W%d_Embryo_%s/4 -- SeqN_Ligation_1/*_Well_%d_Embryo_%s_*.tif',data_dir,raw_data_dir_name,plate,well,plate,well,embryos_in_fov,well,embryos_in_fov));
    elseif cycle==4
        seq_reader = bfGetReader(sprintf('%s/%s/ImageData/Plate %d, Well %d/P%d_W%d_Embryo_%s/5 -- SeqN_Ligation_2/*_Well_%d_Embryo_%s_*.tif',data_dir,raw_data_dir_name,plate,well,plate,well,embryos_in_fov,well,embryos_in_fov));
    elseif cycle==5
        seq_reader = bfGetReader(sprintf('%s/%s/ImageData/Plate %d, Well %d/P%d_W%d_Embryo_%s/6 -- SeqN-2_Ligation_2/*_Well_%d_Embryo_%s_*.tif',data_dir,raw_data_dir_name,plate,well,plate,well,embryos_in_fov,well,embryos_in_fov));
    elseif cycle==6
        seq_reader = bfGetReader(sprintf('%s/%s/ImageData/Plate %d, Well %d/P%d_W%d_Embryo_%s/7 -- SeqN-3_Ligation_2/*_Well_%d_Embryo_%s_*.tif',data_dir,raw_data_dir_name,plate,well,plate,well,embryos_in_fov,well,embryos_in_fov));
    elseif cycle==7
        seq_reader = bfGetReader(sprintf('%s/%s/ImageData/Plate %d, Well %d/P%d_W%d_Embryo_%s/8 -- SeqN-4_Ligation_2/*_Well_%d_Embryo_%s_*.tif',data_dir,raw_data_dir_name,plate,well,plate,well,embryos_in_fov,well,embryos_in_fov));
    end
    
    current_file = seq_reader.getCurrentFile;
    disp(sprintf('%s ----- Loading seq data %s ... ', time, current_file));
    
    xlen(cycle+1) = seq_reader.getSizeX; ylen(cycle+1) = seq_reader.getSizeY; zlen(cycle+1) = seq_reader.getSizeZ; 
  
    num_channels = seq_reader.getSizeC;
    
    % load DAPI
    
    dapi_stacks{cycle+1} = zeros(ylen(cycle+1),xlen(cycle+1),zlen(cycle+1),'uint16');
    
    for z=1:zlen(cycle+1)
        dapi_stacks{cycle+1}(:,:,z) = readPlane(seq_reader,series,z,dapi_channel,timepoint);
    end
    
    %if length(split(embryos_in_fov,"+"))>1
        c = sprintf('cycle_%d',cycle);
        dapi_stacks{cycle+1} = imcrop_xyz(dapi_stacks{cycle+1},bounderies(p_w_e).(c));
        xlen(cycle+1) = bounderies(p_w_e).(c)(4); 
        ylen(cycle+1) = bounderies(p_w_e).(c)(5); 
        zlen(cycle+1) = bounderies(p_w_e).(c)(6);
        
    %end
    
    max_dims = [max(xlen) max(ylen) max(zlen)];
    
    cycle
    
    figure; imshow(capImage(max(dapi_stacks{cycle+1},[],3),99,'prc'),[])
    

    figure; imshow(capImage(max(dapi_stacks{cycle+1},[],3),99,'prc'),[])
    figure; imshow(capImage(max(dapi_offset_stacks{1},[],3),99,'prc'),[])
    figure; imshowpair(capImage(max(dapi_offset_stacks{1},[],3),99,'prc'),capImage(max(dapi_stacks{cycle+1},[],3),99,'prc'))
    
    disp(sprintf('%s: loaded DAPI for cycle %d',sec2time(toc),cycle))
    
    % offset DAPI

    [dapi_offset_stacks{cycle+1} offsets{cycle+1}] = get_offset_xyz(dapi_offset_stacks{1},dapi_stacks{cycle+1},max_dims);
    figure; imshow(capImage(max(dapi_offset_stacks{cycle+1},[],3),99,'prc'),[])
    figure; imshowpair(capImage(max(dapi_offset_stacks{1},[],3),99,'prc'),capImage(max(dapi_offset_stacks{cycle+1},[],3),99,'prc'))
    
    disp(sprintf('%s: offset DAPI for cycle %d',sec2time(toc),cycle))
    
    % write DAPI
    %for i=1:num_embryos
        offset_dir = sprintf('%s/processed/offset/',embryo_dir);
        write_3d_tif(sprintf('%s/cy%02d_dapi.tif',offset_dir,cycle),imcrop_xyz(dapi_offset_stacks{cycle+1},embryo_bounds));
    %end
    
    disp(sprintf('%s: wrote DAPI for cycle %d',sec2time(toc),cycle))
    
    % loop through seq channels
    disp(sprintf('%s: loading seq data %d',sec2time(toc),cycle))
    channels = 1:num_channels; channels(channels==dapi_channel) = [];
    for channel=[channels]
        
        % load seq
        if length(split(embryos_in_fov,"+"))>1
            seq_stack = zeros(seq_reader.getSizeY,seq_reader.getSizeX,seq_reader.getSizeZ,'uint16');
            for z=1:seq_reader.getSizeZ
                seq_stack(:,:,z) = readPlane(seq_reader,series,z,channel,timepoint);
            end
            c = sprintf('cycle_%d',cycle);
            seq_stack = imcrop_xyz(seq_stack,bounderies(p_w_e).(c));
        else
            seq_stack = zeros(ylen(cycle+1),xlen(cycle+1),zlen(cycle+1),'uint16');
            for z=1:zlen(cycle+1)
                seq_stack(:,:,z) = readPlane(seq_reader,series,z,channel,timepoint);
            end
            
        end 
        
        disp(sprintf('%s: loaded cycle %d, channel %d',sec2time(toc),cycle, channel))
        
        % offset seq
        
        offset_seq_stack = apply_offset_xyz(seq_stack,offsets{cycle+1},max_dims);
        figure; imshowpair(capImage(max(dapi_offset_stacks{2},[],3),99,'prc'),capImage(max(offset_seq_stack,[],3),99,'prc'))
        
        disp(sprintf('%s: offset cycle %d, channel %d',sec2time(toc),cycle, channel))
        
        % write seq
        
        %for i=1:num_embryos
            offset_dir = sprintf('%s/processed/offset/',embryo_dir);
            write_3d_tif(sprintf('%s/cy%02d_ch%02d.tif',offset_dir,cycle,channel),imcrop_xyz(offset_seq_stack,embryo_bounds));
        %end
        
        disp(sprintf('%s: wrote cycle %d, channel %d',sec2time(toc),cycle, channel))
        
    end
    
end
% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Embryo Ex-Seq (non-ExM) registration of FoV

function[] = process_fov()
%% define parameters

expt = 'Mechno_Transcriptome_trial1_Jan2021';

%Well 1, Embryo 3, and Well 2 Embryos 1 & 2 look good!
well = 2;
embryos = [1]; 
embryo = 1;
%% set up environment

clearvars -except expt well embryos first_embryo
tic
home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\ExSeq_Zack';
data_dir = '\\sodium\broad_thechenlab\ehsan\DataRepository\SplintR';
cd(home_dir);
addpath(genpath('scripts/'))
%% Load hyb DAPI

% all files have 1 series, 1 timepoint for now

series = 1;
timepoint = 1;

num_embryos = size(embryos,2);
hyb_reader = bfGetReader(sprintf('%s/%s/Mechanics, then ExSeq (Jan 2021)/Data/Well %d, Embryo %d/Rolony Detection/*Embryo_%d*.tif',data_dir,expt,well,embryo,embryo));
%hyb_reader = bfGetReader(sprintf('\\\\sodium/broad_thechenlab/ehsan/Ehsan_Anu_ExSeq_trial_2/without_ab_staining/R0_detection_hyb/*well%d*embryo_%d*.tif',well,first_embryo));
num_hyb_channels = hyb_reader.getSizeC;
hyb_dapi_channel = 2;
%hyb_dapi_channel = 3;
xlen(1) = hyb_reader.getSizeX; ylen(1) = hyb_reader.getSizeY; zlen(1) = hyb_reader.getSizeZ;

dapi_stacks{1} = zeros(ylen(1),xlen(1),zlen(1),'uint16');
    
for z=1:zlen(1)
    dapi_stacks{1}(:,:,z) = readPlane(hyb_reader,series,z,hyb_dapi_channel,timepoint);
end

%% set region of interest (per embryo)

figure; imshow(capImage(max(dapi_stacks{1},[],3),99,'prc'),[])

%% set embryo bounds, manual for now
%[xmin ymin zmin width height depth]
%embryo_bounds = [...
%    500 1000 1 700 700 218; ... % embryo 1
%    1025 475 1 700 700 218 ...  % embryo 2
%    ];

%%wel1_2 embryo 1
xmin = 645;
xmax = 1488;
ymin = 669;
ymax = 1494;
zmin = 1;
zmax = 218;

%%well_1  embryo 3
xmin = 800;
xmax = 1476;
ymin = 790;
ymax = 1474;
zmin = 1;
zmax = 348;

width = xmax-xmin;
height = ymax-ymin;
depth = zmax;
embryo_bounds = [xmin ymin zmin width height depth]; 

%embryo_bounds = [...
%    400 750 1 700 700 313; ... % well 1, embryo 1
%    1100 850 1 700 700 313 ...  % well 1, embryo 2
%    ];

%embryo_bounds = [...
%    650 1 1 700 700 258; ... % well 1, embryo 3
%    150 400 1 700 700 258; ...  % well 1, embryo 4
%    1025 1345 1 700 700 258 ...  % well 1, embryo 5
%    ];

%embryo_bounds = [725 675 1 700 700 244]; % well 3, embryo 1
%embryo_bounds = [650 725 1 700 700 230]; % well 3, embryo 2
%embryo_bounds = [475 900 1 700 700 191]; % well 3, embryo 2

% code to test your crop

crop_test = imcrop_xyz(dapi_stacks{1},embryo_bounds(1,:));
figure; imshow(capImage(max(crop_test,[],3),99,'prc'),[])

%figure; 
%for z=1:size(crop_test,3)
%    imshow(capImage(crop_test(:,:,z),99,'prc'),[])
%end



%% make embryo directories

for i=1:num_embryos

    embryo_dir = sprintf('%s/well_%d/embryo_%d',expt,well,embryos(i));
    %embryo_dir = sprintf('processed/%s/well%d_embryo%d',expt,well,embryos(i));
    mkdir(embryo_dir)


    fig_dir = sprintf('%s/figure',embryo_dir);
    mkdir(fig_dir)
    processed_dir = sprintf('%s/processed',embryo_dir);
    mkdir(processed_dir)    
    
    offset_dir = sprintf('%s/processed/offset',embryo_dir);
    mkdir(offset_dir)

    reg_dir = sprintf('%s/reg',embryo_dir);
    if ~exist(reg_dir, 'dir') mkdir(reg_dir), end


    
end


for i=1:num_embryos
    dlmwrite(sprintf('%s/well_%d/embryo_%d/processed/bounds.txt',expt,well,embryos(i)),embryo_bounds(i,:));
end



%% save hyb images

hyb_probe_channel = 1;
%hyb_stain_channel = 2;

hyb_probe_stack = zeros(ylen(1),xlen(1),zlen(1),'uint16');
%hyb_stain_stack = zeros(ylen(1),xlen(1),zlen(1),'uint16');
    
for z=1:zlen(1)
    hyb_probe_stack(:,:,z) = readPlane(hyb_reader,series,z,hyb_probe_channel,timepoint);
    %hyb_stain_stack(:,:,z) = readPlane(hyb_reader,series,z,hyb_stain_channel,timepoint);
end

% load hyb probe and stain

for i=1:num_embryos
    
    write_3d_tif(sprintf('%s/well_%d/embryo_%d/processed/offset/hyb_probe.tif',expt,well,embryos(i)),imcrop_xyz(hyb_probe_stack,embryo_bounds(i,:)));
    %write_3d_tif(sprintf('%s/hyb_stain.tif',offset_dir),imcrop_xyz(hyb_stain_stack,embryo_bounds(i,:)));
    write_3d_tif(sprintf('%s/well_%d/embryo_%d/processed/offset/hyb_dapi.tif',expt,well,embryos(i)),imcrop_xyz(dapi_stacks{1},embryo_bounds(i,:)));
    
end

clearvars hyb_probe_stack hyb_stain_stack

disp(sprintf('%s: wrote hyb and assoc. images',sec2time(toc)))

%% load data

% set channel and cycle info

dapi_channel = 5;
num_cycles = 4; %Round
%offsets = {};
%dapi_offset_stacks = {};
dapi_offset_stacks{1} = dapi_stacks{1};

figure; imshow(capImage(max(dapi_offset_stacks{1},[],3),99,'prc'),[])

for cycle=1:num_cycles
    
    % load reader
    seq_reader = bfGetReader(sprintf('%s/%s/Mechanics, then ExSeq (Jan 2021)/Data/Well %d, Embryo %d/R%d/*Well_%d_*Embryo_%d*.tif',data_dir,expt,well,embryo,cycle,well,embryo));
    %seq_reader = bfGetReader(sprintf('\\\\sodium/broad_thechenlab/ehsan/Ehsan_Anu_ExSeq_trial_2/without_ab_staining/R%d/*well%d*embryo_%d*.tif',cycle,well,first_embryo));
    %seq_reader = bfGetReader(sprintf('data/%s/R%d/*well%d*embryo_%d*.tif',expt,cycle,well,first_embryo));
    xlen(cycle+1) = seq_reader.getSizeX; ylen(cycle+1) = seq_reader.getSizeY; zlen(cycle+1) = seq_reader.getSizeZ; 

    max_dims = [max(xlen) max(ylen) max(zlen)];
    num_channels = seq_reader.getSizeC;
    
    % load DAPI
    
    dapi_stacks{cycle+1} = zeros(xlen(cycle+1),ylen(cycle+1),zlen(cycle+1),'uint16');
    for z=1:zlen(cycle+1)
        dapi_stacks{cycle+1}(:,:,z) = readPlane(seq_reader,series,z,dapi_channel,timepoint);
    end
    
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
    for i=1:num_embryos
        offset_dir = sprintf('%s/well_%d/embryo_%d/processed/offset/',expt,well,embryos(i));
       %offset_dir = sprintf('processed/%s/well%d_embryo%d/offset',expt,well,embryos(i));
        write_3d_tif(sprintf('%s/cy%02d_dapi.tif',offset_dir,cycle),imcrop_xyz(dapi_offset_stacks{cycle+1},embryo_bounds(i,:)));
    end
    
    disp(sprintf('%s: wrote DAPI for cycle %d',sec2time(toc),cycle))
    
    % loop through seq channels
    
    channels = 1:num_channels; channels(channels==dapi_channel) = [];
    for channel=[channels]
        
        % load seq
        
        seq_stack = zeros(xlen(cycle+1),ylen(cycle+1),zlen(cycle+1),'uint16');
        for z=1:zlen(cycle+1)
            seq_stack(:,:,z) = readPlane(seq_reader,series,z,channel,timepoint);
        end
        
        disp(sprintf('%s: loaded cycle %d, channel %d',sec2time(toc),cycle, channel))
        
        % offset seq
        
        offset_seq_stack = apply_offset_xyz(seq_stack,offsets{cycle+1},max_dims);
        figure; imshowpair(capImage(max(dapi_offset_stacks{2},[],3),99,'prc'),capImage(max(offset_seq_stack,[],3),99,'prc'))
        
        disp(sprintf('%s: offset cycle %d, channel %d',sec2time(toc),cycle, channel))
        
        % write seq
        
        for i=1:num_embryos
            offset_dir = sprintf('%s/%s/well_%d/embryo_%d/processed/offset/',home_dir,expt,well,embryos(i));
            %offset_dir = sprintf('processed/%s/well%d_embryo%d/offset',expt,well,embryos(i));
            write_3d_tif(sprintf('%s/cy%02d_ch%02d.tif',offset_dir,cycle,channel),imcrop_xyz(offset_seq_stack,embryo_bounds(i,:)));
        end
        
        disp(sprintf('%s: wrote cycle %d, channel %d',sec2time(toc),cycle, channel))
        
    end
    
end

%% Load full stain stack + DAPI
{
stain_reader = bfGetReader(sprintf('%s/%s/Mechanics, then ExSeq (Jan 2021)/Data/Well %d, Embryo %d/Morphology/*Well_%d_*Embryo_%d*.tif',data_dir,expt,well,embryo,well,embryo));
stain_xlen = stain_reader.getSizeX; stain_ylen = stain_reader.getSizeY; stain_zlen = stain_reader.getSizeZ; 

stain_memb_stack = zeros(stain_xlen,stain_ylen,stain_zlen,'uint16');
stain_dapi_stack = zeros(stain_xlen,stain_ylen,stain_zlen,'uint16');

for z=1:stain_zlen
    stain_memb_stack(:,:,z) = readPlane(stain_reader,series,z,2,timepoint);
    stain_dapi_stack(:,:,z) = readPlane(stain_reader,series,z,3,timepoint);
end

disp(sprintf('%s: loaded stain stack + DAPI',sec2time(toc)))

%% Offset DAPI and apply to stain

[stain_dapi_offset_stack stain_offsets] = get_offset_xyz(dapi_offset_stacks{1},stain_dapi_stack,max_dims);
figure; imshowpair(capImage(max(dapi_offset_stacks{1},[],3),99,'prc'),capImage(max(stain_dapi_stack,[],3),99,'prc'))
disp(sprintf('%s: offset stain DAPI',sec2time(toc)))

offset_stain_memb_stack = apply_offset_xyz(stain_memb_stack,stain_offsets,max_dims);
figure; imshowpair(capImage(max(stain_dapi_offset_stack,[],3),99,'prc'),capImage(max(offset_stain_memb_stack,[],3),99,'prc'))
figure; imshowpair(capImage(max(dapi_offset_stacks{1},[],3),99,'prc'),capImage(max(offset_stain_memb_stack,[],3),99,'prc'))


disp(sprintf('%s: applied offset to stain stack',sec2time(toc)))

offset_dir = sprintf('%s/%s/well_%d/embryo_%d/processed/offset/',home_dir,expt,well,embryos(i));
write_3d_tif(sprintf('%s/stain_dapi.tif',offset_dir),imcrop_xyz(stain_dapi_offset_stack,embryo_bounds(i,:)));
write_3d_tif(sprintf('%s/stain_memb.tif',offset_dir),imcrop_xyz(offset_stain_memb_stack,embryo_bounds(i,:)));

disp(sprintf('%s: saved stain stacks',sec2time(toc)))
}


%%%%%% dapi and membrane stain will be used for segmentation (CellProfiler)
%after saving cellprofiler output into tiff II needed to open and save them in Fiji again

%% Load seq stack

num_channels = 4;
num_cycles = 4;

bounds = dlmread(sprintf('%s/processed/bounds.txt',embryo_dir));

stack = zeros(bounds(5),bounds(4),bounds(6),num_channels,num_cycles,'uint16');

for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/processed/offset/cy%02d_ch%02d.tif',embryo_dir,cycle,channel);
        stack(:,:,:,channel,cycle) = read_3d_tif(filename,bounds(5),bounds(4),bounds(6));
    end
end

% Load hyb + DAPI stacks

%hyb_stack = zeros(bounds(5),bounds(4),bounds(6),'uint16');
%dapi_stack = zeros(bounds(5),bounds(4),bounds(6),'uint16');

hyb_stack = uint16(read_3d_tif(sprintf('%s/processed/offset/hyb_probe.tif',embryo_dir),bounds(5),bounds(4),bounds(6)));
dapi_stack = uint16(read_3d_tif(sprintf('%s/processed/offset/hyb_dapi.tif',embryo_dir),bounds(5),bounds(4),bounds(6)));

figure; imshowpair(capImage(max(hyb_stack,[],3),99,'prc'),capImage(max(dapi_stack,[],3),99,'prc'))

%seg_stack = read_3d_tif(sprintf('%s/processed/offset/seg.tif',embryo_dir),bounds(5),bounds(4),bounds(6));
%stain_stack = read_3d_tif(sprintf('%s/processed/offset/seg_mem.tif',embryo_dir),bounds(5),bounds(4),bounds(6));

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load cropped stain stack + segmentation
%well1_embryo2 

seg_stack = read_3d_tif(sprintf('%s/well_%d/embryo_%d/seg_nuclei_memb/stain_memb_exseq_well%d_embryo%d_segmentaedMemb.tiff',expt,well,embryo,well,embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
seg_mem_stack = read_3d_tif(sprintf('%s/%s/well_%d/embryo_%d/processed/offset/stain_memb.tif',home_dir,expt,well,embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
seg_dapi_stack = read_3d_tif(sprintf('%s/%s/well_%d/embryo_%d/processed/offset/stain_dapi.tif',home_dir,expt,well,embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));

%max_dims = [676, 684, 348];%bounds.txt : 800,790,1,676,684,348
max_dims = [embryo_bounds(4), embryo_bounds(5), embryo_bounds(6)]
[offset_seg_dapi_stack seg_offsets] = get_offset_xyz(dapi_stack,seg_dapi_stack,max_dims);

figure; imshowpair(capImage(max(dapi_stack,[],3),99,'prc'),capImage(max(offset_seg_dapi_stack,[],3),99,'prc'))

offset_seg_stack = apply_offset_xyz(seg_stack,seg_offsets,max_dims);
offset_seg_mem_stack = apply_offset_xyz(seg_mem_stack,seg_offsets,max_dims);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

write_3d_tif(sprintf('%s/seg.tif',offset_dir),imcrop_xyz(offset_seg_stack,embryo_bounds(i,:)));
write_3d_tif(sprintf('%s/seg_dapi.tif',offset_dir),imcrop_xyz(offset_seg_dapi_stack,embryo_bounds(i,:)));
write_3d_tif(sprintf('%s/seg_mem.tif',offset_dir),imcrop_xyz(offset_seg_mem_stack,embryo_bounds(i,:)));
figure; imshowpair(capImage(max(dapi_stack,[],3),99,'prc'),capImage(max(seg_stack,[],3),99,'prc'))

disp(sprintf('%s: Loaded stack',sec2time(toc)))

%% Hyb + DAPI registration

comb_stack = cat(4,dapi_stack,seg_stack,stain_stack);
%comb_stack = dapi_stack;

min_overlap = 0.5;
reg_stack = zeros(size(stack),'uint16');
curr_cycle = zeros(size(dapi_stack),'uint16');
reg_stack(:,:,:,:,1) = stack(:,:,:,:,1);
reg_offsets = {};

% Register hyb to cycle 1

[reg_comb_stack reg_hyb_stack hyb_stats hyb_offset] = offset_cycle(max(reg_stack(:,:,:,:,1),[],4),hyb_stack,comb_stack,min_overlap);
reg_dapi_stack = reg_comb_stack(:,:,:,1);
reg_seg_stack = reg_comb_stack(:,:,:,2);
reg_stain_stack = reg_comb_stack(:,:,:,3);

figure; imshowpair(capImage(max(reg_dapi_stack,[],3),99,'prc'),capImage(max(reg_seg_stack,[],3),99,'prc'))

disp(sprintf('%s: Registered hyb and DAPI to cycle 1: %.02f',sec2time(toc),hyb_stats));

%% Registration

for cycle=2:num_cycles
    
    curr_cycle = max(stack(:,:,:,:,cycle),[],4);  
    [reg_stack(:,:,:,:,cycle) reg_cycle stats reg_offsets{cycle}] = offset_cycle(max(reg_stack(:,:,:,:,1),[],4),curr_cycle,stack(:,:,:,:,cycle),min_overlap);

    disp(sprintf('%s: Cycle %d registration: %.03f',sec2time(toc),cycle,stats));
    
end

%% Color correct cycle 1 to itself (b/c hyb is bad) and save offsets

base_stack = max(reg_stack(:,:,:,[1],1),[],4);

reg_stack2 = reg_stack;

color_stats = {};
color_offsets = {};

cycle = 1;
for channel=2:num_channels
    
    [reg_stack2(:,:,:,channel,cycle) tmp color_stats{cycle,channel} color_offsets{cycle,channel}] = offset_cycle(reg_stack(:,:,:,1,cycle),reg_stack(:,:,:,channel,cycle),reg_stack(:,:,:,channel,cycle),min_overlap);  
    disp(sprintf('%s: Cycle %d, channel %d color correction: %.03f',sec2time(toc),cycle,channel,color_stats{cycle,channel}));
    
end

%%

%figure; imshowpair(capImage(max(max(reg_stack(:,:,:,:,1),[],3),[],4),99,'prc'),capImage(max(max(reg_stack2(:,:,:,:,1),[],3),[],4),99,'prc'))
%figure; imshowpair(capImage(max(max(reg_stack(:,:,:,1,1),[],3),[],4),99,'prc'),capImage(max(max(reg_stack(:,:,:,2,1),[],3),[],4),99,'prc'))
%figure; imshowpair(capImage(max(max(reg_stack2(:,:,:,:,1),[],3),[],4),99,'prc'),capImage(max(max(reg_stack2(:,:,:,2,1),[],3),[],4),99,'prc'))



%% Secondary color correction

color_stats = {};
color_offsets = {};

for cycle=2:num_cycles
    
    for channel=1:num_channels
        [reg_stack2(:,:,:,channel,cycle) tmp color_stats{cycle,channel} color_offsets{cycle,channel}] = offset_cycle(max(reg_stack2(:,:,:,:,1),[],4),reg_stack(:,:,:,channel,cycle),reg_stack(:,:,:,channel,cycle),min_overlap);
        disp(sprintf('%s: Cycle %d, channel %d color correction: %.03f',sec2time(toc),cycle,channel,color_stats{cycle,channel}));
    end

end

%%

figure; imshowpair(capImage(max(max(reg_stack2(:,:,:,2,1),[],3),[],4),99,'prc'),capImage(max(max(reg_stack2(:,:,:,2,2),[],3),[],4),99,'prc'))

%% Channel quantile normalization and zero matrix

norm_stack = zeros(size(reg_stack2));

for cycle=1:num_cycles
    cycle
    tmp_stack = reshape(reg_stack2(:,:,:,:,cycle),bounds(5)*bounds(4)*bounds(6),num_channels);
    quantile_norm = quantilenorm(double(tmp_stack));
    norm_stack(:,:,:,:,cycle) = reshape(quantile_norm,bounds(5),bounds(4),bounds(6),num_channels);
    
    %norm_stack(:,:,:,:,cycle) = norm_stack(:,:,:,:,cycle) - mean(reshape(norm_stack(:,:,:,:,cycle),[],1));
    
end

%norm_stack(norm_stack<0) = 0;

disp(sprintf('%s: Performed quantile normaliation and zeroed matrix',sec2time(toc)))

%% Save visualizations

% All cycles

f = figure('visible','off');
p = tight_subplot(num_channels,num_cycles,[0.001 0.001],[0.001 0.001],[0.001 0.001]);

for cycle=1:num_cycles
    for channel=1:num_channels   
        axes(p((channel-1)*num_cycles+cycle)); imshow(capImage(max(norm_stack(:,:,:,channel,cycle),[],3),99,'prc'),[]);
        title(sprintf('cycle %d, channel %d',cycle,channel))
    end
end

xdim = size(norm_stack,2)/max(size(norm_stack,1),size(norm_stack,2));
ydim = size(norm_stack,1)/max(size(norm_stack,1),size(norm_stack,2));

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 (num_cycles)*2*xdim, num_channels*2.25*ydim];

fig_dir = sprintf('%s/well_%d/embryo_%d/figure/all_stacks_xy',expt,well, embryo);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/embryo%d.png',fig_dir,embryo));
%saveas(fig,sprintf('%s/well%d_embryo%d.png',fig_dir,well,embryo));

%% Save visualizations

% All cycles

f = figure('visible','off');
p = tight_subplot(num_channels,num_cycles,[0.001 0.001],[0.001 0.001],[0.001 0.001]);

for cycle=1:num_cycles
    for channel=1:num_channels   
        axes(p((channel-1)*num_cycles+cycle)); imshow(capImage(squeeze(max(norm_stack(:,:,:,channel,cycle),[],1))',99,'prc'),[]);
        title(sprintf('cycle %d, channel %d',cycle,channel))
    end
end

xdim = size(norm_stack,2)/max(size(norm_stack,2),size(norm_stack,3));
zdim = size(norm_stack,3)/max(size(norm_stack,2),size(norm_stack,3));

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 (num_cycles)*2*xdim, num_channels*3*zdim];

fig_dir = sprintf('%s/well_%d/embryo_%d/figure/all_stacks_xz',expt,well, embryo);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/embryo%d.png',fig_dir,embryo));
%saveas(fig,sprintf('%s/well%d_embryo%d.png',fig_dir,well,embryo));

%% Cycle overlap

f = figure('visible','off');
p = tight_subplot(2,num_cycles,[0.001 0.001],[0.001 0.001],[0.001 0.001]);

for cycle=1:num_cycles
     
    if cycle == 1
        init_corr(cycle) = corr(double(reshape(hyb_stack,[],1)),double(reshape(max(stack(:,:,:,:,cycle),[],4),[],1)));
        reg_corr(cycle) = corr(double(reshape(hyb_stack,[],1)),reshape(max(norm_stack(:,:,:,:,cycle),[],4),[],1));
        
        axes(p(cycle)); imshowpair(capImage(max(hyb_stack,[],3),95,'prc'),capImage(max(max(stack(:,:,:,:,cycle),[],4),[],3),95,'prc')); 
        title(sprintf('Hyb - Cycle %d: %.03f',cycle,init_corr(cycle)));
        
        axes(p(num_cycles+cycle)); imshowpair(capImage(max(hyb_stack,[],3),95,'prc'),capImage(max(max(norm_stack(:,:,:,:,cycle),[],4),[],3),95,'prc'));
        title(sprintf('Hyb - Cycle %d: %.03f',cycle,reg_corr(cycle)))
        
    else
        init_corr(cycle) = corr(double(reshape(max(stack(:,:,:,:,cycle-1),[],4),[],1)),double(reshape(max(stack(:,:,:,:,cycle),[],4),[],1)));
        reg_corr(cycle) = corr(double(reshape(max(norm_stack(:,:,:,:,cycle-1),[],4),[],1)),double(reshape(max(norm_stack(:,:,:,:,cycle),[],4),[],1))); 
        
        axes(p(cycle)); imshowpair(capImage(max(max(stack(:,:,:,:,cycle-1),[],4),[],3),95,'prc'),capImage(max(max(stack(:,:,:,:,cycle),[],4),[],3),95,'prc')); 
        title(sprintf('Cycle %d - %d: %.03f',cycle-1,cycle,init_corr(cycle)));
        
        axes(p(num_cycles+cycle)); imshowpair(capImage(max(max(norm_stack(:,:,:,:,cycle-1),[],4),[],3),95,'prc'),capImage(max(max(norm_stack(:,:,:,:,cycle),[],4),[],3),95,'prc'));
        title(sprintf('Cycle %d - %d: %.03f',cycle-1,cycle,reg_corr(cycle)))
         
    end
   
end

xdim = size(norm_stack,2)/max(size(norm_stack,1),size(norm_stack,2));
ydim = size(norm_stack,1)/max(size(norm_stack,1),size(norm_stack,2));

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 (num_cycles)*2*xdim, (2)*2.25*ydim];
fig_dir = sprintf('%s/well_%d/embryo_%d/figure/cycle_overlap_xy',expt, well, embryo);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/embryo%d.png',fig_dir,embryo));
%saveas(fig,sprintf('%s/well%d_embryo%d.png',fig_dir,well,embryo));

disp(sprintf('%s: Saved visualizations',sec2time(toc)))

%% Cycle overlap

f = figure('visible','off');
p = tight_subplot(2,num_cycles,[0.001 0.001],[0.001 0.001],[0.001 0.001]);

for cycle=1:num_cycles
    
    if cycle == 1
        init_corr(cycle) = corr(double(reshape(hyb_stack,[],1)),double(reshape(max(stack(:,:,:,:,cycle),[],4),[],1)));
        reg_corr(cycle) = corr(double(reshape(hyb_stack,[],1)),reshape(max(norm_stack(:,:,:,:,cycle),[],4),[],1));
        
        axes(p(cycle)); imshowpair(...
            capImage(squeeze(max(hyb_stack,[],1)),95,'prc')',...
            capImage(squeeze(max(max(stack(:,:,:,:,cycle),[],4),[],1)),95,'prc')'...
            ); 
        title(sprintf('Hyb - Cycle %d: %.03f',cycle,init_corr(cycle)));
        
        axes(p(num_cycles+cycle)); imshowpair(...
            capImage(squeeze(max(hyb_stack,[],1)),95,'prc')',...
            capImage(squeeze(max(max(norm_stack(:,:,:,:,cycle),[],4),[],1)),95,'prc')'...
            );
        title(sprintf('Hyb - Cycle %d: %.03f',cycle,reg_corr(cycle)))
        
    else
        init_corr(cycle) = corr(double(reshape(max(stack(:,:,:,:,cycle-1),[],4),[],1)),double(reshape(max(stack(:,:,:,:,cycle),[],4),[],1)));
        reg_corr(cycle) = corr(double(reshape(max(norm_stack(:,:,:,:,cycle-1),[],4),[],1)),double(reshape(max(norm_stack(:,:,:,:,cycle),[],4),[],1))); 
        
        axes(p(cycle)); imshowpair(...
            capImage(squeeze(max(max(stack(:,:,:,:,cycle-1),[],4),[],1)),95,'prc')',...
            capImage(squeeze(max(max(stack(:,:,:,:,cycle),[],4),[],1)),95,'prc')'...
            ); 
        title(sprintf('Cycle %d - %d: %.03f',cycle-1,cycle,init_corr(cycle)));
        
        axes(p(num_cycles+cycle)); imshowpair(...
            capImage(squeeze(max(max(norm_stack(:,:,:,:,cycle-1),[],4),[],1)),95,'prc')',...
            capImage(squeeze(max(max(norm_stack(:,:,:,:,cycle),[],4),[],1)),95,'prc')'...
            );
        title(sprintf('Cycle %d - %d: %.03f',cycle-1,cycle,reg_corr(cycle)))
         
    end
   
end

xdim = size(norm_stack,2)/max(size(norm_stack,2),size(norm_stack,3));
zdim = size(norm_stack,3)/max(size(norm_stack,2),size(norm_stack,3));

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 (num_cycles)*2*xdim, (2)*3*zdim];


fig_dir = sprintf('%s/well_%d/embryo_%d/figure/cycle_overlap_xz',expt, well, embryo);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(fig,sprintf('%s/embryo%d.png',fig_dir,embryo));
%saveas(fig,sprintf('%s/well%d_embryo%d.png',fig_dir,well,embryo));

disp(sprintf('%s: Saved visualizations',sec2time(toc)))

%% Save stacks

% Save normalized stack

for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/cy%02d_ch%02d.tif',reg_dir,cycle,channel);
        write_3d_tif(filename, norm_stack(:,:,:,channel,cycle));
    end
end

% Save hyb + DAPI and segmentation

write_3d_tif(sprintf('%s/hyb_dapi.tif',reg_dir),reg_dapi_stack);
write_3d_tif(sprintf('%s/hyb_probe.tif',reg_dir),reg_hyb_stack);
write_3d_tif(sprintf('%s/seg.tif',reg_dir),reg_seg_stack);
write_3d_tif(sprintf('%s/seg_mem.tif',reg_dir),reg_stain_stack);

disp(sprintf('%s: Saved normalized stacks',sec2time(toc)))


%% Load stain segmentation
bounds = dlmread(sprintf('%s/well_%d/embryo_%d/processed/bounds.txt',expt, well,embryo));
stain_seg_stack = read_3d_tif(sprintf('%s/well_%d/embryo_%d/seg_nuclei_memb/stain_memb_exseq_well%d_embryo%d_segmentaedMemb.tiff',expt,well,embryo,well,embryo),bounds(5),bounds(4),bounds(6));

%% Register stain segmentation to stain

%reg_stain_seg_stack = get_offset_xyz(reg_dapi_stack,double(stain_seg_stack>0),size(stain_seg_stack));
%crop_stain_seg_stack = imcrop_xyz(reg_stain_seg_stack,[1 1 1 size(reg_dapi_stack,1) size(reg_dapi_stack,2) 225]);

%figure; imshowpair(max(crop_stain_seg_stack,[],3)>0,capImage(max(reg_dapi_stack,[],3),95,'prc'))
%%%

figure; for z=1:size(stain_stack,3)
    imshow(capImage(stain_stack(:,:,z),99,'prc'),[])
end

%figure; 
%for z=1:size(stain_seg_stack,3)
%    imshow(reg_stain_seg_stack(:,:,z),[])
%end

misc_dir = sprintf('%s/well_%d/embryo_%d/figure/misc',expt,well,embryo);
if ~exist(misc_dir, 'dir') mkdir(misc_dir), end

video = VideoWriter(sprintf('%s/well_%d/embryo_%d/figure/misc/embryo%d_stain_comp2',expt,well,embryo,embryo),'MPEG-4');
video.FrameRate = 2.5;
open(video);

figure; 
for z=1:size(stain_stack,3)
    imshowpair(capImage(reg_dapi_stack(:,:,z),99,'prc'),stain_seg_stack(:,:,z)>0)
    writeVideo(video,getframe(gcf))
end

close(video)



%% Load stacks

bounds = dlmread(sprintf('%s/processed/bounds.txt',embryo_dir));
norm_stack = zeros(bounds(5),bounds(4),bounds(6),num_channels,num_cycles);

% Load normalized stack

for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/reg/cy%02d_ch%02d.tif',embryo_dir,cycle,channel);
        norm_stack(:,:,:,channel,cycle) = read_3d_tif(filename,bounds(5),bounds(4),bounds(6));
    end
end

% Load hyb + DAPI and segmentation

reg_hyb_stack = read_3d_tif(sprintf('%s/reg/hyb_probe.tif',embryo_dir),bounds(5),bounds(4),bounds(6));
reg_dapi_stack = read_3d_tif(sprintf('%s/reg/hyb_dapi.tif',embryo_dir),bounds(5),bounds(4),bounds(6));

disp(sprintf('%s: Loaded normalized stacks',sec2time(toc)))

%% Deconvolve with Gaussian filter

deconv_stack = zeros(size(norm_stack));

for cycle=1:num_cycles
    for channel=1:num_channels
        for z=1:size(deconv_stack,3)
            high_pass_filter = imgaussfilt(norm_stack(:,:,z,channel,cycle),2);
            high_pass_image = norm_stack(:,:,z,channel,cycle) - high_pass_filter;
            deconv_stack(:,:,z,channel,cycle) = high_pass_image;
        end
    end
end

deconv_stack(deconv_stack<0) = 0;

disp(sprintf('%s: Deconvolved images with Gaussian filter',sec2time(toc)))



%% Peak calling

all_peaks = [];
channel_peaks = {};
cycle = 1;

for channel=1:num_channels
        
    [Maxima,MaxPos,Minima,MinPos]=MinimaMaxima3D(deconv_stack(:,:,:,channel,cycle),1,0);
    
    channel_peaks_pos{channel} = MaxPos;
    channel_peaks_max{channel} = Maxima;

    disp(sprintf('%s: Found %d 3D peaks in cycle %d, channel %d',sec2time(toc),length(channel_peaks_pos{channel}),cycle,channel));
          
end

peaks_table = table;
peaks_table.x = [channel_peaks_pos{1}(:,1); channel_peaks_pos{2}(:,1); channel_peaks_pos{3}(:,1); channel_peaks_pos{4}(:,1)];
peaks_table.y = [channel_peaks_pos{1}(:,2); channel_peaks_pos{2}(:,2); channel_peaks_pos{3}(:,2); channel_peaks_pos{4}(:,2)];
peaks_table.z = [channel_peaks_pos{1}(:,3); channel_peaks_pos{2}(:,3); channel_peaks_pos{3}(:,3); channel_peaks_pos{4}(:,3)];
peaks_table.val = [channel_peaks_max{1}; channel_peaks_max{2}; channel_peaks_max{3}; channel_peaks_max{4}];
peaks_table.channel = [repmat(1,length(channel_peaks_max{1}),1); repmat(2,length(channel_peaks_max{2}),1); repmat(3,length(channel_peaks_max{3}),1); repmat(4,length(channel_peaks_max{4}),1);];

writetable(peaks_table,sprintf('%s/peaks.txt',embryo_dir))

%% Read peaks

peaks_table = readtable(sprintf('%s/peaks.txt',embryo_dir));

head(peaks_table)

ch1 = (peaks_table(peaks_table.val == '1',:));


%% Peak call video

out_dir = sprintf('%s/well_%d/embryo_%d/figure/peak_calls',expt,well,embryo);
if ~exist(out_dir, 'dir') mkdir(out_dir), end

video = VideoWriter(sprintf('%s/well%d_embryo%d',out_dir,well,embryo),'MPEG-4');
video.FrameRate = 2.5;
open(video);

z_buffer = 2;

%thresh = [150 150 150 200]; well1_embryo3
thresh = [200 250 300 400];%well2_embryo1

for z=1:size(norm_stack,3)
    
    f = figure('visible','off');
    p = tight_subplot(1,4,[0.001 0.001],[0.001 0.001],[0.001 0.001]);
    
    for channel=1:num_channels
    
        peaks = peaks_table{peaks_table.channel == channel & peaks_table.val > thresh(channel),1:3};
        %peaks = channel_peaks_pos{channel}(channel_peaks_max{channel}>thresh(channel),:);
        z_channel_peaks = peaks(ismember(peaks(:,3),(z+[-z_buffer:+z_buffer])),:);

        axes(p(channel));
        imshowpair(capImage(deconv_stack(:,:,z,channel,1),99,'prc'),capImage(reg_dapi_stack(:,:,z),95,'prc')); hold on;
        plot(z_channel_peaks(:,2),z_channel_peaks(:,1),'Marker','.','MarkerEdgeColor','r','LineStyle','none');
        title(sprintf('Peak calls (n = %d)',size(peaks,1)))
        
    end
   
    %fig.Units = 'pixels';
    f.Position = [0 0 size(norm_stack,2).*0.5.*(4) size(norm_stack,1).*0.75];    
    writeVideo(video,getframe(gcf))
    
end
close(video); %close the file
close all

disp(sprintf('%s: Saved peak call video',sec2time(toc)));
%}
%%

%thresh = [100 100 100 150];%; well1_embryo3

thresh = [100 100 100 150];%; well2_embryo1
all_peaks = [];
for channel=1:num_channels
    all_peaks = [all_peaks; peaks_table{peaks_table.channel == channel & peaks_table.val > thresh(channel),1:3}];
end

pad = 2;
mat = zeros(size(all_peaks,1),num_channels,num_cycles);
purity = zeros(size(all_peaks,1),4);
consensus = zeros(size(all_peaks,1),4);

xleft = max(all_peaks(:,1)-pad,1); xright = min(all_peaks(:,1)+pad,size(deconv_stack,1));
yleft = max(all_peaks(:,2)-pad,1); yright = min(all_peaks(:,2)+pad,size(deconv_stack,2));
zleft = max(all_peaks(:,3)-pad,1); zright = min(all_peaks(:,3)+pad,size(deconv_stack,3));

for i=1:size(all_peaks,1); disp(i)
    
    peak = all_peaks(i,:);
    peak_mat = squeeze(sum(sum(sum(deconv_stack(xleft(i):xright(i),yleft(i):yright(i),zleft(i):zright(i),:,:),1),2),3));
    
    mat(i,:,:) = peak_mat;
    purity(i,:) = max(peak_mat.^2,[],1)./sum(peak_mat.^2,1);
    [tmp consensus(i,:)] = max(peak_mat,[],1);
    
end

    

%%
{
thresh = 0.75;
purity_score = sum(purity>thresh,2);

figure; histogram(sum(purity>thresh,2))
xlabel('rounds w/ purity score > 0.75'); xticks([0:4])
ylabel('spots')
}
%%
{
out_dir = sprintf('figures/%s/purity_video',expt);
if ~exist(out_dir, 'dir') mkdir(out_dir), end

video = VideoWriter(sprintf('%s/embryo%d',out_dir,embryo),'MPEG-4');
video.FrameRate = 2.5;
open(video);

z_buffer = 2;

for z=1:size(norm_stack,3)
    
    f = figure;
    
    z_peaks = all_peaks(ismember(all_peaks(:,3),(z+[-z_buffer:+z_buffer])),:);
    z_purity = purity_score(ismember(all_peaks(:,3),(z+[-z_buffer:+z_buffer])));
        
    imshowpair(capImage(max(deconv_stack(:,:,z,channel,1),[],4),99,'prc'),capImage(reg_dapi_stack(:,:,z),95,'prc')); hold on;
    scatter(z_peaks(:,2),z_peaks(:,1),15,z_purity,'filled');
    colormap(jet); caxis([0 4])

    %fig.Units = 'pixels';
    %f.Position = [0 0 size(deconv_stack,2) size(deconv_stack,1)];    
    writeVideo(video,getframe(gcf))
    
end
close(video); %close the file
close all
%}

%% Add in barcode info

valid_barcodes = [...
    "0012212","3023023","3011033","2010130",...
    "0232030","3302023","3122303","2112010","3310103","0102320","0110113","1103130","3231202",...
    "2021200","1221220","0133213","0031203",...
    "2033212","2311012","1312121","2002211","0322110",...
    "3032231","0330132","3101021","1120321","1331011","3203300","0301003","2123022","1022102","2131002",...
    "1132331","0020020","3321221","0211210","3212321","1202120","2303131","3220033","2332113","1300300",...
    "2320210","0223310","1233013","3113223","1030202","0313311","2201322"...
    ]';

tmp = char(valid_barcodes);
valid_barcodes = string([tmp(:,1) tmp(:,2) tmp(:,3) tmp(:,4) ]);

genes = [...
    "Gata3","Krt18","Krt8","Id2",...
    "Lgals1","Sox17","Cubn","Foxq1","Sparc","Col4a1","Cdx2","Nanog","Bmyc",...
    "Uap1","Eomes","Gata6","Spic",...
    "Amot","Eno1","Sox2","Klf4","Gsc",...
    "Bhmt","Jam2","Tdgf1","Elf3","Fgd1","Satb1","Morc1","Pou5f1","Pak1","Fgf4",...
    "Tead4","Fgfr2","Otx2","Apoe","Lama1","Lrp2","Fgfr1","Fn1","Slc7a3","Fbxo15",...
    "Tdh","Dppa2","Dkk1","B3gnt5","B3galnt1","St6galnac4","Acaa2_1"...
    ];

[sort_valid_barcodes sort_order] = sort(valid_barcodes);
sort_genes = genes(sort_order);

%%

sel_consensus = consensus(purity_score>=4,:);
barcodes = strrep(string(num2str(4-sel_consensus)),' ','');
[C, ia, ic] = unique(barcodes);
counts = accumarray(ic,1);

on_target = ismember(C,sort_valid_barcodes);
on_target_prc = sum(counts(on_target))./sum(counts);
on_target_genes = sort_genes(ismember(sort_valid_barcodes,C));

figure;
hold on
for i = 1:length(C)
    h=bar(i,counts(i));
    if on_target(i) == 1
        set(h,'FaceColor','r');
    else
        set(h,'FaceColor','b');
    end
end
hold off

xlim([1 size(C,1)])

%xticks([find(counts>sum(counts)/256*3)])
%xticklabels(C(counts>sum(counts)/256*3))

xticks(find(on_target))
xticklabels(strcat(on_target_genes'," - ",C(on_target)," (",string(counts(on_target)),")"))
set(gca,'fontsize',8)

%xticks(find(~on_target & counts>sum(counts)/256*3))
%xticklabels(strcat("N/A (",C(on_target),")"))

xtickangle(90)
title(sprintf("stringent on target %%: %.02f",on_target_prc*100))

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 5];

fig_dir = sprintf('%s/well_%d/embryo_%d/figure/stringent_on_target',expt,well,embryo);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

%saveas(fig,sprintf('%s/embryo%d.png',fig_dir,embryo));
saveas(fig,sprintf('%s/well%d_embryo%d.png',fig_dir,well,embryo));

genes_stringent = on_target_genes;
counts_stringent = counts(on_target);

%%

sel_consensus = consensus(purity_score>=3,:);
barcodes = strrep(string(num2str(4-sel_consensus)),' ','');
[C, ia, ic] = unique(barcodes);
counts = accumarray(ic,1);

on_target = ismember(C,sort_valid_barcodes);
on_target_prc = sum(counts(on_target))./sum(counts);
on_target_genes = sort_genes(ismember(sort_valid_barcodes,C));

figure;
hold on
for i = 1:length(C)
    h=bar(i,counts(i));
    if on_target(i) == 1
        set(h,'FaceColor','r');
    else
        set(h,'FaceColor','b');
    end
end
hold off

xlim([1 size(C,1)])

%xticks([find(counts>sum(counts)/256*3)])
%xticklabels(C(counts>sum(counts)/256*3))

xticks(find(on_target))
xticklabels(strcat(on_target_genes'," - ",C(on_target)," (",string(counts(on_target)),")"))
set(gca,'fontsize',8)

%xticks(find(~on_target & counts>sum(counts)/256*3))
%xticklabels(strcat("N/A (",C(on_target),")"))

xtickangle(90)
title(sprintf("lenient on target %%: %.02f",on_target_prc*100))

fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 10 5];

fig_dir = sprintf('%s/well_%d/embryo_%d/figure/lenient_on_target',expt,well,embryo);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

%saveas(fig,sprintf('%s/embryo%d.png',fig_dir,embryo));
saveas(fig,sprintf('%s/well%d_embryo%d.png',fig_dir,well,embryo));

genes_lenient = on_target_genes;
counts_lenient = counts(on_target);

%% Save gene counts

%on_target_genes = sort_genes(ismember(sort_valid_barcodes,C));

gene_counts = table;
gene_counts.gene = sort_genes';

gene_counts.counts_stringent = zeros(size(gene_counts,1),1);
gene_counts.counts_stringent(ismember(sort_genes',genes_stringent')) = counts_stringent;

gene_counts.counts_lenient = zeros(size(gene_counts,1),1);
gene_counts.counts_lenient(ismember(sort_genes',genes_lenient')) = counts_lenient;

fig_dir = sprintf('%s/well_%d/embryo_%d/figure/gene_counts_table',expt,well,embryo);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

writetable(gene_counts,sprintf('%s/gene_counts_well%d_embryo%d.txt',fig_dir,well,embryo))

%% 

sel_consensus = consensus(purity_score>=3,:);
barcodes = strrep(string(num2str(4-sel_consensus)),' ','');

sel_purity = purity_score(purity_score>=3);
on_target_purity = sel_purity(ismember(barcodes,sort_valid_barcodes));

sel_peaks = all_peaks(purity_score>=3,:);
on_target_peaks = sel_peaks(ismember(barcodes,sort_valid_barcodes),:);

barcode_dict = strings(4000,1);
barcode_dict(str2double(sort_valid_barcodes)) = sort_genes';
on_target_genes = barcode_dict(str2double(barcodes(ismember(barcodes,sort_valid_barcodes))));

gene_xyz_table = table;
gene_xyz_table.gene = on_target_genes;
gene_xyz_table.x = on_target_peaks(:,1);
gene_xyz_table.y = on_target_peaks(:,2);
gene_xyz_table.z = on_target_peaks(:,3);
gene_xyz_table.purity = on_target_purity;

fig_dir = sprintf('%s/well_%d/embryo_%d/figure/gene_xyz_table',expt,well,embryo);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

writetable(gene_xyz_table,sprintf('%s/gene_xyz_table_well%d_embryo%d.txt',fig_dir,well,embryo))

%% Load segmentation

seg_stack = read_3d_tif(sprintf('%s/reg/seg.tif',embryo_dir),bounds(5),bounds(4),bounds(6));
seg_stain_stack = read_3d_tif(sprintf('%s/reg/seg_mem.tif',embryo_dir),bounds(5),bounds(4),bounds(6));

figure; imshowpair(capImage(max(seg_stack,[],3),99,'prc'),capImage(max(seg_stain_stack,[],3),99,'prc'))

%reg_seg_stack
%reg_dapi_stack
%reg_stain_stack
%% Get assignments

sel = gene_xyz_table.purity >= 4; % 3 pr 4

gene_idx = sub2ind(size(seg_stack),gene_xyz_table.x(sel),gene_xyz_table.y(sel),gene_xyz_table.z(sel));
cell_assignment = seg_stack(gene_idx);

%% transcripts per cell

assigned = cell_assignment(cell_assignment ~= 0);
[C ia ic] = unique(assigned);
counts = accumarray(ic,1);

num_cells = max(seg_stack(:));
genes_per_cell = zeros(num_cells,1);
genes_per_cell(C) = counts;

%figure; bar(genes_per_cell); hold on;

colors = distinguishable_colors(101); colors(4,:) = [];
fig = figure; hold on;
for i = 1:length(genes_per_cell)
    h=bar(i,genes_per_cell(i));
    set(h,'FaceColor',colors(i,:));
end

xlim([0 num_cells])
xlabel('cell index'); ylabel('# transcripts')
title(sprintf('%.2f%% of on-target transcripts assigned',sum(cell_assignment>0)./size(cell_assignment,1)*100))

saveas(fig,sprintf('%s/well_%d/embryo_%d/figure/misc/transcripts_per_cell_embryo3.png',expt,well,embryo));

%% create matrix

sel = gene_xyz_table.purity >= 4; % 3 pr 4
num_genes = size(valid_barcodes,1);
gene_xyz_table.cell_assignment = zeros(size(gene_xyz_table,1),1);
gene_xyz_table.cell_assignment(sel) = cell_assignment;

cell_counts_mat = zeros(num_cells,num_genes);

for i=1:num_genes
    
    sel_cells = gene_xyz_table.cell_assignment(gene_xyz_table.gene == sort_genes(:,i));
    [C ia ic] = unique(sel_cells(sel_cells>0));
    counts = accumarray(ic,1);
    cell_counts_mat(C,i) = counts;
end

cell_counts_table = array2table(cell_counts_mat);
cell_counts_table.Properties.VariableNames = cellstr(sort_genes');
cell_counts_table.cell_index = (1:size(cell_counts_table,1))';

cell_info = regionprops(seg_stack,'Centroid','Area');
centroids = reshape([cell_info.Centroid]',3,[])';
cell_counts_table.x = centroids(:,2)*.17;
cell_counts_table.y = centroids(:,1)*.17;
cell_counts_table.z = centroids(:,3)*.4;
cell_counts_table.area = [cell_info.Area]';

cell_counts_table = [cell_counts_table(:,end-4:end) cell_counts_table(:,1:end-5)]; 
%cell_counts_table = cell_counts_table(cell_counts_table.area > 100,:);

fig_dir = sprintf('%s/well_%d/embryo_%d/figure/cell_counts_table',expt, well, embryo);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

writetable(cell_counts_table,sprintf('%s/cell_counts_table_well%d_embryo%d.txt',fig_dir,well,embryo))

%% visualiztion

video = VideoWriter(sprintf('%s/figure/misc/embryo%d_cell_assignments',embryo_dir,embryo),'MPEG-4');
video.FrameRate = 5;
open(video);

zbuf = 1;

fig = figure; 
for z=1:size(reg_dapi_stack,3)
    
    %fig = figure;
    imshow(capImage(reg_dapi_stack(:,:,z),99,'prc'),[]); hold on;
    
    z_peaks = gene_xyz_table(any(gene_xyz_table.z == z+[-zbuf:zbuf],2) & gene_xyz_table.cell_assignment > 0,:);
    scatter(z_peaks.y,z_peaks.x,15,colors(z_peaks.cell_assignment,:),'filled');

    for i=1:num_cells
        [tmp L] = bwboundaries(seg_stack(:,:,z)==i,8,'noholes');
        if size(tmp,1)>0
            B{i} = tmp{1};
        end
    end

    for k = 1:length(B)
        if size(B{k},1)>0
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1),'Color',colors(k,:), 'LineWidth', 1); hold on;
        end
    end
    
    fig.Units = 'pixels';
    fig.Position = [0 0 size(reg_dapi_stack,2) size(reg_dapi_stack,1)];    
    
    writeVideo(video,getframe(gcf))
end

close(video)


%% visualize cell centroid positions

figure; 
scatter3(cell_counts_table.x,cell_counts_table.y,cell_counts_table.z,50,colors(1:size(cell_counts_table,1),:),'filled'); hold on;
sel = gene_xyz_table.cell_assignment > 0;
scatter3(gene_xyz_table.x(sel)*.17,gene_xyz_table.y(sel)*.17,gene_xyz_table.z(sel)*.4,5,colors(gene_xyz_table.cell_assignment(sel),:),'.');
%scatter3(gene_xyz_table.x(~sel)*.17,gene_xyz_table.y(~sel)*.17,gene_xyz_table.z(~sel)*.4,5,[.7 .7 .7],'.');

axis equal;

%% visualize cell centroid positions

sel_gene = 'Sox17';
sel_exp = cell_counts_table{:,find(sort_genes==sel_gene)+5};

figure; 
scatter3(cell_counts_table.x,cell_counts_table.y,cell_counts_table.z,50,sel_exp,'filled'); hold on;
colorbar; %caxis([prctile(sel_exp,10) prctile(sel_exp,90)])
title(sel_gene)

sel = gene_xyz_table.gene == sel_gene & gene_xyz_table.cell_assignment>0;
scatter3(gene_xyz_table.x(sel)*.17,gene_xyz_table.y(sel)*.17,gene_xyz_table.z(sel)*.4,10,sel_exp(gene_xyz_table.cell_assignment(sel)),'.');
%scatter3(gene_xyz_table.x(~sel)*.17,gene_xyz_table.y(~sel)*.17,gene_xyz_table.z(~sel)*.4,5,[.7 .7 .7],'.');

axis equal;

%% visualiztion

sel_gene = 'Nanog';
sel_exp = cell_counts_table{:,find(sort_genes==sel_gene)+5};
cm = parula(max(sel_exp)-min(sel_exp)+1);

video = VideoWriter(sprintf('figures/%s/misc/embryo3_Fbxo15_exp',expt),'MPEG-4');
video.FrameRate = 5;
open(video);

zbuf = 1;

fig = figure; 
for z=1:size(reg_dapi_stack,3)
    
    %fig = figure;
    imshow(capImage(reg_dapi_stack(:,:,z),99,'prc'),[]); hold on;
    
    z_peaks = gene_xyz_table(any(gene_xyz_table.z == z+[-zbuf:zbuf],2) & gene_xyz_table.cell_assignment > 0 & gene_xyz_table.gene == sel_gene,:);
    scatter(z_peaks.y,z_peaks.x,15,cm(1+sel_exp(z_peaks.cell_assignment),:),'filled');

    for i=1:num_cells
        [tmp L] = bwboundaries(seg_stack(:,:,z)==i,8,'noholes');
        if size(tmp,1)>0
            B{i} = tmp{1};
        end
    end

    for k = 1:length(B)
        if size(B{k},1)>0
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1),'Color',cm(1+sel_exp(k),:), 'LineWidth', 1); hold on;
        end
    end
    
    fig.Units = 'pixels';
    fig.Position = [0 0 size(reg_dapi_stack,2) size(reg_dapi_stack,1)];    
    
    writeVideo(video,getframe(gcf))
end

close(video)

%% calculate convex hulls

for i=1:num_cells; disp(i)
    if cell_info(i).Area > 100
        [x y z] = ind2sub(size(seg_stack),find(seg_stack==i));
        hulls{i} = convhull(x*.17,y*.17,z*.4);
    end
end

%%

figure; 
for i=1:num_cells; disp(i)
    if cell_info(i).Area > 100
        [x y z] = ind2sub(size(seg_stack),find(seg_stack==i));
        k1 = hulls{i};
        trisurf(k1,x*.17,y*.17,z*.4,'FaceColor',colors(i,:),'FaceAlpha',0.25,'LineStyle','none'); hold on;
    end
end

axis equal;
%%

%sel_gene = 'Fbxo15';
for sel_gene=[sort_genes]

sel_exp = cell_counts_table{:,find(sort_genes==sel_gene)+5};
cm = parula(max(sel_exp)-min(sel_exp)+1);

figure; 
for i=1:num_cells; disp(i)
    if cell_info(i).Area > 100
        [x y z] = ind2sub(size(seg_stack),find(seg_stack==i));
        k1 = hulls{i};
        trisurf(k1,x*.17,y*.17,z*.4,'FaceColor',cm(1+sel_exp(i),:),'FaceAlpha',0.25,'LineStyle','none'); hold on;
    end
end

title(sel_gene)
colormap(parula); colorbar; caxis([min(sel_exp) max(sel_exp)])

axis equal;

fig_dir = sprintf('%s/well_%d/embryo_%d/figure/gene_exp_convex',expt);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

saveas(gcf,sprintf('%s/embryo3_%s.png',fig_dir,sel_gene))

end
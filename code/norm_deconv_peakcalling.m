% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Embryo Ex-Seq (non-ExM) registration of FoV

function[] = norm_deconv_peakcalling(data_dir, home_dir, expt,plate,well,embryo, bounderies)

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
disp(sprintf('%s ----- quantile normaliation ...',time));

%%
embryo_dir = sprintf('%s/plate_%d_well_%d/embryo_%d',expt,plate,well,embryo);
p_w_e = sprintf('plate_%d_well_%d_embryo_%d', plate, well, embryo);
embryo_bounds = dlmread(sprintf('%s/processed/bounds.txt',embryo_dir));
%embryo_bounds = bounderies(p_w_e).hyb;

reg_dir = sprintf('%s/reg',embryo_dir);
if ~exist(reg_dir, 'dir') mkdir(reg_dir), end

%% Load seq stack

num_channels = 4;
num_cycles = 7;

%embryo_bounds = dlmread(sprintf('%s/processed/bounds.txt',embryo_dir));
stack = zeros(embryo_bounds(5),embryo_bounds(4),embryo_bounds(6),num_channels,num_cycles,'uint16');

for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/processed/offset/cy%02d_ch%02d.tif',embryo_dir,cycle,channel);
        stack(:,:,:,channel,cycle) = read_3d_tif(filename,embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
    end
end

hyb_dapi = uint16(read_3d_tif(sprintf('%s/processed/offset/hyb_dapi.tif',embryo_dir),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6)));
hyb_stack = uint16(read_3d_tif(sprintf('%s/processed/offset/hyb_probe.tif',embryo_dir),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6)));

figure; imshow(capImage(max(hyb_dapi,[],3),99,'prc'))
%% 
%so we want to register in Z for channels 1,3,4 and X,Y,Z for channel2 -
%Channel 2's offset in xy should be consistent. 

% Register hyb to cycle 1
reg_offsets = {};
max_dims = [embryo_bounds(4),embryo_bounds(5),embryo_bounds(6)]
for cyc = 1:num_cycles
    for channel = 1:4
        [reg_stack_temp reg_offsets{cyc,channel}] = get_offset_xyz(hyb_stack,stack(:,:,:,channel,cyc),max_dims);
         tempoffset = reg_offsets{cyc,channel};

    end
end

offset_ch2 = median(vertcat(reg_offsets{:,2}));
for cyc = 1:num_cycles
    
    for channel = 1:4
    tempoffset = reg_offsets{cyc,channel};
        if channel == 2
              tempoffset(1) = offset_ch2(1); 
             tempoffset(2) = offset_ch2(2); 
         else
             tempoffset(1) = 0; 
             tempoffset(2) = 0; 
         end
        reg_stack_temp = apply_offset_xyz(stack(:,:,:,channel,cyc),tempoffset,max_dims);
        reg_stack(:,:,:,channel,cyc) = reg_stack_temp;


    end
end


figure; imshowpair(capImage(max(max(reg_stack(:,:,:,:,1),[],4),[],3),99,'prc'),capImage(max(hyb_stack,[],3),99,'prc'))
figure; imshowpair(capImage(max(max(reg_stack(:,:,:,:,1),[],4),[],3),99,'prc'),capImage(max(max(reg_stack(:,:,:,:,4),[],4),[],3),99,'prc'))
      
%disp(sprintf('%s: Registered all channels to hyb: %.02f',sec2time(toc),hyb_stats));
disp(sprintf('%s: Registered all channels to hyb %.02f',sec2time(toc)));
%% Channel quantile normalization and zero matrix
clear stack;

norm_stack = zeros(size(reg_stack));

%norm  correction
for cycle=1:num_cycles
    cycle
    %tmp_stack = reshape(reg_stack(:,:,:,:,cycle),embryo_bounds(5)*embryo_bounds(4)*embryo_bounds(6),num_channels-1);
    tmp_stack = reshape(reg_stack(:,:,:,:,cycle),embryo_bounds(5)*embryo_bounds(4)*embryo_bounds(6),num_channels);
    quantile_norm = quantilenorm(double(tmp_stack));
    %norm_stack(:,:,:,:,cycle) = reshape(quantile_norm,embryo_bounds(5),embryo_bounds(4),embryo_bounds(6),num_channels-1);
    norm_stack(:,:,:,:,cycle) = reshape(quantile_norm,embryo_bounds(5),embryo_bounds(4),embryo_bounds(6),num_channels);
        
end

disp(sprintf('%s: Performed quantile normaliation and zeroed matrix',sec2time(toc)))


%% Save stacks
disp(sprintf('%s: Saving quantile normaliation and zeroed matrix',sec2time(toc)))
% Save normalized stack with color correction

for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/cy%02d_ch%02d.tif',reg_dir,cycle,channel);
        write_3d_tif(filename, norm_stack(:,:,:,channel,cycle));
    end
end

% Load normalized stack


%for cycle=1:num_cycles
%    for channel=1:num_channels
%        filename = sprintf('%s/reg/cy%02d_ch%02d.tif',embryo_dir,cycle,channel);
%        norm_stack(:,:,:,channel,cycle) = read_3d_tif(filename,embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
%    end
%end

%%Try dogshed
addpath(genpath('\\sodium\broad_thechenlab\Fei\Spot Detection'))
%L = dogshed_fun(norm_stack(:,:,:,1,1));

%save_img(L, sprintf('watershed_channel%02d.tif',1));

%L = dogshed_fun(norm_stack(:,:,:,2,1));

%save_img(L, sprintf('watershed_channel%02d.tif',2));

%L = dogshed_fun(norm_stack(:,:,:,3,1));

%save_img(L, sprintf('watershed_channel%02d.tif',3));

%L = dogshed_fun(norm_stack(:,:,:,4,1));
%save_img(L, sprintf('watershed_channel%02d.tif',4));
%% Deconvolve with Gaussian filter
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- deconvolving ...',time));
deconv_stack = zeros(size(norm_stack));

FILTER_SIZE_XY = 6; % for [fxf] filter 
FILTER_SIZE_Z = 6; 
STDXY1 = 1;%4 
STDXY2 = 4;%5 % btw, spots are black not white if std1 > std2
STDZ1 = 1; 
STDZ2 = 4; 
gauss1 = gauss_filter_fun_xyz(FILTER_SIZE_XY, FILTER_SIZE_Z, STDXY1, STDZ1);
gauss2 = gauss_filter_fun_xyz(FILTER_SIZE_XY, FILTER_SIZE_Z, STDXY2, STDZ2);
dgauss = gauss1 - gauss2;

for cycle=1:num_cycles
    for channel=1:num_channels

            img_filtered = imfilter(norm_stack(:,:,:,channel,cycle), dgauss, 'symmetric', 'conv', 'same'); 

            deconv_stack(:,:,:,channel,cycle) = img_filtered;

    end
end

deconv_stack(deconv_stack<0) = 0;

disp(sprintf('%s: Deconvolved images with Gaussian filter',sec2time(toc)))

disp(sprintf('%s: Saving deconvolved images',sec2time(toc)))

deconv_dir = sprintf('%s/deconv',embryo_dir);
 
for cycle=1:num_cycles
     for channel=1:num_channels
         filename = sprintf('%s/cy%02d_ch%02d.tif',deconv_dir,cycle,channel);
         write_3d_tif(filename, deconv_stack(:,:,:,channel,cycle));
     end
end

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


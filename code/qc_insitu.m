function[] = qc_insitu(home_dir, expt,plate,well,embryo, bounderies)

%% set up environment
tic
cd(home_dir);
addpath(genpath('scripts/'))
 
%%
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s -----...',time));

%%
embryo_dir = sprintf('%s/%s/plate_%d_well_%d/embryo_%d',home_dir, expt,plate,well,embryo);
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

% Load normalized stack


for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/reg/cy%02d_ch%02d.tif',embryo_dir,cycle,channel);
        norm_stack(:,:,:,channel,cycle) = read_3d_tif(filename,embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
    end
end

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

all_stacks_xy_dir = sprintf('%s/figure/all_stacks_xy',embryo_dir);
if ~exist(all_stacks_xy_dir, 'dir') mkdir(all_stacks_xy_dir), end

saveas(fig,sprintf('%s/embryo%d.png',all_stacks_xy_dir,embryo));

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

all_stacks_xz_dir = sprintf('%s/figure/all_stacks_xz',embryo_dir);
if ~exist(all_stacks_xz_dir, 'dir') mkdir(all_stacks_xz_dir), end

saveas(fig,sprintf('%s/embryo%d.png',all_stacks_xz_dir,embryo));

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
cycle_overlap_xy_dir = sprintf('%s/figure/cycle_overlap_xy',embryo_dir);
if ~exist(cycle_overlap_xy_dir, 'dir') mkdir(cycle_overlap_xy_dir), end

saveas(fig,sprintf('%s/embryo%d.png',cycle_overlap_xy_dir));
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


cycle_overlap_xz_dir = sprintf('%s/figure/cycle_overlap_xz',embryo_dir);
if ~exist(cycle_overlap_xz_dir, 'dir') mkdir(cycle_overlap_xz_dir), end

saveas(fig,sprintf('%s/embryo%d.png',cycle_overlap_xz_dir));
disp(sprintf('%s: Saved visualizations',sec2time(toc)))



%% Peak call video

%% Read peaks
deconv_dir = sprintf('%s/deconv',embryo_dir);
for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/cy%02d_ch%02d.tif',deconv_dir,cycle,channel);
        deconv_stack(:,:,:,channel,cycle) = read_3d_tif(filename,embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
    end
end


peaks_table = readtable(sprintf('%s/peaks.txt',embryo_dir));

reg_dapi_stack = uint16(read_3d_tif(sprintf('%s/processed/offset/hyb_dapi.tif',embryo_dir),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6)));
reg_hyb_stack = uint16(read_3d_tif(sprintf('%s/processed/offset/hyb_probe.tif',embryo_dir),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6)));


out_dir = sprintf('%s/figure/peak_calls',embryo_dir);
if ~exist(out_dir, 'dir') mkdir(out_dir), end

video = VideoWriter(sprintf('%s/embryo%d',out_dir, embryo),'MPEG-4');
video.FrameRate = 2.5;
open(video);

z_buffer = 2;
threshold_data_transcriptomeOnly_May2022_auto;
p_w_e = sprintf('plate_%d_well_%d_embryo_%d', plate, well, embryo);
peak_thresh = threshold(p_w_e).cutoff; %determined through histogram

for z=1:size(norm_stack,3)
    
    f = figure('visible','off');
    p = tight_subplot(1,4,[0.001 0.001],[0.001 0.001],[0.001 0.001]);
    
    for channel=1:num_channels
    
        peaks = peaks_table{peaks_table.channel == channel & peaks_table.val > peak_thresh(channel),1:3};
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


cycle = 1;
channel = 1;
for channel = 1:num_channels
    
    CONVOLUTION_KERNEL = ones(5,5,5);
    
    peaks = peaks_table{peaks_table.channel == channel & peaks_table.val > peak_thresh(channel),1:3};
    
    tmp_stack = deconv_stack(:,:,:,channel,cycle);
    tmp_peaks = zeros(size(tmp_stack), 'uint16');
    
    peaks_idx = sub2ind(size(tmp_stack), peaks(:,1), peaks(:,2), peaks(:,3));
    tmp_peaks(peaks_idx) = 1;
    
    tmp_peaks = convn(tmp_peaks, CONVOLUTION_KERNEL, 'same');

    filename = sprintf('%s/cy%02d_ch%02d_peaks.tif',embryo_dir,cycle,channel);
    write_3d_tif(filename, tmp_peaks);
    
end


tic
expt = 'TRANSCRIPTOMICS_May_2022';
home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\ExSeq_Zack';
addpath(home_dir);
plate = 102;
well = 3;
embryo = 27; 

zMin = 0; 
zMax = 150;

for z=1:size(deconv_stack,3)
    
    for channel=1:num_channels
    
        peaks = peaks_table{peaks_table.channel == channel & peaks_table.val > peak_thresh(channel),1:3};
        %peaks = channel_peaks_pos{channel}(channel_peaks_max{channel}>thresh(channel),:);
        z_channel_peaks = peaks(ismember(peaks(:,3),(z+[-z_buffer:+z_buffer])),:);
       
         figure;
         colormap(gray);
         imagesc(deconv_stack(:,:,z,channel,1)); colorbar; caxis([zMin, zMax]); hold on;
         plot(z_channel_peaks(:,2),z_channel_peaks(:,1),'Marker','o','MarkerEdgeColor','r','LineStyle','none');
         title(sprintf('Peak overlaps',size(peaks,1)))
         
         saveas(fig,sprintf('%s/embryo%d_channel_%d_%d.png',embryo_dir,embryo,channel,z));
    end
   

    
end


disp(sprintf('%s: Saved peak call video',sec2time(toc)));
end


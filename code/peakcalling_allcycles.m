function[] = peakcalling_allcycles(home_dir, expt,plate,well,embryo, bounderies)

channel_peaks = {};
num_channels = 4;
num_cycles = 7;

embryo_dir = sprintf('%s/%s/plate_%d_well_%d/embryo_%d',home_dir,expt,plate,well,embryo);
embryo_bounds = dlmread(sprintf('%s/processed/bounds.txt',embryo_dir));
%% %% Read peaks
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- deconv_stacks ...',time));

deconv_stack = zeros(embryo_bounds(5),embryo_bounds(4),embryo_bounds(6),num_channels,num_cycles);

for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/deconv/cy%02d_ch%02d.tif',embryo_dir,cycle,channel);
        deconv_stack(:,:,:,channel,cycle) = read_3d_tif(filename,embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
    end
end


%% Peak calling
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Peak calling ...',time));

channel_peaks = {};
cycle = 4;

for cycle=1:num_cycles
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

writetable(peaks_table,sprintf('%s/peaks_cycle%d.txt',embryo_dir, cycle))

end

% for cycle=1:num_cycles
%     for channel=1:num_channels
%         
%         [Maxima,MaxPos,Minima,MinPos]=MinimaMaxima3D(deconv_stack(:,:,:,channel,cycle),1,0);
% 
%         channel_peaks_pos{cycle}{channel} = MaxPos;
%         channel_peaks_max{cycle}{channel} = Maxima;
% 
%         disp(sprintf('%s: Found %d 3D peaks in cycle %d, channel %d',sec2time(toc),length(channel_peaks_pos{cycle}{channel}),cycle,channel));
%     end          
% end
% 
% time = datestr(datetime(now,'ConvertFrom','datenum'));
% disp(sprintf('%s ----- Writing Peaks ...',time));
% 
% for cycle=1:num_cycles
%     peaks_table = table;
% 
%     peaks_table.x = [channel_peaks_pos{cycle}{1}(:,1); channel_peaks_pos{cycle}{2}(:,1); channel_peaks_pos{cycle}{3}(:,1); channel_peaks_pos{cycle}{4}(:,1)];
%     peaks_table.y = [channel_peaks_pos{cycle}{1}(:,2); channel_peaks_pos{cycle}{2}(:,2); channel_peaks_pos{cycle}{3}(:,2); channel_peaks_pos{cycle}{4}(:,2) ];
%     peaks_table.z = [channel_peaks_pos{cycle}{1}(:,3); channel_peaks_pos{cycle}{2}(:,3); channel_peaks_pos{cycle}{3}(:,3); channel_peaks_pos{cycle}{4}(:,3)];
% 
%     peaks_table.val = [channel_peaks_max{cycle}{1}; channel_peaks_max{cycle}{2}; channel_peaks_max{cycle}{3}; channel_peaks_max{cycle}{4}];
% 
%     peaks_table.channel = [repmat(1,length(channel_peaks_max{cycle}{1}),1); repmat(2,length(channel_peaks_max{cycle}{2}),1); repmat(3,length(channel_peaks_max{cycle}{3}),1); repmat(4,length(channel_peaks_max{cycle}{4}),1);];
%     peaks_table.cycle = [repmat(cycle,length(channel_peaks_max{cycle}{1}),1); repmat(cycle,length(channel_peaks_max{cycle}{2}),1); repmat(cycle,length(channel_peaks_max{cycle}{3}),1); repmat(cycle,length(channel_peaks_max{cycle}{4}),1);];
% 
%     writetable(peaks_table,sprintf('%s/peaks_cycle%d.txt',embryo_dir,cycle))
end
end
function[] = visualize_spot(stack,spot,intensity_mat,match_umi,peak,region_size)

    num_channels = size(stack,4);
    num_cycles = size(stack,5);

    [B L] = bwboundaries(max(spot,[],3),8,'noholes');
    
    figure('visible','off');
    p = tight_subplot(num_channels+1,num_cycles,[0.001 0.001],[0.001 0.001],[0.001 0.001]);

    for cycle=1:num_cycles
    
        cap = prctile(reshape(stack(:,:,:,:,cycle),[],1),99.9);
        
        for channel=1:num_channels+1
            if channel<5
                axes(p((channel-1)*num_cycles+cycle)); imshow(stack(peak(1)+[-region_size(1):region_size(1)],peak(2)+[-region_size(2):region_size(2)],peak(3),channel,cycle),[0 cap]); hold on;
                colormap(jet)
                if str2num(match_umi(cycle)) == channel & intensity_mat(channel,cycle+1) == max(intensity_mat(:,cycle+1))
                    title(sprintf('%.03f',intensity_mat(channel,cycle+1).^2),'Color','r')
                elseif str2num(match_umi(cycle)) == channel
                    title(sprintf('%.03f',intensity_mat(channel,cycle+1).^2),'Color','b')
                else
                    title(sprintf('%.03f',intensity_mat(channel,cycle+1).^2))
                end
                %title(sprintf('%.02e',intensity_mat(channel,cycle+1)))
            else
                axes(p((channel-1)*num_cycles+cycle)); imshow(max(stack(peak(1)+[-region_size(1):region_size(1)],peak(2)+[-region_size(2):region_size(2)],peak(3),:,cycle),[],4),[0 cap]); hold on;
                colormap(jet)
            end
            boundary = B{1};
            plot(boundary(:,2),boundary(:,1),'Color','g', 'LineWidth', 2); hold on;
        end
    end
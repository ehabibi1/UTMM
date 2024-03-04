function[nucleus_volume nucleus_3d_bounds] = nucleus_3d_segmentation(stack) 
   
    xlen = size(stack,1);
    ylen = size(stack,2);
    zlen = size(stack,3);
    
    image = uint16(stack);
    [counts,x] = imhist(uint16(image),max(image(:)));
    [T,EM] = otsuthresh(counts);
    thresh = T*65536;
    
    for z=1:zlen
    
        stack(:,:,z) = stack(:,:,z) > 50;%thresh;
    	stack(:,:,z) = imdilate(stack(:,:,z), strel('disk',4)); 
    	stack(:,:,z) = imfill(stack(:,:,z), 'holes');
    	stack(:,:,z) = imerode(stack(:,:,z), strel('diamond',2));
    	stack(:,:,z) = imerode(stack(:,:,z), strel('diamond',2));
    	stack(:,:,z) = imopen(stack(:,:,z), strel('disk',4));
    
    end
    
    labels = bwlabeln(stack);
    regions = regionprops(labels, 'Area');
    labels_filt = zeros(xlen,ylen,zlen);

    for i = 1:size(regions)
        if regions(i).Area == max([regions.Area])
            labels_filt(labels==i)=i;
        end
    end
        
    nucleus_volume = bwlabeln(labels_filt);
    regions = regionprops(nucleus_volume, 'BoundingBox');
    nucleus_3d_bounds = uint16(regions(1).BoundingBox);
    nucleus_3d_bounds = double([1 1 nucleus_3d_bounds(3) ylen xlen nucleus_3d_bounds(6)]);

    nucleus_volume = imcrop_xyz(nucleus_volume, nucleus_3d_bounds);

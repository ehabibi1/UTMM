
gene = "cdx2";
normalization = "sum";


cellxgene = readtable('cell_counts_table_plate_1_well2_embryo5.txt');

info = imfinfo('plate_1_well_2_embryo_5.tif')

w = [info.Width];
h = [info.Height];

label_matrix = read_3d_tif('plate_1_well_2_embryo_5.tif',h(1),w(1) , length(w));



info = imfinfo('plate_1_well_2_embryo_5_stain_dapi.tif')

w = [info.Width];
h = [info.Height];
dapi_matrix = read_3d_tif('plate_1_well_2_embryo_5_stain_dapi.tif',h(1),w(1) , length(width));


B = imgaussfilt3(dapi_matrix,10);
T = multithresh(B,2);
Bin = imbinarize(B,T(1));

dapiidx = find(Bin>0);



% outmat = zeros(size(dapi_matrix));
% 
% for j = 1:height(cellxgene)
%    idx =  find(cellxgene.Properties.VariableNames == gene)
%    geneval = cellxgene{j,idx};
%    if normalization == "sum"
%    geneval = (geneval*10000 +1)./sum(cellxgene{j,6:end});
%    end
%    labelidx = find(label_matrix == j);
%    
%      temp = zeros(size(label_matrix));
%  
%      temp(labelidx) = 1;
%      temp = imerode(temp, true(3));
%      temp = imerode(temp, true(3));
%      labelidx = find(temp>0);
%      
%    outind = intersect(labelidx,dapiidx);
%    
%    outmat(outind) = geneval;
%    
% end
% figure
% imagesc(squeeze(max(outmat,[],3)))
% outmat = rescale(outmat,0,255) ;
% write_3d_tif([char(gene) '.tif'], outmat);    
% 
% 
% %find boundaries between cells
% 
% downsampleL = imresize3(label_matrix,0.5);
% outbound = zeros(size(downsampleL));
% for j = 1:height(cellxgene)
%      temp = zeros(size(downsampleL));
%      labelidx = find(downsampleL == j);
%      
%      temp(labelidx) = 1;
%      temp = temp -  imerode(temp, true(3));
%      temp = imdilate(temp, true(3));
%      outbound = outbound+temp;
% end
% write_3d_tif('boundaries.tif', outbound);    

for i = 17:width(cellxgene)
gene = string(cellxgene.Properties.VariableNames{i});
figure(100)
cmap = winter(256);
genevalarray = [];
for j = 1:height(cellxgene)
   idx =  find(cellxgene.Properties.VariableNames == gene)
   geneval = cellxgene{j,idx};
   if normalization == "sum"
   geneval = (geneval*10000 +1)./sum(cellxgene{j,6:end});
   end
   genevalarray = [genevalarray geneval];
end  
downsampleL = imresize3(label_matrix,0.5);
for j = 1:height(cellxgene)
   labelidx = find(downsampleL == j);
   
     temp = zeros(size(downsampleL));
 
     temp(labelidx) = 1;
     temp = imerode(temp, true(7));

     [nX,nY,nZ] = size(temp);
    [X,Y,Z] = meshgrid(1:nY,1:nX,1:nZ); 
    [f v] = isosurface(X,Y,Z,temp,0.5);
   figure(100)
    p = patch('Faces',f,'Vertices',v);
    reducepatch(p,40);
    
    mappedcolor = floor((genevalarray(j) - min(genevalarray))./range(genevalarray)*255)+1;
    
    p.FaceColor = cmap(mappedcolor,:);
    
    p.EdgeColor = cmap(mappedcolor,:);
    alpha(0.3)   
end
colormap(winter)
colorbar
view(3)
title(gene)
saveas(gcf, [char(gene) '_volume.png']);
saveas(gcf, [char(gene)  '_volume.fig']);
close(100)

end


%%look at pseudotime

pseudotime = readtable('pseudotime.txt');
pt_reorder = pseudotime([1 9:16 2:8],:);

figure;

ps = pt_reorder.Var4;
downsampleL = imresize3(label_matrix,0.5);
for j = 1:height(cellxgene)
   labelidx = find(downsampleL == j);
   
     temp = zeros(size(downsampleL));
 
     temp(labelidx) = 1;
     temp = imerode(temp, true(7));

     [nX,nY,nZ] = size(temp);
    [X,Y,Z] = meshgrid(1:nY,1:nX,1:nZ); 
    [f v] = isosurface(X,Y,Z,temp,0.5);
   figure(100)
    p = patch('Faces',f,'Vertices',v);
    reducepatch(p,40);
    
    mappedcolor = floor((ps(j) - min(ps))./range(ps)*255)+1;
    
    p.FaceColor = cmap(mappedcolor,:);
    
    p.EdgeColor = cmap(mappedcolor,:);
    alpha(0.3)   
end
colormap(winter)
colorbar
view(3)
title('Pseudotime')
saveas(gcf, ['ps_volume.png']);
saveas(gcf, ['ps_volume.fig']);
close(100)

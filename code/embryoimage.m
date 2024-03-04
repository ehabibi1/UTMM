function imagevoxels = embryoimage(embryotable,geneids)
 
dimensionsx = max(embryotable.x) + 10;
dimensionsy = max(embryotable.y) + 10;
dimensionsz = round(max(embryotable.z)) + 10;
 
imagevoxels = zeros(dimensionsx,dimensionsy,dimensionsz);
embryotable.z = round(embryotable.z);
se4 = strel('sphere', 3)
 
for gene = 1:length(geneids)
    
    geneid = geneids{gene};
    
    geneidx = find(embryotable.gene == string(geneid));
    
    %             temp = zeros(size(imagevoxels));
    %
    %             temp([embryotable.x(j)],[embryotable.y(j)],[embryotable.z(j)]) = 1;
    %             temp = imdilate(temp,se4);
    %             imagevoxels = imagevoxels + temp*gene;
    %
    for idx = 1:length(geneidx)
        j = geneidx(idx);
        imagevoxels ([embryotable.x(j)-2:1:embryotable.x(j)+2],[embryotable.y(j)-2:1:embryotable.y(j)+2],[embryotable.z(j)-2:1:embryotable.z(j)+2]) = gene;
    end
end
 
 
end
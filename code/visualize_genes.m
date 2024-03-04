% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Embryo Ex-Seq (non-ExM) registration of embryo

function[] = process_spots(expt,well,embryo,genes,thresh,img_suffix,out_dir)
%% Manually set parameters and load datasets

%expt = '201211_exseq_trial2';
expt = '210120_without_ab_staining';            % comment out if you want to batch run
well = 1;                                       % comment out if you want to batch run
embryo = 2;                                     % comment out if you want to batch run
out_dir = 'X:\\zchiang/embryo_exseq/figures';   % comment out if you want to batch run

%% Load data

colors = distinguishable_colors(101); colors(4,:) = [];
embryo_dir = sprintf('processed/%s/well%d_embryo%d',expt,well,embryo);

bounds = dlmread(sprintf('%s/bounds.txt',embryo_dir));
dapi_stack = read_3d_tif(sprintf('%s/offset/cy01_dapi.tif',embryo_dir),bounds(5),bounds(4),bounds(6));
gene_xyz_table = readtable(sprintf('figures/%s/gene_xyz_table/gene_xyz_table_well%d_embryo%d.txt',expt,well,embryo));

%% Max projection

genes = ["Gata6","Bhmt","Cdx2"];         % make sure you use double quotes, comment out if you want to batch run
thresh = 4;                                 % 3 = lenient, 4 = stringent, comment out if you want to batch run
img_suffix = "Gata6_Sox2_Cdx2_stringent"   % will be appended to wellX_embryoY_, comment out if you want to batch run

figure; imshow(capImage(max(dapi_stack,[],3),99,'prc'),[]); hold on;

for i=1:size(genes,2);
    sel_gene = genes(i);
    sel_locs = gene_xyz_table(gene_xyz_table.gene == sel_gene & gene_xyz_table.purity >= thresh,:);
    scatter(sel_locs.y,sel_locs.x,10,colors(i,:),'filled');
end

legend(genes)

fig_dir = sprintf('%s/%s/vis_output_max_project',out_dir,expt);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

%saveas(gcf,sprintf('%s/well%d_embryo%d_%s.png',fig_dir,well,embryo,img_suffix));

%% 3D video

genes = ["Fbxo15","Nanog","Cdx2"];         % make sure you use double quotes, comment out if you want to batch run
thresh = 3;                                % 3 = lenient, 4 = stringent, comment out if you want to batch run
img_suffix = "Fbxo15_Nanog_Cdx2_lenient"   % will be appended to wellX_embryoY_, comment out if you want to batch run

fig_dir = sprintf('%s/%s/vis_output_3d_video',out_dir,expt);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

video = VideoWriter(sprintf('%s/well%d_embryo%d_%s',fig_dir,well,embryo,img_suffix),'MPEG-4');
video.FrameRate = 5;
open(video);

figure; 
buf = 1;

for z=1:bounds(6)
    
    imshow(capImage(dapi_stack(:,:,z),99,'prc'),[]); hold on;

    for i=1:size(genes,2);
        sel_gene = genes(i);
        sel_locs = gene_xyz_table(gene_xyz_table.gene == sel_gene & any(gene_xyz_table.z == z+[-buf:buf],2) & gene_xyz_table.purity >= thresh,:);
        scatter(sel_locs.y,sel_locs.x,10,colors(i,:),'filled');
    end

    legend(genes)
    
    writeVideo(video,getframe(gcf))
        
end

close(video)





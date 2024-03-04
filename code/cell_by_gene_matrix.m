% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Embryo Ex-Seq (non-ExM) registration of FoV

function[] = cell_by_gene_matrix(home_dir, expt, species, plate,well, embryo, bounderies, threshold, segmentation, hd_threshold)

%% set up environment
tic
%home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\InSitu_preImpEmbryo';
%data_dir = '\\sodium\broad_thechenlab\ehsan\DataRepository\SplintR';
%SERVER
%home_dir = '/broad/thechenlab/ehsan/analysis/InSitu_preImpEmbryo';
%data_dir = '/broad/thechenlab/ehsan/DataRepository/SplintR';
cd(home_dir);
addpath(genpath('scripts/'))
 
%%
%time = datestr(datetime(now,'ConvertFrom','datenum'));
%disp(sprintf('%s ----- Registering morphology to hyb --> plate %d well %d embryo %d ...',time, plate, well, embryo));

%num_embryos = size(embryos,2);
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Loading  image_bounderies and peak threshold for --> plate %d well %d embryo %d ...',time, plate, well, embryo));


embryo_dir = sprintf('%s/projects/%s/%s/transcriptome/plate_%d_well_%d/embryo_%d',home_dir,species,expt,plate,well,embryo);
embryo_bounds = dlmread(sprintf('%s/processed/bounds.txt',embryo_dir));
%% %% Read peaks
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Reading Peaks ...',time));
peaks_table = readtable(sprintf('%s/peaks.txt',embryo_dir));

num_channels = 4;
num_cycles = 7;

p_w_e = sprintf('plate_%d_well_%d_embryo_%d', plate, well, embryo);
peak_thresh = threshold(p_w_e).cutoff; %determined through histogram

time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- peak_thresh : ',time, peak_thresh));

all_peaks = [];
for channel=1:num_channels
    %all_peaks = [all_peaks; peaks_table{peaks_table.channel == channel & peaks_table.val > thresh(channel),1:3}];
    all_peaks = [all_peaks; peaks_table{peaks_table.channel == channel & peaks_table.val > peak_thresh(channel),1:3}];
end

pad = 2;
mat = zeros(size(all_peaks,1),num_channels,num_cycles);
%purity = zeros(size(all_peaks,1),4);
%consensus = zeros(size(all_peaks,1),4);
purity = zeros(size(all_peaks,1),7);
consensus = zeros(size(all_peaks,1),7);


%load deconv_stack 
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- loading  deconv stack ',time));
deconv_stack = zeros(embryo_bounds(5),embryo_bounds(4),embryo_bounds(6),num_channels,num_cycles);

for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/deconv/cy%02d_ch%02d.tif',embryo_dir,cycle,channel);
        deconv_stack(:,:,:,channel,cycle) = read_3d_tif(filename,embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
    end
end


xleft = max(all_peaks(:,1)-pad,1); xright = min(all_peaks(:,1)+pad,size(deconv_stack,1));
yleft = max(all_peaks(:,2)-pad,1); yright = min(all_peaks(:,2)+pad,size(deconv_stack,2));
zleft = max(all_peaks(:,3)-pad,1); zright = min(all_peaks(:,3)+pad,size(deconv_stack,3));

for i=1:size(all_peaks,1); %disp(i)
    
    peak = all_peaks(i,:);
    peak_mat = squeeze(sum(sum(sum(deconv_stack(xleft(i):xright(i),yleft(i):yright(i),zleft(i):zright(i),:,:),1),2),3));
    
    mat(i,:,:) = peak_mat;
    purity(i,:) = max(peak_mat.^2,[],1)./sum(peak_mat.^2,1);
    [tmp consensus(i,:)] = max(peak_mat,[],1);
    
end

fig = figure;
histogram(mat(:,1,1));
hold on
histogram(mat(:,1,2));
histogram(mat(:,1,3));
histogram(mat(:,1,4));
histogram(mat(:,1,5));
histogram(mat(:,1,6));
histogram(mat(:,1,7));

title('Figure 1')

mkdir(sprintf('%s/figure/misc',embryo_dir))
saveas(fig,sprintf('%s/figure/misc/Figure1_plate_%d_well_%d_embryo_%d',embryo_dir,plate,well,embryo));


time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- quantilenorm',time));
Nmat = zeros(size(mat));
for j = 1:7
Nmattem = quantilenorm(squeeze(mat(:,:,j)));
for i = 1:4
Nmat(:,i,j) = Nmattem(:,i);
end
end

for i=1:size(all_peaks,1); % disp(i)
    
    peak_mat =  squeeze(Nmat(i,:,:));
    purity(i,:) = max(peak_mat.^2,[],1)./sum(peak_mat.^2,1);
    [tmp consensus(i,:)] = max(peak_mat,[],1);
    
end

fig = figure;
histogram(Nmat(:,1,1))
hold on
histogram(Nmat(:,2,1))
histogram(Nmat(:,3,1))

title('Figure 2')

mkdir(sprintf('%s/figure/misc',embryo_dir))
saveas(fig,sprintf('%s/figure/misc/Figure2_plate_%d_well_%d_embryo_%d',embryo_dir,plate,well,embryo));

thresh = 0.7;
purity_score = sum(purity>thresh,2);

fig = figure;
histogram(purity(:,1))
hold on
histogram(purity(:,2))
histogram(purity(:,3))
histogram(purity(:,4))
legend('purity r1', 'purity r2','purity r3','purity r7')

title('Figure 3_1')

mkdir(sprintf('%s/figure/misc',embryo_dir))
saveas(fig,sprintf('%s/figure/misc/Figure3_1_plate_%d_well_%d_embryo_%d',embryo_dir,plate,well,embryo));

fig = figure;
histogram(purity(:,3))
hold on
histogram(purity(:,5))
histogram(purity(:,6))
histogram(purity(:,7))
legend('purity r3', 'purity r5','purity r6','purity r7')

title('Figure 3_2')

mkdir(sprintf('%s/figure/misc',embryo_dir))
saveas(fig,sprintf('%s/figure/misc/Figure3_2_plate_%d_well_%d_embryo_%d',embryo_dir,plate,well,embryo));



fig = figure;
hold on;
scatter(Nmat(:,1,1), Nmat(:,4,1));

scatter(mat(:,1,1), mat(:,4,1));

title('Figure 4')

mkdir(sprintf('%s/figure/misc',embryo_dir))
saveas(fig,sprintf('%s/figure/misc/mat_vs_Nmat_Figure4_4_plate_%d_well_%d_embryo_%d',embryo_dir,plate,well,embryo));

%% Add in barcode info

%  valid_barcodes = [...
%      "0012212","3023023","3011033","2010130",...
%      "0232030","3302023","3122303","2112010","3310103","0102320","0110113","1103130","3231202",...
%      "2021200","1221220","0133213","0031203",...
%      "2033212","2311012","1312121","2002211","0322110",...
%      "3032231","0330132","3101021","1120321","1331011","3203300","0301003","2123022","1022102","2131002",...
%      "1132331","0020020","3321221","0211210","3212321","1202120","2303131","3220033","2332113","1300300",...
%      "2320210","0223310","1233013","3113223","1030202","0313311","2201322"...
%      "3131211","0101032","2112123","0123112","0110022","0032001","0331302","1312213","0131131","2330100","1203103","2232201","0223203","1301212","2120113","3200130",...
%      ]';

valid_barcodes = [...
    "0012212","3023023","3011033","2010130","0232030","3302023","3122303","2112010","3310103","0102320", ...
    "0110133","1103130","3231202","2021200","1221220","0133213","0031203","2033212","2311012","1312121", ...
    "2002211","0322110","3032231","0330132","3101021","1120321","1331011","3203300","0301003","2123022", ...
    "1022102","2131002","1132331","0020020","3321221","0211210","3212321","1202120","2303131","3220033", ...
    "2332113","1300300","2320210","0223310","1233013","3113223","1030202","0313311","2201322","3131211", ...
    "0101032","2112123","0123112","0110022","0032001","0331302","1312213","0131131","2330100","1203103", ...
    "2232201","0223203","1301212","2120113","3200130","3233020","0303200","0302012","0220230","0011330", ...
    "0330010","1100123","3033200","3023301","3131322","3012120","2310201","3022322","3202010","1022021", ...
    "1212310","3031121","3112001","3031313","2302122","0120011","0113030","1323100","0221321","0120332", ...
    "0112313","2001001","0030312","3020100","2123231","3312312","0130323","2033003","3223002","2312232", ...
    "3113012","3001232","3220201","3321300","3232003","1121302","3010010","2012101","3231031","2203030", ...
    "1013131","1220303","2232320","1033120","3212200","2013021","2120030","1331122","3120122","0011002", ...
    "3301030","3221332","2231103","0230113","3003031","3310021","0320202","0221012","0212011","1010122", ...
    "3101200","3102132","2020121","3030022","1020223","3213313","2321102","1302332","1023113","3130013", ...
    "0022300","2110331","3323133","1220131","0023232","2321020","2231010","2103011","2022331","3313202", ...
    ]';
%%
% 
% # ITEM1
% # ITEM2
% 
%2 -- SeqN-1_Ligation_1 = BC Position 1
%3 -- SeqN-1_Ligation_2 = BC Position 2
%4 -- SeqN_Ligation_1 = BC Position 3
%5 -- SeqN_ligation_2 = BC Position 4
%6 -- SeqN-2_Ligation_2 = BC Position 5
%7 -- SeqN-3_Ligation_2 = BC Position 6
%8 -- SeqN-4_Ligation_2 = BC_Position 7
tmp = char(valid_barcodes);
valid_barcodes = string([tmp(:,1) tmp(:,2) tmp(:,3) tmp(:,4) tmp(:,5) tmp(:,6) tmp(:,7)]);

%  genes = [...
%      "Gata3","Krt18","Krt8","Id2",...
%      "Lgals1","Sox17","Cubn","Foxq1","Sparc","Col4a1","Cdx2","Nanog","Bmyc",...
%      "Uap1","Eomes","Gata6","Spic",...
%      "Amot","Eno1","Sox2","Klf4","Gsc",...
%      "Bhmt","Jam2","Tdgf1","Elf3","Fgd1","Satb1","Morc1","Pou5f1","Pak1","Fgf4",...
%      "Tead4","Fgfr2","Otx2","Apoe","Lama1","Lrp2","Fgfr1","Fn1","Slc7a3","Fbxo15",...
%      "Tdh","Dppa2","Dkk1","B3gnt5","B3galnt1","St6galnac4","Acaa2_1",...
%      "Yap1","Taz","Prkcz","Cldn4","Tfap2c","Lats2","Epcam","Cdh1","Gata2","Stat3","Ddah1","S100a11","Rbms1","Dppa1","Csrp1","Atp6v0a4"...
%      ];

genes = [...
   "Gata3","Krt18","Krt8","Id2","Lgals1","Sox17","Cubn","Foxq1","Sparc","Col4a1", ...
    "Cdx2","Nanog","Bmyc","Uap1","Eomes","Gata6","Spic","Amot","Eno1","Sox2", ...
    "Klf4","Gsc","Bhmt","Jam2","Tdgf1","Elf3","Fgd1","Satb1","Morc1","Pou5f1", ...
    "Pak1","Fgf4","Tead4","Fgfr2","Otx2","Apoe","Lama1","Lrp2","Fgfr1","Fn1", ...
    "Slc7a3","Fbxo15","Tdh","Dppa2","Dkk1","B3gnt5","B3galnt1","St6galnac4","Acaa2_1","Yap1", ...
    "Taz","Prkcz","Cldn4","Tfap2c","Lats2","Epcam","Cdh1","Gata2","Stat3","Ddah1", ...
    "S100a11","Rbms1","Dppa1","Csrp1","Atp6v0a4","Bmp4","Cebpa","Dab2","Esrrb","Gata4", ...
    "Klf2","Klf5","Pdgfa","Pdgfra","Pecam1","Tfap2a","Tspan8","Utf1","Gcnt1","Ugcg", ...
    "Carm1","Ctnna1","Iqgap1","Cdc42","Ank","Elf5","Rbpj","Hand1","Sox21","Amotl2", ...
    "Rhoa","Tsc2","Fat1","Rhob","Gja1","Llgl1","Myo9a","Fscn1","Igf2r","Xab2", ...
    "Zfp740","Dnmt3b","Dnmt3a","Dnmt1","Tet1","Tet2","Tet3","Rxra","Zfp820","Rock1", ...
    "Rock2","Sox7","Jup","Tbx3","Dpagt1","Ctnnb1","Rac1","Myc","Cited1","Tgif1", ...
    "Etv5","Hmgb1","Calm2","Lrpap1","Hhex","Egln1","Ugp2","Zfp57","Tcf3","Wwtr1", ...
    "Cldn7","Dppa4","Piezo1","Pmm2","Klf6","Nr5a2","Zfp296","Fabp3","Clic4","Myh9", ...
    "Tmem62","Foxa2","F2r","Gnai2","Mbnl3","Grhl2","Itgb1","Pard6b","Ezr","Nf2", ...
    ];
    
    
    [sort_valid_barcodes sort_order] = sort(valid_barcodes);
sort_genes = genes(sort_order);

%%

sel_consensus = consensus(purity_score>=6,:);
barcodes = strrep(string(num2str(4-sel_consensus)),' ','');



x = ['0' '1' '2' '3'];

possible_barcodes = [];


tic
for i=x
    for j=x
        for k=x
            for l=x
                for m=x
                    for n=x
                        for o=x
                            bc = strcat(i,j,k,l,n,m,o);
                            possible_barcodes = [possible_barcodes; bc];
                            
                        end
                    end
                end
            end
        end
    end                          
end
toc

match_bc_hmD_1 = containers.Map('KeyType','char','ValueType','any');

for bc=1:length(valid_barcodes) %disp(bc)
    value = [];
    for j = 1:length(possible_barcodes)
        y1 = string(regexp(possible_barcodes(j,:), '.', 'match'));
        y2 = regexp(valid_barcodes(bc,:), '.', 'match');
        hd = hamming(y1, y2);
        if  hd <= hd_threshold
            tmp = char(valid_barcodes(bc,:));
            value = [value; possible_barcodes(j,:)];
            match_bc_hmD_1(tmp) = struct('barcodes',[value]);
        end
    end
end

tic
barcodes_with_hd = barcodes;
for ii=1:length(valid_barcodes) %disp(ii)
    [c,ia] = ismember(barcodes_with_hd,match_bc_hmD_1(char(valid_barcodes(ii))).barcodes);
    barcodes_with_hd(c) = valid_barcodes(ii);
end
toc

% % barcodes_with_hd = barcodes;
% % 
% % time = datestr(datetime(now,'ConvertFrom','datenum'));
% % disp(sprintf('%s ----- adding barcodes with hamming distance of 1',time));
% % 
% % for i = 1:length(valid_barcodes) disp(i)
% %     for j = 1:length(barcodes)
% %         y1 = regexp(barcodes(j,:), '.', 'match');
% %         y2 = regexp(valid_barcodes(i,:), '.', 'match');
% %         hd = hamming(y1, y2);
% %         if  hd == 1
% %             %strcat(barcodes(j,:)," -->",valid_barcodes(i,:))
% %             barcodes_with_hd(j,:) = valid_barcodes(i,:);
% %         end
% %     end
% % end
% % 
% % %[C, ia, ic] = unique(barcodes_with_hd);

[C, ia, ic] = unique(barcodes_with_hd);
counts = accumarray(ic,1);

on_target = ismember(C,sort_valid_barcodes);
on_target_prc = sum(counts(on_target))./sum(counts);
on_target_genes = sort_genes(ismember(sort_valid_barcodes,C));

%%

% collect on target and off target indices
on_target_idx = (1:size(C,1))';
off_target_idx = (1:size(C,1))';

on_target_vals = zeros(size(C,1),1);
off_target_vals = zeros(size(C,1),1);



for ii = 1:length(C)    
    if on_target(ii) == 1
        on_target_vals(ii) = counts(ii);
    else
        off_target_vals(ii) = counts(ii);
    end   
        
%         on_target_idx = [on_target_idx; ii];
%     else
%         off_target_idx = [off_target_idx; ii];
%     end
end

% on_target_vals = counts(on_target_idx);
% off_target_vals = counts(off_target_idx);

fig=figure;
hold on;
h1 = bar(on_target_idx, on_target_vals)
set(h1, 'FaceColor', 'r')
h2 = bar(off_target_idx, off_target_vals)
set(h2, 'FaceColor', 'b')

hold off

xlim([1 size(C,1)])

xticks(find(on_target))
xticklabels(strcat(on_target_genes'," - ",C(on_target)," (",string(counts(on_target)),")"))
set(gca,'fontsize',8)

title('Figure 5')

mkdir(sprintf('%s/figure/misc',embryo_dir))
saveas(fig,sprintf('%s/figure/misc/on_target_genes_Figure5_plate_%d_well_%d_embryo_%d_hd_threshold_%d',embryo_dir,plate,well,embryo, hd_threshold));

%%

% figure;
% hold on
% for i = 1:length(C)
%     h=bar(i,counts(i));
%     if on_target(i) == 1
%         set(h,'FaceColor','r');
%     else
%         set(h,'FaceColor','b');
%     end
% end
% hold off
% 
% xlim([1 size(C,1)])
% 
% 
% xticks(find(on_target))
% xticklabels(strcat(on_target_genes'," - ",C(on_target)," (",string(counts(on_target)),")"))
% set(gca,'fontsize',8)

% xticks(find(~on_target))
% xticklabels(strcat(C(~on_target)," (",string(counts(~on_target)),")"))
% set(gca,'fontsize',8)

%xticks(find(~on_target & counts>sum(counts)/256*3))
%xticklabels(strcat("N/A (",C(on_target),")"))

%xtickangle(90)
%title(sprintf("stringent on target %%: %.02f",on_target_prc*100))

%fig = gcf;
%fig.PaperUnits = 'inches';
%fig.PaperPosition = [0 0 10 5];

%fig_dir = sprintf('%s/well_%d/embryo_%d/figure/stringent_on_target',expt,well,embryo);
%if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

%saveas(fig,sprintf('%s/embryo%d.png',fig_dir,embryo));
%saveas(fig,sprintf('%s/well%d_embryo%d.png',fig_dir,well,embryo));

genes_stringent = on_target_genes;
counts_stringent = counts(on_target);


%%
%%Visualize_spots
%all_barcodes = strrep(string(num2str(4-consensus)),' ','');
%for peak=100:1000:size(all_barcodes,1)
%main_title = sprintf('peak %03d, consensus %s, purity %d/4',peak,all_barcodes(peak),purity_score(peak));
%fig = vis_peak(deconv_stack,squeeze(Nmat(peak,:,:)),all_peaks(peak,:),10,2,main_title);
%spot_dir = sprintf('%s/figure/spot_check/plate_%d_well_%d_embryo_%d/',embryo_dir,plate,well,embryo);
%if ~exist(spot_dir, 'dir') mkdir(spot_dir), end
%saveas(fig,sprintf('%s/figure/spot_check/plate_%d_well_%d_embryo_%d/peak%03d.png',embryo_dir,plate,well,embryo,peak));
%end



%%

%sel_consensus = consensus(purity_score>=6,:);
%barcodes = strrep(string(num2str(4-sel_consensus)),' ','');
%[C, ia, ic] = unique(barcodes);
%counts = accumarray(ic,1);

%on_target = ismember(C,sort_valid_barcodes);
%on_target_prc = sum(counts(on_target))./sum(counts);
%on_target_genes = sort_genes(ismember(sort_valid_barcodes,C));

%figure;
%hold on
%for i = 1:length(C)
%    h=bar(i,counts(i));
%    if on_target(i) == 1
%        set(h,'FaceColor','r');
%    else
%        set(h,'FaceColor','b');
%    end
%end
%hold off

%xlim([1 size(C,1)])

%xticks([find(counts>sum(counts)/256*3)])
%xticklabels(C(counts>sum(counts)/256*3))

%xticks(find(on_target))
%xticklabels(strcat(on_target_genes'," - ",C(on_target)," (",string(counts(on_target)),")"))
%set(gca,'fontsize',8)

% xticks(1:length(C))
% xticklabels(strcat(C," (",string(counts),")"))
% set(gca,'fontsize',8)

%xticks(find(~on_target & counts>sum(counts)/256*3))
%xticklabels(strcat("N/A (",C(on_target),")"))

%xtickangle(90)

%xticks(find(~on_target & counts>sum(counts)/256*3))
%xticklabels(strcat("N/A (",C(on_target),")"))


%title(sprintf("lenient on target %%: %.02f",on_target_prc*100))

%fig = gcf;
%fig.PaperUnits = 'inches';
%fig.PaperPosition = [0 0 10 5];

%fig_dir = sprintf('%s/plate_%d_well_%d/embryo_%d/figure/lenient_on_target',expt,plate,well,embryo);
%if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

%saveas(fig,sprintf('%s/embryo%d.png',fig_dir,embryo));
%saveas(fig,sprintf('%s/plate_%d_well%d_embryo%d.png',fig_dir,plate,well,embryo));

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

fig_dir = sprintf('%s/figure/gene_counts_table',embryo_dir);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

writetable(gene_counts,sprintf('%s/gene_counts_plate_%d_well%d_embryo_%d_hd_threshold_%d.txt',fig_dir,plate,well,embryo,hd_threshold))



%% 
%sel_consensus = consensus(purity_score>=6,:);
%barcodes = strrep(string(num2str(4-sel_consensus)),' ','');

sel_purity = purity_score(purity_score>=6);
%on_target_purity = sel_purity(ismember(barcodes,sort_valid_barcodes));
on_target_purity = sel_purity(ismember(barcodes_with_hd,sort_valid_barcodes));


sel_peaks = all_peaks(purity_score>=6,:);
%on_target_peaks = sel_peaks(ismember(barcodes,sort_valid_barcodes),:);
on_target_peaks = sel_peaks(ismember(barcodes_with_hd,sort_valid_barcodes),:);

barcode_dict = strings(4000,1);
barcode_dict(str2double(sort_valid_barcodes)) = sort_genes';
%on_target_genes = barcode_dict(str2double(barcodes(ismember(barcodes,sort_valid_barcodes))));
on_target_genes = barcode_dict(str2double(barcodes_with_hd(ismember(barcodes_with_hd,sort_valid_barcodes))));

gene_xyz_table = table;
gene_xyz_table.gene = on_target_genes;
gene_xyz_table.x = on_target_peaks(:,1);
gene_xyz_table.y = on_target_peaks(:,2);
gene_xyz_table.z = on_target_peaks(:,3);
gene_xyz_table.purity = on_target_purity;

fig_dir = sprintf('%s/figure/gene_xyz_table',embryo_dir);
if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

writetable(gene_xyz_table,sprintf('%s/gene_xyz_table_plate_%d_well%d_embryo_%d_hd_threshold_%d.txt',fig_dir,plate,well,embryo,hd_threshold))

% filter reads by disttance 
geneids = unique(gene_xyz_table.gene);

gene_xyz_distancefiltered_table = table();

for gene = 1:length(geneids)
            geneid = geneids{gene};
            geneidx = find(gene_xyz_table.gene == string(geneid));

            x = gene_xyz_table.x(geneidx);
            y = gene_xyz_table.y(geneidx);
            z = gene_xyz_table.z(geneidx);

            coords = [x y z];

            dists = pdist(coords,'cityblock');

            squaredists=squareform(dists);
            squaredists = triu(squaredists,1);
            closeinds = find(squaredists<5 & squaredists>0);
            [row, col] = ind2sub(size(squaredists),closeinds);


            coords(row,:) = [];

            geneidx(row) =[];
            new_data_gene = gene_xyz_table(geneidx,:);
            gene_xyz_distancefiltered_table = [gene_xyz_distancefiltered_table; new_data_gene];
end

writetable(gene_xyz_distancefiltered_table,sprintf('%s/gene_xyz_table_plate_%d_well%d_embryo_%d_hd_threshold_%d_distancefiltered.txt',fig_dir,plate,well,embryo,hd_threshold))

time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Loading  segmentation mask ...', time));


if segmentation == 1 
    %% Load segmentation
    ReReg_dir = sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/ReReg',home_dir,species,expt, plate, well, embryo);
   
    if exist(ReReg_dir, 'dir')
        time = datestr(datetime(now,'ConvertFrom','datenum'));
        disp(sprintf('%s ----- Loading  segmentation mask from ReReg dir ...', time));
        segmentation_dir = ReReg_dir;
    else
        time = datestr(datetime(now,'ConvertFrom','datenum'));
        disp(sprintf('%s ----- Loading  segmentation mask  ...', time));
        segmentation_dir = sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/',home_dir,species,expt, plate, well, embryo);;
    end 
    
    seg_stack = read_3d_tif(sprintf('%s/p%dw%de%d-membSeg-relabeled.tif',segmentation_dir,plate,well,embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
    reg_stain_Ecad = read_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Ecad.tif',segmentation_dir, plate,well,embryo),embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
    figure; imshowpair(capImage(max(seg_stack,[],3),99,'prc'),capImage(max(reg_stain_Ecad,[],3),99,'prc'))
    %% Get assignments
    labels = unique(seg_stack(:));
    
    A = labels(labels~=0);
    num_labels = length(A);
    A_tmp = (1:length(A))';
    for z = 1:size(seg_stack,3)
        seg_stack_corrected(:,:,z) = changem(seg_stack(:,:,z),A_tmp',A');
    end 
     
    %sel = gene_xyz_table.purity >= 6; % 3 pr 4
    sel = gene_xyz_distancefiltered_table.purity >= 6; % 3 pr 4
    %gene_idx = sub2ind(size(seg_stack_corrected),gene_xyz_table.x(sel),gene_xyz_table.y(sel),gene_xyz_table.z(sel));
    gene_idx = sub2ind(size(seg_stack_corrected),gene_xyz_distancefiltered_table.x(sel),gene_xyz_distancefiltered_table.y(sel),gene_xyz_distancefiltered_table.z(sel));

    cell_assignment = seg_stack_corrected(gene_idx);

    %% transcripts per cell

    assigned = cell_assignment(cell_assignment ~= 0);
    [C ia ic] = unique(assigned);
    counts = accumarray(ic,1);
    
    num_cells = max(seg_stack_corrected(:));
    genes_per_cell = zeros(num_cells,1);
   
    genes_per_cell(C) = counts;
    figure; bar(genes_per_cell); hold on;
    

     
    colors = distinguishable_colors(101); colors(4,:) = [];
    fig = figure; hold on;
    for i = 1:length(genes_per_cell)
        h=bar(i,genes_per_cell(i));
        set(h,'FaceColor',colors(i,:));
    end

    xlim([0 num_cells])
    xlabel('cell index'); ylabel('# transcripts')
    title(sprintf('%.2f%% of on-target transcripts assigned',sum(cell_assignment>0)./size(cell_assignment,1)*100))

    mkdir(sprintf('%s/figure/misc',embryo_dir))
    saveas(fig,sprintf('%s/figure/misc/transcripts_per_cell_plate_%d_well_%d_embryo_%d_hd_threshold_%d_distancefiltered.png',embryo_dir,plate,well,embryo,hd_threshold));

%% create matrix
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- generating cell by gene matrix ...'));


    sel = gene_xyz_distancefiltered_table.purity >= 6; % 3 pr 4
    num_genes = size(valid_barcodes,1);
    gene_xyz_distancefiltered_table.cell_assignment = zeros(size(gene_xyz_distancefiltered_table,1),1);
    gene_xyz_distancefiltered_table.cell_assignment(sel) = cell_assignment;

    cell_counts_mat = zeros(length(C),num_genes);

    for i=1:num_genes

        sel_cells = gene_xyz_distancefiltered_table.cell_assignment(gene_xyz_distancefiltered_table.gene == sort_genes(:,i));
        [C ia ic] = unique(sel_cells(sel_cells>0));
        counts = accumarray(ic,1);
        cell_counts_mat(C,i) = counts;
    end

    cell_counts_table = array2table(cell_counts_mat);
    cell_counts_table.Properties.VariableNames = cellstr(sort_genes');
    cell_counts_table.cell_index = (1:size(cell_counts_table,1))';
   

    cell_info = regionprops(seg_stack_corrected,'Centroid','Area');
    centroids = reshape([cell_info.Centroid]',3,[])';
    cell_counts_table.x = centroids(:,2)*.17;
    cell_counts_table.y = centroids(:,1)*.17;
    cell_counts_table.z = centroids(:,3)*.4;
    cell_counts_table.area = [cell_info.Area]';
    
    cell_counts_table = [cell_counts_table(:,end-4:end) cell_counts_table(:,1:end-5)]; 
    %cell_counts_table = cell_counts_table(cell_counts_table.area > 100,:);
     cell_counts_table.seg_label = string(A);

    fig_dir = sprintf('%s/figure/cell_counts_table',embryo_dir);
    if ~exist(fig_dir, 'dir') mkdir(fig_dir), end

    writetable(cell_counts_table,sprintf('%s/cell_counts_table_plate_%d_well%d_embryo_%d_hd_threshold_%d_distancefiltered.txt',fig_dir,plate,well,embryo, hd_threshold))

    disp('DONE!...');
    
end
end
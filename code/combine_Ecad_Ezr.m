function[] = merge_ecad_ezr(plate, well, embryo,mergeAb)
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- merging Ecad & Ezr plate %d well %d embryo %d ...',time, plate, well, embryo));

home_dir = 'W:\ehsan\analysis\ExSeq_Zack\';
cd(home_dir);
addpath(genpath('scripts/'))
addpath('W:\ehsan\bucket\ehsan')

embryo_dir = sprintf('//sodium/broad_thechenlab/ehsan/analysis/ExSeq_Zack/TRANSCRIPTOMICS_May_2022/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/',plate,well,embryo);
addpath(genpath('\\sodium\broad_thechenlab\ehsan\analysis\ExSeq_Zack\'));

time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Reading Ecad plate %d well %d embryo %d ...',time, plate, well, embryo));

ecad=FastTiff(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Ecad.tif',embryo_dir,plate, well, embryo));

disp(sprintf('%s ----- Reading Ezr plate %d well %d embryo %d ...',time, plate, well, embryo));

ezr=FastTiff(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Ezr.tif',embryo_dir,plate, well, embryo));

%ecad = FastTiff('W:\ehsan\analysis\ExSeq_Zack\TRANSCRIPTOMICS_May_2022\segmentation\reg_to_hyb_dapi\plate_102\well_3\embryo_3\plate_102_well_3_embryo_3_stain_Ecad.tif');
%ezr = FastTiff('W:\ehsan\analysis\ExSeq_Zack\TRANSCRIPTOMICS_May_2022\segmentation\reg_to_hyb_dapi\plate_102\well_3\embryo_3\plate_102_well_3_embryo_3_stain_Ezr.tif');
median_ecad = median(ecad(:));
median_ezr = median(ezr(:));
ecadn = ecad - median_ecad;
ezrn = ezr - median_ezr;


ecadn(ecadn<0) = 0;
ezrn(ezrn<0) = 0;

ecadn_ezrn_tmp = (ecadn/prctile(ecadn(:),99) + ezrn/prctile(ezrn(:),99));%prctile(ezrn(:),95)

%subplot(2, 2, 1)
%histogram(log2(peaks_table.val(peaks_table.channel==1)), 'BinLimits', [-5, 10]);
%prctile(ecadn(:),95)

ecadn_ezrn_tmp(ecadn_ezrn_tmp>1) = 1;
ecadn_ezrn = 2^16 * ecadn_ezrn_tmp;

%ecadn_ezrn = 2^16 * ecadn_ezrn_tmp/2;

write_3d_tif(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_EcadEzr_merged.tif',embryo_dir,plate, well, embryo),ecadn_ezrn);
%write_3d_tif(sprintf('%s/p%dw%de%d-memb.tif',embryo_dir,plate, well, embryo),ecadn_ezrn);
if mergeAb == 1
    merged_file = sprintf('%s/plate_%d_well_%d_embryo_%d_stain_EcadEzr_merged.tif',embryo_dir,plate, well, embryo);
    merged_file_copy = sprintf('%s/p%dw%de%d-memb.tif',embryo_dir,plate, well, embryo);
    copyfile(merged_file, merged_file_copy);
    
elseif mergeAb == 0
    merged_file = sprintf('%s/plate_%d_well_%d_embryo_%d_stain_Ecad.tif',embryo_dir,plate, well, embryo);
    merged_file_copy = sprintf('%s/p%dw%de%d-memb.tif',embryo_dir,plate, well, embryo); 
    copyfile(merged_file, merged_file_copy);
end
    
    
    

end



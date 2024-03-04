%%
tic
expt = 'TRANSCRIPTOMICS_May_2022';
%addpath(genpath('\\sodium\brod_thechenlab\ehsan/analysis/ExSeq_Zack/TRANSCRIPTOMICS_May_2022/'))
%server
%addpath(genpath('/broad/thechenlab/ehsan/analysis/ExSeq_Zack/TRANSCRIPTOMICS_May_2022/'))
home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\ExSeq_Zack';
data_dir = '\\sodium\broad_thechenlab\ehsan\DataRepository\SplintR';
%SERVER
%home_dir = '/broad/thechenlab/ehsan/analysis/ExSeq_Zack';
%data_dir = '/broad/thechenlab/ehsan/DataRepository/SplintR';
cd(home_dir);
%addpath(genpath(home_dir));
addpath(home_dir);

addpath(genpath('scripts/'))
image_bounderies_transcriptomeOnly_May2022;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Processing data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plate = 1;
% well = 2;
% embryos = [1 2 3 4 5 6 7]; 
% embryos_in_fovs = ["1" "2" "3+4+5" "6" "7"];
% 
% for i=1:length(embryos_in_fovs)
%     embryos_in_fov = embryos_in_fovs(i);
%     numberOfembryo = split(embryos_in_fov,"+");
%     if length(numberOfembryo) == 1
%         embryo = str2num(str2mat(numberOfembryo));
%         process_fov_v1(expt,plate,well,embryos_in_fov,embryo, bounderies);
%         process_reg_hyb_to_morpho(expt,plate,well,embryos_in_fov,embryo, bounderies);
%         norm_deconv_peakcalling(expt,plate,well,embryo, bounderies);
%     else
%         for embryo=str2num(str2mat(numberOfembryo))
%             process_fov_v1(expt,plate,well,embryos_in_fov,embryo, bounderies);
%             process_reg_hyb_to_morpho(expt,plate,well,embryos_in_fov,embryo, bounderies);
%             norm_deconv_peakcalling(expt,plate,well,embryo, bounderies);
%         end
%     end
% end


plate = 102;
well = 3;
embryos = [3:6 10:35]; 
for embryo=embryos
    plot_peak_threshold(expt, plate, well, embryo)
end

plate = 102;
well = 4;
embryos = [1:22]; 
for embryo=embryos
    plot_peak_threshold(expt, plate, well, embryo)
end




tic
expt = '20220615';
home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\InSitu_preImpEmbryo\';
cd(home_dir);
addpath(home_dir)
addpath(genpath('scripts/'))
species = 'mouse';
plate = 201;
well = 1;
embryos=[9,11:13];
image_bounderies_mechanotranscriptomics_2022;
ReReg = 0;
for embryo=embryos
    embryo
    %chech_seg(home_dir, expt, species,plate,well,embryo, bounderies, ReReg)
    reRegister_process_reg_hyb_to_morpho_server_finereg(home_dir, expt,species, plate,well,embryo, bounderies)
end

plate = 102;
well = 4;
embryos=[8:12,15,18,20];
for embryo=embryos
    embryo
    reRegister_process_reg_hyb_to_morpho_server_finereg(home_dir, expt,species, plate,well,embryo, bounderies)
end

plate = 102;
well = 3;
embryos=[7,9,17,27,33];
for embryo=embryos
    embryo
    reRegister_process_reg_hyb_to_morpho_server_finereg(home_dir, expt,species, plate,well,embryo, bounderies)
end


raw_data_dir_name = 'MECHANOTRANSCRIPTOMICS_July_2022';
expt = '20220726';
home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\InSitu_preImpEmbryo';
data_dir = '\\sodium\broad_thechenlab\ehsan\DataRepository\SplintR';
species = 'mouse';
cd(home_dir);
addpath(genpath('scripts/'));
addpath(genpath('/broad/thechenlab/Fei/Spot Detection'));
image_bounderies_mechanotranscriptomics_2022;
plate = 202;
well = 1;
embryos_in_fov = "11";
embryos=[11];
for embryo=embryos
    embryo
    process_reg_hyb_to_morpho_server_finereg(data_dir, home_dir, raw_data_dir_name, expt, species,plate,well,embryos_in_fov,embryo, bounderies, 1, 3, 4);
end

embryos_in_fov = "12+13";
embryos=[12:13];
for embryo=embryos
    embryo
    process_reg_hyb_to_morpho_server_finereg(data_dir, home_dir, raw_data_dir_name, expt, species,plate,well,embryos_in_fov,embryo, bounderies, 1, 3, 4);
end

embryos_in_fov = "14";
embryos=[14];
for embryo=embryos
    embryo
    process_reg_hyb_to_morpho_server_finereg(data_dir, home_dir, raw_data_dir_name, expt, species,plate,well,embryos_in_fov,embryo, bounderies, 1, 3, 4);
end






tic
expt = '20220726';
home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\InSitu_preImpEmbryo\';
cd(home_dir);
addpath(home_dir)
addpath(genpath('scripts/'))
species = 'mouse';
plate = 203;
well = 1;

%well 4 embryo [12:15, 22] needs to be re-registered
embryos = [3,4,5,6,7,8:10,12:18,20:22,23,24,25]
%,11,19
for embryo=embryos
    
    embryo
    segmentation_dir = sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/',home_dir,species,expt, plate, well, embryo);
    
    temp_dir = sprintf('%s/temp', segmentation_dir);
    mkdir(temp_dir);
    %dapi=FastTiff(sprintf('%s/plate_%d_well_%d_embryo_%d_stain_dapi.tif',embryo_dir,plate, well, embryo));
    dapi_file=sprintf('%s/plate_%d_well_%d_embryo_%d_stain_dapi.tif',segmentation_dir,plate, well, embryo);
    %write_3d_tif(sprintf('%s/p%dw%de%d_1.tif',embryo_dir,plate, well, embryo),dapi);
    dapi_resave=sprintf('%s/p%dw%de%d_1.tif',segmentation_dir,plate, well, embryo);
    copyfile(dapi_file, dapi_resave);
    preprocess_dapi(home_dir, species, expt, plate, well, embryo);
    param_path = 'W:\ehsan\tools\StarryNite-master\StarryNite-master\example_parameter_files\newmatlab\';
    SD_red_40x_2 = matlab.desktop.editor.openDocument(fullfile(param_path,'SD_red_40x_2.txt'));
    %waitfor(SD_red_40x_2,'Closed',1);  % this line waits for the editor property 'Opened' to change from a 1 to a 0 when the file is closed
    
    %combine_Ecad_Ezr(plate, well, embryo, 1);%0 Ecad ---- run on server
    preprocess_memb(home_dir, species, expt, plate, well, embryo);
    
    mydlg1 = warndlg('After modifying the SD_red_40x_2 file, press OK to continue');
    waitfor(mydlg1);
    %disp('This prints after you close the warning dialog.');
    %%disp(' ----- cd  StarryNite -----');
    poolobj = gcp('nocreate');
    delete(poolobj);
    cd 'W:\ehsan\tools\StarryNite-master\StarryNite-master';
    addpath(genpath(pwd));
    %disp(' ----- cd  distribution_lineaging -----');
    cd 'distribution_lineaging';
    lineage_launcher_v2;
    
    poolobj = gcp('nocreate');
    delete(poolobj);
    
    mydlg2 = warndlg('After StarryNite run finished, press OK to continue');
    waitfor(mydlg2);
    
    modify_StarryNite_xml(home_dir, species, expt, plate, well, embryo);
    cd 'W:\ehsan\tools\AceTree-master'; 
    system( 'powershell -command "java -jar AceTree.jar"' );
    
    mydlg2 = warndlg('After running AceTree, press OK to continue');
    waitfor(mydlg2);
    
    cd(temp_dir);
    to_nuclei_csv(home_dir, species, expt, plate, well, embryo);
    poolobj = gcp('nocreate');
    delete(poolobj);
    
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
expt = '20220726';
home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\InSitu_preImpEmbryo\';
cd(home_dir);
addpath(home_dir)
addpath(genpath('scripts/'))
species = 'mouse';
cd 'W:\ehsan\tools\3DMMS\3DMMS_public-master';
addpath(genpath(pwd))
plate = 203;
well = 1;
embryos = [14]%4:10,12:18,20:25
for embryo=embryos
    poolobj = gcp('nocreate');
    delete(poolobj);
    
    embryo
    
    I3DMMS_payman(home_dir, expt, species, plate, well, embryo)
    %memb_postSeg(home_dir, expt, plate, well, embryo)
    
end

home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\InSitu_preImpEmbryo\';
species = 'mouse';
cd(home_dir);
addpath(genpath('scripts/'))
expt = '20220726';
plate = 203;
well = 1;
embryos = [6]; %4:10,12:16
for embryo=embryos
    
    embryo_dir = sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/ReReg', home_dir,species, expt, plate, well, embryo);
    
    if exist(embryo_dir, 'dir')
        
        disp(sprintf('embryo %d was reregistered ...',embryo));
        membSeg = FastTiff(sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/ReReg/p%dw%de%d-membSeg.tif', home_dir,species, expt, plate, well, embryo,plate, well, embryo));
        [membSeg_corLabel table] = relabel_segMask(membSeg);
        write_3d_tif(sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/ReReg/p%dw%de%d-membSeg-relabeled.tif', home_dir,species, expt, plate, well, embryo,plate, well, embryo),membSeg_corLabel);
        filename = sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/ReReg/p%dw%de%d_membSeg_OldNewLabels.txt',home_dir,species, expt, plate, well, embryo,plate, well, embryo);
        writetable(table, filename);
        
    else
        disp(sprintf('embryo %d was NOT re-registered ...',embryo));
        membSeg = FastTiff(sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/p%dw%de%d-membSeg.tif', home_dir,species, expt, plate, well, embryo,plate, well, embryo));
        [membSeg_corLabel table] = relabel_segMask(membSeg);
        write_3d_tif(sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/p%dw%de%d-membSeg-relabeled.tif', home_dir,species, expt, plate, well, embryo,plate, well, embryo),membSeg_corLabel);
        filename = sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/p%dw%de%d_membSeg_OldNewLabels.txt',home_dir,species, expt, plate, well, embryo,plate, well, embryo);
        writetable(table, filename);
    end 
    
    

end






tic
plate = 102;
well = 3;
embryos = [3:6 10:35]; 
for embryo=embryos
    cell_by_gene_matrix(expt, plate, well, embryo, bounderies, threshold, 0)
end
toc





%javaaddpath('/opt/MATLAB/R2019b/java/jar/ij-1.53c.jar');
%javaaddpath('/opt/MATLAB/R2019b/java/jar/mij.jar');
javaaddpath 'C:\Program Files\MATLAB\R2017a\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2017a\java\ij.jar'
addpath('C:\Users\ehabibi\Downloads\fiji-win64\Fiji.app\scripts');

home_dir = 'W:\ehsan\analysis\InSitu_preImpEmbryo\projects';
cd('W:\ehsan\analysis\InSitu_preImpEmbryo\mechanics_transcriptome_alignment_pipeline\')
addpath(genpath(pwd))
species = 'mouse';
doDisplacement = false;
dropUncertainMeasures = true;


expt = '20220615';
exptxlsfile = 'June 2022 Embryo Mechanics Metadata.xlsx';
exptxlsfile
rounds = ["1" "5"]; %"2" "2.1" "3" "4"
sampleIds = [1:3,5:8,10,11,13,14,19:26];%2 and 4 have issues 1:3,5:8,10,11,13,14,19:26
for sampleId=sampleIds
    for round=rounds
    sampleId
    round
    viz = bafun({VizOpts.Visualize,VizOpts.Assignment,VizOpts.PatchBounds, VizOpts.CHullCenter, VizOpts.Registered},@bitor);
    alignment_pipeline_mechanotranscriptome(home_dir, species, expt, exptxlsfile,sampleId, viz, doDisplacement, dropUncertainMeasures, round)
    end
end


expt = '20220726';
exptxlsfile = 'July 2022 Embryo Mechanics Metadata_v2.xlsx';
exptxlsfile
rounds = ["1" "4"];% "2" "3" "3.1"
sampleIds = [1,3,13,15,17,21,23,26,27,29]; 
for sampleId=sampleIds
    for round=rounds
    sampleId
    viz = bafun({VizOpts.Visualize,VizOpts.Assignment,VizOpts.PatchBounds, VizOpts.CHullCenter, VizOpts.Registered},@bitor);
    alignment_pipeline_mechanotranscriptome(home_dir, species, expt, exptxlsfile,sampleId, viz, doDisplacement, dropUncertainMeasures, round)
    end
end
dropUncertainMeasures = false;

round = 1;
sampleIds = [18]; %2 and 4 have issues 1:3,5:8,10,11,13,14,19:26
for sampleId=sampleIds
    sampleId
    viz = bafun({VizOpts.Visualize,VizOpts.Assignment,VizOpts.PatchBounds, VizOpts.CHullCenter, VizOpts.Registered},@bitor);
    alignment_pipeline_mechanotranscriptome(home_dir, species, expt, exptxlsfile,sampleId, viz, doDisplacement, dropUncertainMeasures, round)
end


expt = '20220726';
exptxlsfile = 'July 2022 Embryo Mechanics Metadata.xlsx';
exptxlsfile
doDisplacement = false;
round = 3;
sampleIds = [1,3,13,15,17,21,23,26]; 

dropUncertainMeasures = false;

for sampleId=sampleIds
    sampleId
    viz = bafun({VizOpts.Visualize,VizOpts.Assignment,VizOpts.PatchBounds, VizOpts.CHullCenter, VizOpts.Registered},@bitor);
    alignment_pipeline_mechanotranscriptome(home_dir, species, expt, exptxlsfile,sampleId, viz, doDisplacement, dropUncertainMeasures, round)
end

dropUncertainMeasures = true;

for sampleId=sampleIds
    sampleId
    viz = bafun({VizOpts.Visualize,VizOpts.Assignment,VizOpts.PatchBounds, VizOpts.CHullCenter, VizOpts.Registered},@bitor);
    alignment_pipeline_mechanotranscriptome(home_dir, species, expt, exptxlsfile,sampleId, viz, doDisplacement, dropUncertainMeasures, round)
end

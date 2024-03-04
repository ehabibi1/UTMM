function[] = preprocess_dapi(home_dir, species, expt, plate,well,embryo)
%MATLAB 2017a : https://www.mathworks.com/matlabcentral/fileexchange/47545-mij-running-imagej-and-fiji-within-matlab
javaaddpath 'C:\Program Files\MATLAB\R2017a\java\mij.jar'
javaaddpath 'C:\Program Files\MATLAB\R2017a\java\ij.jar'

addpath('C:\Users\ehabibi\Downloads\fiji-win64\Fiji.app\scripts');

%ImageJ 2021
MIJ.start() %2017
IJ=ij.IJ(); 

%macro_path = sprintf('%s/%s/segmentation/reg_to_hyb_dapi',home_dir,expt); 
macro_path = home_dir;

segmentation_dir = sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/',home_dir,species,expt, plate, well, embryo);

dapi_path = sprintf('%s/p%dw%de%d_1.tif',segmentation_dir,plate,well,embryo);
 
args=...
        strcat('image_path=',dapi_path);
    
IJ.runMacroFile(java.lang.String(fullfile(macro_path,'dapi_preprocess.ijm'))...
    ,java.lang.String(args));

IJ.run('Close All');

end

%You need to add paths to ij.jar file and such before running the code, e.g.:
javaaddpath('C:\Program Files\MATLAB\R2017a\java\mij.jar');
addpath('\\sodium\broad_thechenlab\ehsan\tools\Fiji.app/scripts/')
addpath('\\sodium\broad_thechenlab\ehsan\tools\SlicerMatlabBridge\');
addpath('\\sodium\broad_thechenlab\ehsan\tools\SlicerMatlabBridge\MatlabCommander\commandserver\');
javaaddpath('C:\Users\ehabibi\Downloads\fiji-win64\Fiji.app\jars\ij-1.53c.jar');
%javaaddpath('/opt/MATLAB/R2019b/java/jar/ij-1.53c.jar');
addpath('\\sodium\broad_thechenlab\ehsan\analyis\payman');

viz = bafun({VizOpts.Visualize,VizOpts.Assignment,VizOpts.PatchBounds, VizOpts.CHullCenter, VizOpts.Registered},@bitor)
alignment_pipeline('mechano_transcriptome_project','joint_measurement_may_2021','Embryo Samples Mechanics Metadata.xlsx', 14, viz, false, false)



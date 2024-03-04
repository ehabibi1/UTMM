addpath('W:\ehsan\tools\Fiji.app\scripts')
javaaddpath('W:\ehsan\tools\Fiji.app\jars\ij-1.53c.jar')
addpath('W:\ehsan\analysis\SlicerMatlabBridge');
addpath('W:\ehsan\analysis\mechano_transcriptome_project');
addpath('W:\ehsan\analysis\');
addpath('W:\ehsan\bucket\ehsan\SlicerMatlabBridge\MatlabCommander\commandserver')
viz = bafun({VizOpts.Visualize,VizOpts.Assignment,VizOpts.PatchBounds, VizOpts.CHullCenter, VizOpts.Registered},@bitor)
alignment_pipeline('mechano_transcriptome_project','joint_measurement_may_2021','Embryo Samples Mechanics Metadata.xlsx', 14, viz, false, false)
% Author: Zachary Chiang, Buenrostro Lab, Harvard University
% Embryo Ex-Seq (non-ExM) registration of FoV

function[] = deconv(expt,plate,well,embryo, bounderies)

%% set up environment
tic
home_dir = '\\sodium\broad_thechenlab\ehsan\analysis\ExSeq_Zack';
data_dir = '\\sodium\broad_thechenlab\ehsan\DataRepository\SplintR';
cd(home_dir);
addpath(genpath('scripts/'))
 
%%
embryo_dir = sprintf('%s/plate_%d_well_%d/embryo_%d',expt,plate,well,embryo);
p_w_e = sprintf('plate_%d_well_%d_embryo_%d', plate, well, embryo);
embryo_bounds = bounderies(p_w_e).hyb;

reg_dir = sprintf('%s/reg',embryo_dir);
if ~exist(reg_dir, 'dir') mkdir(reg_dir), end

%% Load stacks

num_channels = 4;
num_cycles = 7;

%embryo_bounds = dlmread(sprintf('%s/processed/bounds.txt',embryo_dir));
norm_stack = zeros(embryo_bounds(5),embryo_bounds(4),embryo_bounds(6),num_channels,num_cycles);

% Load normalized stack

time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- Load normalized stack ...',time));
for cycle=1:num_cycles
    for channel=1:num_channels
        filename = sprintf('%s/reg/cy%02d_ch%02d.tif',embryo_dir,cycle,channel);
        norm_stack(:,:,:,channel,cycle) = read_3d_tif(filename,embryo_bounds(5),embryo_bounds(4),embryo_bounds(6));
    end
end

%%Try dogshed
addpath(genpath('\\sodium\broad_thechenlab\Fei\Spot Detection'))
%L = dogshed_fun(norm_stack(:,:,:,1,1));

%save_img(L, sprintf('watershed_channel%02d.tif',1));

%L = dogshed_fun(norm_stack(:,:,:,2,1));

%save_img(L, sprintf('watershed_channel%02d.tif',2));

%L = dogshed_fun(norm_stack(:,:,:,3,1));

%save_img(L, sprintf('watershed_channel%02d.tif',3));

%L = dogshed_fun(norm_stack(:,:,:,4,1));
%save_img(L, sprintf('watershed_channel%02d.tif',4));
%% Deconvolve with Gaussian filter
time = datestr(datetime(now,'ConvertFrom','datenum'));
disp(sprintf('%s ----- deconvolving ...',time));
deconv_stack = zeros(size(norm_stack));

FILTER_SIZE_XY = 6; % for [fxf] filter 
FILTER_SIZE_Z = 6; 
STDXY1 = 1;%4 
STDXY2 = 4;%5 % btw, spots are black not white if std1 > std2
STDZ1 = 1; 
STDZ2 = 4; 
gauss1 = gauss_filter_fun_xyz(FILTER_SIZE_XY, FILTER_SIZE_Z, STDXY1, STDZ1);
gauss2 = gauss_filter_fun_xyz(FILTER_SIZE_XY, FILTER_SIZE_Z, STDXY2, STDZ2);
dgauss = gauss1 - gauss2;

for cycle=1:num_cycles
    for channel=1:num_channels

            img_filtered = imfilter(norm_stack(:,:,:,channel,cycle), dgauss, 'symmetric', 'conv', 'same'); 

            deconv_stack(:,:,:,channel,cycle) = img_filtered;

    end
end


deconv_stack(deconv_stack<0) = 0;

disp(sprintf('%s: Deconvolved images with Gaussian filter',sec2time(toc)));

disp(sprintf('%s: Saving deconvolved images',sec2time(toc)));

deconv_dir = sprintf('%s/deconv',embryo_dir);
 
for cycle=1:num_cycles
     for channel=1:num_channels
         filename = sprintf('%s/cy%02d_ch%02d.tif',deconv_dir,cycle,channel);
         write_3d_tif(filename, deconv_stack(:,:,:,channel,cycle));
     end
end
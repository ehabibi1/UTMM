%%
expt = 'Mechanics_Transcriptome_May2021';
addpath(genpath('\\sodium\brod_thechenlab\ehsan/analysis/ExSeq_Zack/Mechanics_Transcriptome_May2021/'))
%server
%addpath(genpath('/broad/thechenlab/ehsan/analysis/ExSeq_Zack/Mechanics_Transcriptome_May2021/'))

image_bounderies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Processing data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
plate = 1;
well = 3;
embryos_in_fovs = ["2+3"];

for i=1:length(embryos_in_fovs)
    embryos_in_fov = embryos_in_fovs(i);
    numberOfembryo = split(embryos_in_fov,"+");
    if length(numberOfembryo) == 1
        embryo = str2num(str2mat(numberOfembryo));
        process_fov_v1_WS(expt,plate,well,embryos_in_fov,embryo, bounderies);
        process_reg_hyb_to_morpho_WS(expt,plate,well,embryos_in_fov,embryo, bounderies);
        norm_deconv_peakcalling_WS(expt,plate,well,embryo, bounderies);
    else
        for embryo=str2num(str2mat(numberOfembryo))
            process_fov_v1_WS(expt,plate,well,embryos_in_fov,embryo, bounderies);
            process_reg_hyb_to_morpho_WS(expt,plate,well,embryos_in_fov,embryo, bounderies);
            norm_deconv_peakcalling_WS(expt,plate,well,embryo, bounderies);
        end
    end
end

plate = 1;
well = 2;
embryos = [1 2 3 4 5 6 7]; 
embryos_in_fovs = ["1" "2" "3+4+5" "6" "7"];

for i=1:length(embryos_in_fovs)
    embryos_in_fov = embryos_in_fovs(i);
    numberOfembryo = split(embryos_in_fov,"+");
    if length(numberOfembryo) == 1
        embryo = str2num(str2mat(numberOfembryo));
        process_fov_v1(expt,plate,well,embryos_in_fov,embryo, bounderies);
        process_reg_hyb_to_morpho(expt,plate,well,embryos_in_fov,embryo, bounderies);
        norm_deconv_peakcalling(expt,plate,well,embryo, bounderies);
    else
        for embryo=str2num(str2mat(numberOfembryo))
            process_fov_v1(expt,plate,well,embryos_in_fov,embryo, bounderies);
            process_reg_hyb_to_morpho(expt,plate,well,embryos_in_fov,embryo, bounderies);
            norm_deconv_peakcalling(expt,plate,well,embryo, bounderies);
        end
    end
end

plate = 1;
well = 3;
embryos = [1 2 3 4 5 6 7 8]; 
embryos_in_fovs = ["1" "2+3" "4" "5+6" "7" "8"];

for i=1:length(embryos_in_fovs)
    embryos_in_fov = embryos_in_fovs(i);
    numberOfembryo = split(embryos_in_fov,"+");
    if length(numberOfembryo) == 1
        embryo = str2num(str2mat(numberOfembryo));
        process_fov_v1(expt,plate,well,embryos_in_fov,embryo, bounderies);
        process_reg_hyb_to_morpho(expt,plate,well,embryos_in_fov,embryo, bounderies);
        norm_deconv_peakcalling(expt,plate,well,embryo, bounderies);
    else
        for embryo=str2num(str2mat(numberOfembryo))
            process_fov_v1(expt,plate,well,embryos_in_fov,embryo, bounderies);
            process_reg_hyb_to_morpho(expt,plate,well,embryos_in_fov,embryo, bounderies);
            norm_deconv_peakcalling(expt,plate,well,embryo, bounderies);
        end
    end
end

plate = 1;
well = 1;
for embryo=1:8
    plot_peak_threshold(expt, plate, well, embryo)
end

plate = 1;
well = 2;

for embryo=1:7
    plot_peak_threshold(expt, plate, well, embryo)
end



plate = 1;
well = 3;
for embryo=1:8
    plot_peak_threshold(expt, plate, well, embryo)
end

plate = 2;
well = 4;

for embryo=1:7
    plot_peak_threshold(expt, plate, well, embryo)
end

plate = 2;
well = 5;
for embryo=1:7
    plot_peak_threshold(expt, plate, well, embryo)
end



%plate 1 well 1 embryo 1
plate 1 well 2 embryo 5
plate 1 well 1 embryo 4
plate 1 well 1 embryo 2
plate 1 well 3 embryo 1 2 cell

plate 1 well 3 embryo 3 8 cell 
plate_2_well_5_embryo_2 2 cell 

plate = 1;
well = 1;
tic
for embryo=2:8
    cell_by_gene_matrix(expt, plate, well, embryo)
end
toc


plate = 1;
well = 2;
tic
for embryo=1:6%redo plate_1_well_2_embryo_3
    cell_by_gene_matrix(expt, plate, well, embryo)
end
toc

plate = 1;
well = 3;
tic
for embryo=[1:6 8]
    cell_by_gene_matrix(expt, plate, well, embryo)
end
toc

plate = 2;
well = 4;
tic
for embryo=1:7
    cell_by_gene_matrix(expt, plate, well, embryo)
end
toc

plate = 2;
well = 5;
tic
for embryo=1:5
    cell_by_gene_matrix(expt, plate, well, embryo)
end
toc







plate = 2;
well = 5; 
embryos_in_fov = "7";
embryo = 7;
process_reg_hyb_to_mito_WS(expt,plate,well,embryos_in_fov,embryo, bounderies);


plate = 2;
well = 4; 
embryos_in_fov = "1";
embryo = 1;
process_reg_hyb_to_mito_WS(expt,plate,well,embryos_in_fov,embryo, bounderies);


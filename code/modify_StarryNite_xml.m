function[] = modify_StarryNite_xml(home_dir, species, expt,plate, well, embryo)

    file_path = sprintf('%s/projects/%s/%s/transcriptome/segmentation/reg_to_hyb_dapi/plate_%d/well_%d/embryo_%d/temp/embryo_%d',home_dir,species,expt, plate, well, embryo, embryo);
    data_name = sprintf('p%dw%de%d',plate,well,embryo);
    fid = fopen(sprintf('%s/%s_emb1.xml',file_path,data_name), 'r');
    f = fread(fid, '*char')';
    fclose(fid);
    
    string = '<start index="1"/>' + "\n" + '<end index="2"/>';
    string = compose(string);
    f = strrep(f, '<end index="1"/>', string); % to remove
    
    % save into new txt file 
    fid = fopen(sprintf('%s/%s_emb1_edited.xml',file_path,data_name), 'w');
    fprintf(fid, '%s', f);
    fprintf(fid, '%s', f);
end
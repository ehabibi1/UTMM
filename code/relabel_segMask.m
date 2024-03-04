function [seg_stack_corrected labels_table] = relabel_segMask(stackA)
    
    labels = unique(stackA(:));
    
    A = labels(labels~=0);
    num_labels = length(A);
    A_tmp = (1:length(A))';
    for z = 1:size(stackA,3)
        seg_stack_corrected(:,:,z) = changem(stackA(:,:,z),A_tmp',A');
    end 
    
    labels_table = table(A_tmp, A);
    labels_table.Properties.VariableNames = {'NewLabel' 'OriginalLable'}
end
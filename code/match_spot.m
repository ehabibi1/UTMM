function[match_ind match_val hamming_dist close_umis] = match_spot(curr_umi,umis,intensity_mat,hamming_thresh)

    % set consensus UMI
    
    %curr_umi = char(consensus_umis(i,:));

    % fill out Hamming distance matrix
    
    hamming = zeros(length(umis),1);
    for j=1:length(curr_umi)
        hamming = hamming + (umis(:,j)~=curr_umi(j));
    end

    % find UMIs with Hamming distance less than threshold
    
    close_umi_indices = find(hamming<=hamming_thresh); k=1;
    while size(close_umi_indices,1) == 0
        close_umi_indices = find(hamming<=hamming_thresh+k); k=k+1;
    end
    close_umis = umis(close_umi_indices,:);
    
    % calculate probabilities of close UMIs
    
    probs = ones(size(close_umis,1),1)';
    for j=1:length(curr_umi)
        bases = str2num(close_umis(:,j));
        probs = probs.*(intensity_mat(bases,j+1).^2)';
    end

    % find most likely UMI(s)
    
    [min_val min_ind] = min(-log10(probs));
    
    match_ind = close_umi_indices(min_ind);
    match_val = min_val;
    hamming_dist = hamming(close_umi_indices(min_ind));

end
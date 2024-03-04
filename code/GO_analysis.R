# Create a data frame containing the counts of genes in your cluster that fall into each pathway
# Load the necessary libraries
library(dplyr)
library(tidyr)
library(tibble)

pathway_df <- data.frame(
  pathway_id = character(),
  category_1 = character(),
  category_2 = character(),
  pathway_name = character(),
  genes = list(),
  stringsAsFactors = FALSE
)

# Retrieve all pathway IDs in the "Environmental Information Processing" category
pathway_ids <- rownames(as.data.frame(keggList("pathway", "mmu")))
# Initialize list to store pathways
pathways <- list()
CLASS_interests <- c("Genetic Information Processing","Environmental Information Processing","Cellular Processes")


for (pid in pathway_ids) {
  pathway <- keggGet(pid)
  
  tryCatch({
    
    if(is.null(pathway[[1]]$CLASS)){
      print(paste0(pid, " is a M, ... skip"))
    }else if (strsplit(pathway[[1]]$CLASS,"; ")[[1]][1] %in% CLASS_interests){
      # Extract the class and name of the pathway
      class <- pathway[[1]]$CLASS[[1]]
      pathway_name <- pathway[[1]]$NAME
      
      sub_categories <- strsplit(class, "; ")[[1]]
      
      # Check if the pathway has any genes
      if (length(pathway[[1]]$GENE) == 0) {
        genes <- list()
      } else {
        # Extract the genes from the pathway
        genes <- unlist(lapply(pathway[[1]]$GENE[grep(";", pathway[[1]]$GENE)], function(x) strsplit(x, "; ")[[1]][1]))
      }
      
      # Add the pathway to the data frame
      pathway_df <- rbind(pathway_df, data.frame(pathway_id = pid, category_1=sub_categories[1], category_2 = sub_categories[2], pathway_name = pathway_name, genes = genes))
    }
    
  }, error = function(e) {
    if (grepl("HTTP 404", e$message)) {
      message(paste0("Pathway ", pid, " not found, skipping."))
    } else {
      message(paste0("Error for pathway ", pid, ": ", e$message))
    }
  })
}

# Convert the genes column to a comma-separated string
pathway_df$genes <- sapply(pathway_df$genes, paste, collapse = ",")


# Retrieve all pathway IDs in the "Environmental Information Processing" category
pathway_ids <- rownames(as.data.frame(keggList("pathway", "mmu")))
# Initialize list to store pathways
pathways <- list()


CLASS_interests <- c("Genetic Information Processing","Environmental Information Processing","Cellular Processes")


# Loop through each pathway ID and extract category and gene list
for (pid in pathway_ids){
  tryCatch({
    
    # Get pathway information
    pathway <- keggGet(pid)
    
    # Check if the pathway has a CLASS attribute
    if(is.null(pathway[[1]]$CLASS)){
      print(paste0(pid, " is a M, ... skip"))
    }else if (strsplit(pathway[[1]]$CLASS,"; ")[[1]][1] %in% CLASS_interests){
      print (pid)
      # Get the pathway CLASS
      pathway_class <- pathway[[1]]$CLASS
      
      # Loop through each pathway CLASS and add pathway to corresponding nested list
      for (class in pathway_class) {
        
        # Split the CLASS into sub-categories
        sub_categories <- strsplit(class, "; ")[[1]]
        
        # Initialize nested list with top-level category if it doesn't exist
        if (!(sub_categories[1] %in% names(pathways))) {
          pathways[[sub_categories[1]]] <- list()
        }
        
        # Initialize nested list with second-level category if it doesn't exist
        if (length(sub_categories) > 1) {
          if (!(sub_categories[2] %in% names(pathways[[sub_categories[1]]]))) {
            pathways[[sub_categories[1]]][[sub_categories[2]]] <- list()
          }
        }
        
        # Get pathway information
        genes <- unlist(lapply(pathway[[1]]$GENE[grep(";", pathway[[1]]$GENE)], function(x) strsplit(x, "; ")[[1]][1]))
        pathway_name <- pathway[[1]]$NAME
        
        # Add pathway to corresponding nested list
        if (length(sub_categories) == 1) {
          pathways[[sub_categories[1]]][[pid]] <- list(name = pathway_name, genes = genes)
        } else {
          pathways[[sub_categories[1]]][[sub_categories[2]]][[pid]] <- list(name = pathway_name, genes = genes)
        }
      }
    }
  }, error = function(e) {
    if (grepl("HTTP 404", e$message)) {
      message(paste0("Pathway ", pid, " not found, skipping."))
    } else {
      message(paste0("Error for pathway ", pid, ": ", e$message))
    }
  })
}

enrigh_go <- function(df, Clusters=NULL, Genes = NULL){
  # Define function to summarize each row
  summarize_row <- function(row) {
    # Get the row name
    name <- rownames(row)
    
    # Get the non-zero values
    values <- row[row != 0]
    
    # Get the column names for the non-zero values
    genes <- names(values)
    
    # Combine the values and genes into a string
    if (length(genes) == 0) {
      summary_str <- paste0(values, collapse = "")
    } else {
      summary_str <- paste0(length(values), " ", paste(genes, collapse = ","))
    }
    
    # Return the summary as a named vector
    return(c(summary_str, name = name))
  }
  
  if(!is.null(Clusters)){
    genes_in_cluster <- rownames(subset(df, df$cluster %in% Clusters))
    tmp <- pathway_df[which (pathway_df$genes %in% genes_in_cluster),]
  }else if(!is.null(Genes)){
    genes_in_cluster <- Genes
    tmp <- pathway_df[which (pathway_df$genes %in% genes_in_cluster),]  
  }
  
  data_table <- table(tmp$pathway_name, tmp$gene)
  
  
  
  # Apply the function to each row of the data frame
  summary_df <- t(apply(data_table, 1, summarize_row))
  
  # Convert the result to a data frame
  summary_df <- data.frame(t(summary_df))
  
  colnames(summary_df) <- "summary"
  # separate the summary column into two columns
  
  summary_df <- separate(summary_df, summary, into = c("count", "genes"), sep = "\\s+", extra = "drop", fill = "right")
  
  summary_df$count <- as.numeric(summary_df$count)
  summary_df <- summary_df[order(summary_df$count, decreasing = T),]
  
  summary_df <- rownames_to_column(summary_df, "pathway")
  return(summary_df)
}

# Loop through each cluster (e.g. cluster 11, cluster 12, ...)

pathwayAnalysis <- function(data, clusters){
  
  cluster_counts_list <- list()  
  for (i in clusters){ # replace ... with the last cluster number
    # Calculate the pathway counts for the current cluster (e.g. using the function you wrote)
    print(i)
    current_cluster_counts <- enrigh_go(data, Clusters = i)
    if(nrow(current_cluster_counts)>0){ 
      # Filter out pathways with zero gene counts
      current_cluster_counts <- current_cluster_counts %>%
        filter(count > 0)
      # Add the current cluster number to the data frame
      current_cluster_counts$cluster <- i
      
      # Store the current cluster counts in the cluster_counts_list
      cluster_counts_list[[i]] <- current_cluster_counts[,c("pathway","count","cluster", "genes")]
    }
  }
  
  # Combine the pathway counts for all clusters into a single data frame
  cluster_counts <- bind_rows(cluster_counts_list)
  
  # Print the resulting data frame
  print(cluster_counts)
  
  cluster_counts <- cluster_counts %>%
    rename(cluster_count = count)
  
  # Group the cluster_counts data frame by the cluster column
  cluster_groups <- cluster_counts %>% group_by(cluster)
  
  # Define a function to perform the enrichment test for each group
  enrichment_test <- function(group) {
    
    # Get the pathway counts for the current group
    group_pathway_counts <- pathway_counts[pathway_counts$pathway %in% group$pathway,]
    
    # Calculate the total number of genes in the pathways
    total_genes <- sum(group_pathway_counts$count)
    total_cluster_genes <- sum(group$cluster_count)
    # Merge the pathway counts with the current group
    merged_counts <- group_pathway_counts %>%
      left_join(group, by = "pathway") %>%
      mutate(count.x = ifelse(is.na(cluster_count), 0, cluster_count),
             non_cluster_count = count - count.x)
    
    # Perform the Fisher's exact test
    
    enriched_counts <- merged_counts %>%
      mutate(pvalue = apply(., 1, function(x) {
        cluster_count <- as.numeric(x[3])
        non_cluster_count <- as.numeric(x[6])
        fisher.test(matrix(c(cluster_count, non_cluster_count, total_cluster_genes - cluster_count, total_genes - non_cluster_count), nrow = 2))$p.value
      }))
    
    return(enriched_counts)
    
  }
  
  # Apply the enrichment test function to each group and get the p-values
  #p_values <- apply(cluster_groups, 1, enrichment_test)
  
  pvalues <- list()
  for(i in clusters){
    group <- subset(cluster_groups, cluster_groups$cluster==i)
    pvalues[[i]] <- enrichment_test(group)
  }
  
  results_all <- do.call("rbind", pvalues)
  # Adjust the p-values for multiple testing using the Benjamini-Hochberg method
  results_all$adjusted_pvalue <- p.adjust(results_all$pvalue, method = "BH")
  results_all <- results_all[order(results_all$pvalue, decreasing = F),]
  results_all <- subset(results_all, results_all$pvalue<1)
  results_all <- results_all[order(results_all$cluster), ]
  
  
  # Adjust the p-values for multiple testing using the Benjamini-Hochberg method
  merged_counts$adjusted_pvalue <- p.adjust(merged_counts$pvalue, method = "BH")
  merged_counts <- merged_counts[order(merged_counts$pvalue, decreasing = F),]
  
  significant_pathways <- merged_counts %>%
    filter(pvalue < 0.05) %>%
    pull(pathway)
  
  
  list_tmp <- list(cluster_counts, results_all, merged_counts)
  names(list_tmp) <- c("cluster_counts", "results_all", "merged_counts")
  return(list_tmp)
}

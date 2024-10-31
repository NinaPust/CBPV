require(tidyverse)
require(ape)

lsd_rates <- function(data, path_iqtree, input_file, output_file_ghost, tree_file, output_file_lsd){
  
  best_tree <- read.tree(tree_file)
  best_tree <- midpoint(best_tree)
 
    phylo_edge <- edge_extract(best_tree) 
    edge_tibble <- left_join(data, phylo_edge, by = "accession_number") 
    tips_n_times <- edge_tibble[, c("label", "Collection_Date")]
    tips_n_times <- tips_n_times %>% drop_na("Collection_Date")
    write.table(tips_n_times[-1, ], file = (paste0(output_file_ghost, 
                                                   "_dates_tree", ".txt")), 
                sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    midpoint_tree <- midpoint(best_tree)
    write.tree(midpoint_tree, file = paste0(output_file_ghost, "_midpoint_tree", ".treefile"))
    
    lsd <- paste0(path_iqtree, " -s ", input_file," --date ", output_file_ghost, "_dates_tree", ".txt", " -te ",output_file_ghost, "_midpoint_tree", ".treefile -pre ", output_file_lsd, "_tree_", " -redo")
    system(lsd) 
    
    
    lsd_tree <- read.tree(paste0(output_file_lsd, "_tree_", ".timetree.nwk"))
    
    
     lsd_edge <- edge_extract(lsd_tree)
    data_edge <- left_join(edge_tibble, lsd_edge [, c("Edge_length", "accession_number")], 
                           by="accession_number")
    data_edge <- data_edge %>%
      rename_with(~paste0("edge.length.lsd"), matches("Edge_length.y")) %>%
      rename_with(~paste0("edge.length.phylo"), matches("Edge_length.x")) %>%
      mutate(edge.length.lsd = ifelse(edge.length.lsd==0.00000e+00, 1e-5, edge.length.lsd)) %>%
      mutate(!!paste0("rates") := !!sym(paste0("edge.length.phylo")) / !!sym(paste0("edge.length.lsd")))
    data_edge <- dplyr::select(data_edge, - Tip)
  
  return(data_edge)
  
}

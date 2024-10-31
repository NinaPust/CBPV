require(ape)
require(tidyverse)
require(MASS)
require(rentrez)

download_align_trimm <- function(data, align_input_file, align_output_file1, align_output_file2, trimmed_output_file, path_to_macse, memory = "1G") {
  accession_numbers <- data$accession_number
  
  batches <- split(accession_numbers, ceiling(seq_along(accession_numbers) / 100))
  
  DNA_bin <- NULL
  for (batch in batches) {
    # Read sequences from GenBank for the current batch
    result <- read.GenBank(batch)
    
    # Combine the results into the main DNAbin object
    if (is.null(DNA_bin)) {
      DNA_bin <- result # Initialize if it's the first batch
    } else {
      DNA_bin <- c(DNA_bin, result) # Concatenate to grow DNAbin object
    }
  }
  
  DNA_nondups <- DNA_bin[!duplicated(DNA_bin)]
  DNA_dups <- DNA_bin[!(names(DNA_bin)) %in% names(DNA_nondups)]
  
  data <- data[!(data$accession_number %in% names(DNA_dups)), ]
  duplicated_accessions <- names(DNA_dups)
  
  write.FASTA(DNA_nondups, file = paste0(align_input_file))
  command<- paste0("java -Xmx", memory, " -jar ", path_to_macse, " -prog alignSequences -max_refine_iter 0 -seq ", align_input_file, " -out_NT ", align_output_file1)
  system(command)

  command2 <- paste0("java -Xmx", memory, " -jar ", path_to_macse, " -prog exportAlignment -align ", align_output_file1, " -codonForInternalStop NNN -codonForFinalStop --- -codonForInternalFS --- -charForRemainingFS - -out_NT ", align_output_file2)
  system(command2)
 
  aligned_DNA <- read.FASTA(paste0(align_output_file2))
  aligned_DNA <- as.matrix(aligned_DNA)
  
  fasttree <- njs(dist.dna(aligned_DNA))
  fasttree$edge.length <- abs(fasttree$edge.length) #absolute because min value was negative
  
  trimmedtree <- trimLongTipBrs(fasttree, pval = 0.001)
  
  aligned_DNA <- aligned_DNA[trimmedtree$tip.label,]
  
  aligned_DNA_trimmed <- trimCols(aligned_DNA, prop = 0.5, codon = T)
  
  data <- data[!(data$accession_number %in% names(aligned_DNA_trimmed)), ]
  
  write.FASTA(aligned_DNA_trimmed, file = paste0(trimmed_output_file))
  
  if(length(duplicated_accessions)>0) {
    warning <- paste0("Samples with accession numbers ", duplicated_accessions, " were removed due to being duplicated.")
    print(warning)}
  else{
    print("There were no duplicated sequences.")
  }
  
  return(data)
  
}

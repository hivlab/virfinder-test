
# split fasta to sequences
insert_newlines_concat <- function(x) {
  x[1] <- paste0(x[1], "\n")
  x[length(x)] <- paste0(x[length(x)], "\n")
  paste(x, collapse = "")
}

#' Subset fasta file by sequence ids
#' @param ids fasta ids without '>', character vector
#' @param infile path to fasta file, character string
#' @param outfile path to output fasta file, character string
#' @example subset_fasta(ids, "fasta/assembly_7.fa", "output/assembly_7_ids.fa")
#' 
subset_fasta <- function(ids, infile, outfile) {
  fa <- readLines(infile)
  i <- grep(">", fa)
  seqlength <- diff(c(i, length(fa) + 1))
  splitvec <- rep(1:length(seqlength), seqlength)
  fasplit <- split(fa, splitvec)
  names(fasplit) <- fa[i]
  fasplit <- lapply(fasplit, insert_newlines_concat)
  fasta_txt <- paste(fasplit[sub(">", "", names(fasplit)) %in% ids], collapse = "")
  fileConn <- file(outfile, "w")
  writeLines(fasta_txt, fileConn)
  close(fileConn)
}


#' Selecting SNPs located in a specific list of genes
#'
#' \code{Map_genestoSNPs} takes a specified gene set and returns a list of SNPs from a genotype data set located within these genes
#'  including n bases at both edges of the coding region. It also generates a file with a list of genes from the gene set
#'  not found in the hg19 (NCBI37.3) version of the human genome. This is useful to check for the correct name.
#'
#' @param geneset A dataframe of one column or an array with genes from gene set, HGNC symbol or NCBI gene ID, no header.
#' @param bfile A string indicating the name of the PLINK binary file with data of the GWAS target sample, without file extension.
#' @param path.to.plink.files String with the full path to the three PLINK binary files (.bed, .bim, .fam)
#' @param extra.kb An array of two numbers indicating the extra bases in kb to add at the 5' and 3' edges of the
#'   coding region for mapping SNPs to genes. Default extra.kb=c(0,0).
#' @param output A string indicating the name of the output file to write the list of genes not identified. This file
#'   is generated in the path defined in \code{path.to.plink.files}.
#'
#' @return A dataframe of one column with the SNPs mapping your genes, identified by chromosome and position (ex. 4:103188709).

Map_genestoSNPs <- function(geneset, bfile, path.to.plink.files, extra.kb=c(0,0), output = "missing_genes.txt") {

  # If list of genes not a dataframe
  if (!is.data.frame(geneset)) {
    geneset <- as.data.frame(geneset)
  }

  # Genes in gene set, first, assuming genenames as Gene_ID (column 1 of hg19)
  pos.genes <- subset(hg19, hg19$V1 %in% geneset[, 1])
  # Check for undetected genes
  if (dim(pos.genes)[1] != 0 & setequal(geneset[, 1], pos.genes$V1) != T){
    N<- length(subset(geneset[, 1], !geneset[,1] %in% pos.genes$V1))
    message("A total of ",N," genes from gene set were not found (written to output file)")
    # Writing missing genes to output file
    missing <- subset(geneset[, 1], !geneset[, 1] %in% pos.genes$V1)
    missing.file.name <- paste0(path.to.plink.files, output)
    write.table(x = missing, file = missing.file.name, append = F, quote = F, row.names = F, col.names = F)
  }
  # Then, assuming genenames as HGCN_Symbol (column 6 of hg19)
  if(dim(pos.genes)[1] == 0) {
    pos.genes <- subset(hg19, hg19$V6 %in% geneset[, 1])
    # Check for undetected genes
    if ( dim(pos.genes)[1] != 0 & setequal(geneset[, 1], pos.genes$V6) != T) {
      N <- length(subset(geneset[, 1], !geneset[, 1] %in% pos.genes$V6))
      message("A total of ",N," genes from gene set were not found (written to output file)")
      # Writing missing genes to output file
      missing <- subset(geneset[, 1], !geneset[, 1] %in% pos.genes$V6)
      missing.file.name <- paste0(path.to.plink.files, output)
      write.table(x = missing, file = missing.file.name, append = F, quote = F, row.names = F, col.names = F)
    }
  } 
  if (dim(pos.genes)[1] == 0) {
    stop("Gene identifiers not in database")
    }

  #Adding extra bases
  if(sum(extra.kb) == 0) {
    gene.limits <- pos.genes[, c(2:4)]
    colnames(gene.limits) <- c("chrom", "start", "end")
    gene.limits$chrom <- paste0("chr", gene.limits$chrom)
  } else if (sum(extra.kb) != 0) {
   forward <- subset(pos.genes, pos.genes$V5 == "+")
   forward$V3 <- forward$V3 - extra.kb[1] * 1000
   forward$V4 <- forward$V4 + extra.kb[2] * 1000
   reverse <- subset(pos.genes, pos.genes$V5 == "-")
   reverse$V3 <- reverse$V3 - extra.kb[2] * 1000
   reverse$V4 <- reverse$V4 + extra.kb[1] * 1000
   gene.limits <- rbind(forward, reverse)
   gene.limits <- gene.limits[, c(2:4)]
   colnames(gene.limits) <- c("chrom", "start", "end")
   gene.limits$chrom <- paste0("chr", gene.limits$chrom)
  }

  # Creating GRanges object with intervals
  int.gr <- as(gene.limits, "GRanges")

  #Preparing input SNPs GRanges
  bfile.path<-paste0(path.to.plink.files,bfile,".bim")
  SNPs <- read.table(bfile.path,header = F)
  SNPs$chrom <- paste0("chr", SNPs$V1)
  SNPs$start <- SNPs$V4
  SNPs$end <- SNPs$V4
  SNPs$names <- SNPs$V2
  SNPs.format <- SNPs[, 7:10]

  # Creating GRanges object with intervals
  snps.gr <- as(SNPs.format, "GRanges")

  # Overlaps
  ov <- GenomicRanges::findOverlaps(int.gr, snps.gr)
  snps.hit <- as.data.frame(ov)[, 2]
  list <- SNPs.format[snps.hit, ]

  # Transforming to input for Calculate_DR2Nagelkerke function
  list$chrom <- gsub(pattern = "chr", replacement = "", x = list$chrom)
  list$V1 <- paste(list$chrom, list$start, sep=":")
  list.pos <- as.data.frame(list$V1)
  colnames(list.pos) <- "V1"
  if (dim(list.pos)[1] == 0){
    message("\nNo SNP maps within your genes!" )
  }

  return(list.pos)
}



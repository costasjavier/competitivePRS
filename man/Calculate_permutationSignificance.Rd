% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Calculate_permutationSignificance.R
\name{Calculate_permutationSignificance}
\alias{Calculate_permutationSignificance}
\title{Calculating significance by permutations}
\usage{
Calculate_permutationSignificance(path.to.plink.files, bfile = "tmp", sumGWAS,
  Perm, Cov = NULL, DR2, cov.names, outfile)
}
\arguments{
\item{path.to.plink.files}{String with the full path to the three PLINK binary files (.bed, .bim, .fam).}

\item{bfile}{A string indicating the name of the PLINK binary file with data of the GWAS target sample, without extension.
The tmp files generated by the \code{Calculate_DR2Nagelkerke} function are used by default. If other files are used, the "snpid" column 
has to be identified as chr:pos (such as 4:103334).}

\item{sumGWAS}{A dataframe with a minimum of these columns "chr.pos", "a1", "p", "logOR", taken from the summary
statistics of the GWAS discovery sample.}

\item{Perm}{A dataframe with the permuted SNP sets, generated by the \code{Generate_permutationSets} function of this package.}

\item{Cov}{A dataframe with a column "IID"  for sample identification and covariates with specific names.}

\item{DR2}{A number corresponding to the value of the increase in pseudo-R2 of Naglekerke of the actual SNP set. It is the
first element of the list generated by the \code{calculate_DR2Nagelkerke} function of this package.}

\item{cov.names}{An array with the names of covariates from the \code{Cov} dataframe to be included in the model.}

\item{outfile}{A string indicating the name of the output file to write the P value. This file
is generated in the path defined in \code{path.to.plink.files}.}
}
\value{
A list with two elements. The first element is a number representing the permuted P value. The second element is
   an array with the pseudo-R2 values of the permuted data sets.
}
\description{
\code{Calculate_permutationSignificance} estimates the significance of the actual PRS based on permutations of
  random sets of SNPs matching the actual SNP set by several properties.The random sets are generated by the
  \code{Generate_permutationSets} function of this package. The result is saved in a file.
}

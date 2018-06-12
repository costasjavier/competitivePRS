## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_knit$set(root.dir = normalizePath(".."))

## ----eval = FALSE--------------------------------------------------------
#  devtools::install_github("costasjavier/competitivePRS",install_vignettes = TRUE)

## ----message = FALSE, warning = FALSE------------------------------------
library(GenomicRanges)
library(rms)
library(competitivePRS)

## ------------------------------------------------------------------------
my.gene.list <- c("COMT", "PRODH", "TTC28", "MTMR3", "FFFF")
my.gene.list

## ------------------------------------------------------------------------
path <- system.file("extdata","exGenelist.txt",package = "competitivePRS")

## ------------------------------------------------------------------------
my.gene.list <- read.table(path, header = F)
head(my.gene.list)

## ------------------------------------------------------------------------
path <- system.file("extdata","",package = "competitivePRS")
path<-paste0(path,"/")

## ------------------------------------------------------------------------
my_SNP_list <- Map_genestoSNPs(geneset = my.gene.list, path.to.plink.files = path, bfile = "exPLINK", extra.kb = c(35,10), output = "my_output_1")
head(my_SNP_list)

## ------------------------------------------------------------------------
path1 <- paste(path,"exSummaryGWAS",sep = "/")
ex_summGWAS <- read.table(path1, header = T)
head(ex_summGWAS)

## ------------------------------------------------------------------------
path2 <- paste(path,"exCovariates",sep="/")
Covariates <-read.table(path2,header = T)
head(Covariates)

## ------------------------------------------------------------------------
PRS_analysis <- Calculate_DR2Nagelkerke(SNPs.list = my_SNP_list, sumGWAS = ex_summGWAS, Cov = Covariates,      path.to.plink.files = path, bfile = "exPLINK", P.threshold = 0.5, amb.remove = T,
cov.names = c("Sex","PCA2"), better.coef = c(0.05,0.5), output.file = "my_output_2")

## ------------------------------------------------------------------------
PRS_analysis$`Pseudo-R2`

## ------------------------------------------------------------------------
PRS_analysis$Input_perm

## ------------------------------------------------------------------------
PRS_analysis$Logistic_model

## ------------------------------------------------------------------------
PRS_analysis$Better_coef

## ------------------------------------------------------------------------
random_sets <- Generate_permutationSets(input.perm = PRS_analysis[[2]], sumGWAS = ex_summGWAS, bfile = "tmp", path.to.plink.files = path, n.perm = 50)

## ------------------------------------------------------------------------
dim(random_sets)

## ------------------------------------------------------------------------
random_sets[1:9,1:6]

## ------------------------------------------------------------------------
result_test <- Calculate_permutationSignificance(path.to.plink.files = path, sumGWAS = ex_summGWAS, 
Perm = random_sets, Cov = Covariates, cov.names = c("Sex","PCA2"), DR2 = PRS_analysis[[1]], outfile = "my_output_3")
                                                

## ------------------------------------------------------------------------
str(result_test)

## ------------------------------------------------------------------------
Generate_Plot(permDR2 = result_test[[2]], DR2 = PRS_analysis[[1]])


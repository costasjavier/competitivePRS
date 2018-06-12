## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
knitr::opts_knit$set(root.dir = normalizePath(".."))

## ----eval = FALSE--------------------------------------------------------
#  devtools::install_github("costasjavier/competitivePRS")

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


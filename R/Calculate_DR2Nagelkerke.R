#' Calculating the increase of pseudo-R2 of Nagelkerke due to PRS
#'
#' \code{Calculate_DR2Nagelkerke} calculates the increase in pseudo-R2 of Naglekerke between a full logistic regression model,
#'    considering covariates and the PRS, and a model considering only covariates. It also creates
#'    an array with the list of SNPs in the PRS model, identified by chromosome and position (ex. 4:103188709).
#'    It also generates a file with the main data of the analysis, including the proccess of SNP clumping, as well as
#'    the value of increase in pseudo-R2 of Nagelkerke, and the OR (and 95 percent C.I.) of the standardized score. Optionally,
#'    it also calculates a better coefficient of determination than Nagelkerke's pseudo-R2, interpetable on the
#'    liability scale and adjusted for ascertainment bias. See Lee et al (2012).
#'
#' @param SNPs.list A dataframe with one column,indicating the SNPs belonging to the gene set, in the format chr:position.
#'   This is the return of function \code{Map_genestoSNPs}.
#' @param sumGWAS A dataframe with a minimum of four columns with names "snpid", "a1", "p","effect", taken from the summary
#'   statistics of the GWAS discovery sample. SNPs in the "snpid" column has to be identified as chr:pos (such as 4:103334).
#'   The "effect" column indicates the effect attributed to allele "a1", as logOR (beta coefficient).
#' @param amb.remove A logical indicating if SNPs with strand ambiguity (i.e. A<->T and C<->G) have to be removed from sumGWAS
#'   before the estimation of PRS. If \code{TRUE}, the sumGWAS has to include an additional column with name "a2" to indicate
#'   the non-effect allele at each SNP. Default \code{FALSE}.
#' @param Cov A dataframe with a column "IID"  for sample identification and covariates with specific names.
#' @param cov.names An array with the names of covariates from the \code{Cov} dataframe to be included in the model.
#' @param P.threshold A number indicating the P threshold of the GWAS discovery sample to include an SNP in the PRS model. Default P.threshold = 1.
#' @param output.file A string indicating the name of the output file to write the main data of the analysis. This file
#'   is generated in the path defined in \code{path.to.plink.files}
#' @param path.to.plink.files String with the full path to the three PLINK binary files (.bed, .bim, .fam)
#' @param bfile A string indicating the name of the PLINK binary file with data of the GWAS target sample, without extension.
#' @param better.coef An array of two numbers, indicating disease prevalence and proportion of cases in the target sample.
#'   Only needed to calculate a better coefficient of determination than Nagelkerke's pseudo-R2, interpetable on the
#'   liability scale and adjusted for ascertainment bias. See Lee et al (2012).
#' @param new.bfiles A logical indicating if new bfiles have to be generated, replacing SNP_ID by chr:pos and removing duplicated
#'   and triallelic SNPs. If \code{TRUE}, new bfiles with the \code{tmp} name will be generated. Default \code{FALSE}.
#'
#' @return A list with four elements. A number representing the increase in pseudo-R2 of Nagelkerke, an array of SNPs
#'   included in the PRS model, the summary of the logistic regression model (generated by the \code{lrs} function at
#'   the \code{rms} package), and a better coefficient of determination, if requested.
#'

Calculate_DR2Nagelkerke <- function(SNPs.list, sumGWAS, amb.remove = FALSE, new.bfiles = FALSE, Cov = NULL, path.to.plink.files, bfile, P.threshold = 1,
                                    cov.names = NULL, better.coef = NULL, output.file = "output") {

  # Basic check of sumGWAS file
  if (!"snpid" %in% colnames(sumGWAS)) stop("No \"snpid\" column in sumGWAS")
  if (!"a1"  %in% colnames(sumGWAS)) stop("No \"a1\" column in sumGWAS")
  if (!"p"  %in% colnames(sumGWAS)) stop("No \"p\" column in sumGWAS")
  if (!"effect"  %in% colnames(sumGWAS)) stop("No \"effect\" column in sumGWAS")
  if (amb.remove == TRUE & !"a2"  %in% colnames(sumGWAS)) stop("No \"a2\" column in sumGWAS, needed for removal of SNPs with strand ambiguity")
  # Basic check of covariates file
  if (!is.null(Cov)){
    if (!("IID" %in% colnames(Cov))) stop("No \"IID\" column in Cov")
    if (length(subset(cov.names,cov.names %in% colnames(Cov))) != length(cov.names)) stop("Some covariate names in cov.names are not in the Cov dataframe")
    if (is.null(cov.names)) stop ("Covariates names no specified")
  }
  if (is.null(Cov) & !is.null(cov.names)) stop("Covariate names specified but covariate file is lacked")
  
  #--------------------------------
  # Calculation of individual PRS
  #--------------------------------

  if (amb.remove == TRUE){
    sumGWAS$a1a2 <- paste(sumGWAS$a1, sumGWAS$a2)
    sumGWAS$a1a2 <- as.factor(toupper(sumGWAS$a1a2))
    sumGWAS <- subset(sumGWAS, sumGWAS$a1a2 != "AT" & sumGWAS$a1a2 != "CG" & sumGWAS$a1a2 != "TA" & sumGWAS$a1a2 != "GC")
  }

  # Opening .bim file
  bimfile <- paste(paste0(path.to.plink.files, bfile), "bim", sep = ".")
  bim.file <- read.table(bimfile, header = F)

  if (new.bfiles == TRUE){
  # Removing unmmapped SNPs and repeated SNPs
  bim.file$chr.pos <- paste(bim.file$V1, bim.file$V4, sep = ":")
  counts <- as.data.frame(table(bim.file$chr.pos))
  to.remove <- subset(counts, counts$Freq != 1)
  remove <- subset(bim.file, bim.file$chr.pos %in% to.remove$Var1)
  SNPstoRemove <- paste0(path.to.plink.files, "SNPstoRemove")
  write.table(x = remove$V2, file = SNPstoRemove, append = F, quote = F, row.names = F, col.names = F)

  # Running PLINK to generate a bfile with unique mapped SNPs
  arg1 <- paste("--bfile", paste0(path.to.plink.files, bfile), sep=" ")
  arg2 <- paste("--exclude", paste0(path.to.plink.files, "SNPstoRemove"), sep=" ")
  arg3 <- paste("--out", paste0(path.to.plink.files, "tmp"), sep=" ")
  system2("plink", args = c(arg1, arg2, "--make-bed", arg3))
  tmp.bimfile <- paste(paste0(path.to.plink.files, "tmp"), "bim", sep=".")
  bim.tmp <- read.table(tmp.bimfile, header = F)
  bim.tmp$V2 <- paste(bim.tmp$V1, bim.tmp$V4, sep=":")
  write.table(x = bim.tmp, file = tmp.bimfile, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
  }
  #Performing LD-clumping using PLINK
  subset.for.clumping <- subset(sumGWAS, sumGWAS$snpid %in% SNPs.list[, 1])
  if (dim(subset.for.clumping)[1] == 0) {
    stop("\nNone of the SNPs in your list are present in the GWAS summary sample data set! Check if format is correct.")
  }
  # Calculate number of SNPs in list not in sumGWAS to include in report
  n <- dim(SNPs.list)[1] - dim(subset.for.clumping)[1]
  For.clumping <- subset.for.clumping[, c("snpid", "p")]
  For.clumping <- subset(For.clumping, For.clumping$p <= P.threshold)
  #Check for existence of SNPs after Pthreshold application
  if (dim(For.clumping)[1] == 0) stop ("No SNP remains for PRS under the selected P threshold")
  colnames(For.clumping) <- c("SNP", "P")
  write.table(For.clumping,paste0(path.to.plink.files, "For.clumping"), append = F, quote = F, sep = "\t", row.names = F, col.names = T)
  if (new.bfiles == TRUE){
  arg1 <- paste("--bfile", paste0(path.to.plink.files, "tmp"), sep=" ")
    } else if (new.bfiles == FALSE) {
  arg1 <- paste("--bfile", paste0(path.to.plink.files, bfile), sep=" ") 
  }
  arg2 <- paste("--clump", paste0(path.to.plink.files, "For.clumping"), sep=" ")
  arg3 <- paste("--out", paste0(path.to.plink.files, "clumped"), sep=" ")
  system2("plink", args=c(arg1, arg2, "--clump-kb 500", "--clump-p1 1", "--clump-p2 1", "--clump-r2 0.1", arg3))

  # In case of SNPs in linkage equilibirum, PLINK does not perform clumping, check
  path<-paste0(path.to.plink.files, "clumped.CLUMPED")
  if (file.exists(path)) {
    clumped <- read.table(path, header = T)
    subset.PRS <- subset(sumGWAS, sumGWAS$snpid %in% clumped$SNP)
  } else {
    subset.PRS <- subset(sumGWAS, sumGWAS$snpid %in% SNPs.list[, 1])
  }

  #Preparing SNP list for permutations.
  pre.input.perm <- subset.PRS[, "snpid"]
  input.perm <- gsub(pattern = "chr", replacement = "", x = pre.input.perm)

  # Calculating the individual scores
  PRS.model <- subset.PRS[, c("snpid", "a1", "effect")]
  name <- paste0(path.to.plink.files, "PRSmodel")
  write.table(PRS.model, file = name, append = F, quote = F, sep = "\t", row.names = F, col.names = F)
  #arg1 as above
  arg2 <- paste("--score", paste0(path.to.plink.files, "PRSmodel"), sep=" ")
  arg3<- paste("--out", paste0(path.to.plink.files, "Scores"), sep=" ")
  system2("plink", args=c(arg1, arg2, arg3))
  path <- paste0(path.to.plink.files, "Scores.PROFILE")
  Scores <- read.table(path, header = T)

  #----------------------------------------------
  # Statistical analyses by logistic regression
  # --------------------------------------------

   # Merge covariables with score, if any
  if (is.null(Cov)) {
      data <- Scores

      # Calculating % missing

      data$percent_missing <- (max(data$CNT) - data$CNT) / max(data$CNT)

      # Additional covariables
      pheno <- data$PHENO
      score <- data$SCORE
      miss <- data$percent_missing
      # Standardized score
      st.score <- scale(score)

      # If there is no missing, there is an error in logistic regression. Check
      max.missing <- max(miss)
      if (max.missing > 0) {
        H0 <- rms::lrm(pheno ~ miss)
        H1 <- rms::lrm(pheno ~ miss + st.score)
        OR <- exp(coef(H1)[3])
        CI <- exp(confint.default(H1))
        ci <- as.numeric(CI[3, 1:2])
        # Increase in Nagelkerke's R2
        DR2 <- H1$stats[10] - H0$stats[10]
      } else {
        H1 <- rms::lrm(pheno ~ st.score)
        OR <- exp(coef(H1)[1])
        CI <- exp(confint.default(H1))
        ci <- as.numeric(CI[2, 1:2])

        # Increase in Nagelkerke's R2
        DR2 <- (H1$stats[10])
      }
  } else if (!is.null(Cov)) {

    # Choosing covariates
    cov.for.used <- colnames(Cov) %in% cov.names
    used.cov <- as.data.frame(Cov[, cov.for.used])
    if (length(cov.names) == 1) colnames(used.cov) <- cov.names#If there is just one covariate, coercing to dataframe implies lost of colnames
    used.cov$IID <- Cov$IID
    used.cov.ord <- used.cov[order(used.cov$IID), ]
    Scores.ord <- Scores[order(Scores$IID), ]
    for (i in 1:dim(used.cov.ord)[2]) {
      assign(colnames(used.cov.ord)[i], used.cov.ord[, i])
    }

    pheno <- Scores.ord$PHENO
    score <- Scores.ord$SCORE

    # Standardized score
    st.score <- scale(score)
    names <- cov.names[1] 
    if (length(cov.names) > 1){
      for (i in 2:length(cov.names)) {
      names <- paste(names, cov.names[i], sep=" + ")
      }
    }
    # If there is no missing, there was an error in logistic regression. Check
    Scores.ord$percent_missing <- (max(Scores.ord$CNT) - Scores.ord$CNT) / max(Scores.ord$CNT)
    miss <- Scores.ord$percent_missing
    max.missing <- max(miss)
    if (max.missing > 0) {
      formula.m0 <- as.formula(paste("pheno ~ miss", names, sep = " + "))
      formula.m1 <- as.formula(paste("pheno ~ st.score + miss", names, sep =  " + "))
      H0 <- rms::lrm(formula.m0)
      H1 <- rms::lrm(formula.m1)
      OR <- exp(coef(H1)[2])
      CI <- exp(confint.default(H1))
      ci <- as.numeric(CI[2, 1:2])

      # Increase in Nagelkerke's R2
      DR2 <- H1$stats[10] - H0$stats[10]
    } else {
      formula.m0 <- as.formula(paste("pheno ~ ", names,sep = "+"))
      H0 <- rms::lrm(formula.m0)
      formula.m1 <- as.formula(paste("pheno ~ st.score", names, sep = "+"))
      H1 <- rms::lrm(formula.m1)
      OR <- exp(coef(H1)[2])
      CI <- exp(confint.default(H1))
      ci <- as.numeric(CI[2,1:2])

      # Increase in Nagelkerke's R2
      DR2 <- H1$stats[10] - H0$stats[10]
    }
  }

  #Calculating a better coefficient of determination. Following function from https://github.com/kn3in/ABC
  if (!is.null(better.coef)) {
   k <- better.coef[1]
   p <- better.coef[2]
   x <- qnorm(1 - k)
   z <- dnorm(x)
   i <- z / k
   cc <- k * (1 - k) * k * (1 - k) / (z^2 * p * (1 - p))
   theta <- i * ((p - k)/(1 - k)) * (i * ((p - k) / ( 1 - k)) - x)
   e <- 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
   h2l <- cc * e * DR2 / (1 + cc * e * theta * DR2)
  }
  if (is.null(better.coef)) {
    h2l <- "n.a."
  }


  # Writing main results to a file
  file.name <- paste0(path.to.plink.files, output.file)
  sink(file.name)
  cat("The number of SNPs from the list absent from the summary GWAS is:", n, "\n\n")
  cat("The number of SNPs in model after clumping and application of P threshold (if any) is:", dim(PRS.model)[1], "\n\n")
  print(H1) 
  cat("\nThe OR for standardized score is",round(OR,3), "95%C.I.", round(ci[1],3), "-", round(ci[2],3), "\n\n")
  cat("The increase in Nagelkerke's R2 between a model with covariates and a model with covariates and PRS is:",  DR2,"\n\n")
  cat("A better value for the coefficient of determination, interpetable on the liability scale and adjusted for ascertainment bias is:",  h2l)
  sink()

  message("\nResults written to file ", file.name)

  output <- list(DR2, input.perm, H1, h2l)
  names(output)<-c("Pseudo-R2", "Input_perm", "Logistic_model", "Better_coef")

  return(output)
}



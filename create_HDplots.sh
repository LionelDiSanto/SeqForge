echo '#------------------------------------------#' > HDplots.R
echo '# Filtering spurious paralogs using HDplot #' >> HDplots.R
echo '#------------------------------------------#' >> HDplots.R
echo '# This script aims to load VCF files filtered with VCFtools and filter them ' >> HDplots.R
echo '# further to remove potential sites resulting from the mapping of paralogous' >> HDplots.R
echo '# reads using HDplots.' >> HDplots.R
echo '' >> HDplots.R
echo '# Reference' >> HDplots.R
echo '# McKinney, G.J., Waples, R.K., Seeb, L.W. and Seeb, J.E. (2017), Paralogs are' >> HDplots.R
echo '# revealed by proportion of heterozygotes and deviations in read ratios in ' >> HDplots.R
echo '# genotyping-by-sequencing data from natural populations. Mol Ecol Resour, 17: 656-669.' >> HDplots.R
echo '# https://doi.org/10.1111/1755-0998.12613' >> HDplots.R
echo '' >> HDplots.R
echo '#------------------------------------------------#' >> HDplots.R
echo '# Double-check necessary libraries are installed #' >> HDplots.R
echo '#------------------------------------------------#' >> HDplots.R
echo 'lib_path <- Sys.getenv("R_LIBS_USER")' >> HDplots.R
echo 'if (!dir.exists(lib_path)) stop("No repository for R packages found. Make sure you have a personal (local) repository for R packages.")' >> HDplots.R
echo 'is_installed <- function(pkg){' >> HDplots.R
echo '  test <- pkg %in% rownames(installed.packages(lib.loc = lib_path))' >> HDplots.R
echo '  if(isTRUE(test)) cat(paste(" ", pkg, "is installed\n")) else cat(paste(" ", pkg, "is NOT installed. Please install it using the command line install.packages()\n"))' >> HDplots.R
echo '  return(test)' >> HDplots.R
echo '}' >> HDplots.R
echo 'lib_to_test <- c("vcfR", "ggplot2", "dplyr", "stringr")' >> HDplots.R
echo 'pac <- sapply(lib_to_test, is_installed)' >> HDplots.R
echo 'out <- sapply(names(pac), function(i) if(isFALSE(pac[i])) stop(paste("The required R package", i, "need to be installed in a local R package repository.")))' >> HDplots.R
echo 'suppressMessages(require(vcfR, warn.conflicts = F))' >> HDplots.R
echo 'suppressMessages(require(ggplot2, warn.conflicts = F))' >> HDplots.R
echo 'suppressMessages(require(dplyr, warn.conflicts = F))' >> HDplots.R
echo 'suppressMessages(require(stringr, warn.conflicts = F))' >> HDplots.R
echo '' >> HDplots.R
echo '#----------#' >> HDplots.R
echo '# Function #' >> HDplots.R
echo '#----------#' >> HDplots.R
echo 'HDplot<-function(vcfData){' >> HDplots.R
echo '  #set up results table' >> HDplots.R
echo '  HDplotTable<-as.data.frame(matrix(NA,nrow=dim(vcfData@gt)[1],ncol=13))' >> HDplots.R
echo '  colnames(HDplotTable)<-c("CHROM","POS","ID","depth_a","depth_b","ratio","num_hets","num_samples","num_called","H_all","H","std","D")' >> HDplots.R
echo '' >> HDplots.R
echo '  #get genotypes from vcf file' >> HDplots.R
echo '  genos<-extract.gt(vcfData, element = "GT", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, ' >> HDplots.R
echo '                    IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)' >> HDplots.R
echo '' >> HDplots.R
echo '  #replace NA genotypes with ./.' >> HDplots.R
echo '  genos[is.na(genos)]<-"./."' >> HDplots.R
echo '' >> HDplots.R
echo '  #get allele reads from vcf file' >> HDplots.R
echo '  reads<-extract.gt(vcfData, element = "AD", mask = FALSE, as.numeric = FALSE, return.alleles = FALSE, ' >> HDplots.R
echo '                    IDtoRowNames = TRUE, extract = TRUE, convertNA = FALSE)' >> HDplots.R
echo '' >> HDplots.R
echo '  #replace reads for samples with missing data with 0' >> HDplots.R
echo '  #reads<-gsub("\\.,\\.","0,0",reads)' >> HDplots.R
echo '  reads[grepl("\\.",reads)]<-"0,0"' >> HDplots.R
echo '  reads[is.na(reads)]<-"0,0"' >> HDplots.R
echo '' >> HDplots.R
echo '  alleleReads<-apply(reads,2,function(x) str_split_fixed(x,",",2))' >> HDplots.R
echo '  alleleReads_1<-alleleReads[1:dim(reads)[1],]' >> HDplots.R
echo '  alleleReads_2<-alleleReads[dim(reads)[1]+1:dim(reads)[1],]' >> HDplots.R
echo '  #convert to numeric format' >> HDplots.R
echo '  alleleReads_1<-apply(alleleReads_1,2, function(x) as.numeric(x))' >> HDplots.R
echo '  alleleReads_2<-apply(alleleReads_2,2, function(x) as.numeric(x))' >> HDplots.R
echo '  #subset to heterozygous genotypes' >> HDplots.R
echo '  #make genotype matrix where heterozygotes are 1 and other genotypes are 0' >> HDplots.R
echo "  hetMatrix<-apply(genos,2,function(x) dplyr::recode(x,'0/0'=0,'1/1'=0,'./.'=0,'0/1'=1,'1/0'=1))" >> HDplots.R
echo "  calledGenos<-apply(genos,2,function(x) dplyr::recode(x,'0/0'=1,'1/1'=1,'0/1'=1,'1/0'=1,.default=NA_real_))" >> HDplots.R
echo '  #multiply read count matrices by heterozygote matrix to get read counts for heterozygous genotypes' >> HDplots.R
echo '  alleleReads_1_het<-alleleReads_1*hetMatrix' >> HDplots.R
echo '  alleleReads_2_het<-alleleReads_2*hetMatrix' >> HDplots.R
echo '  #rows are loci and columns are samples' >> HDplots.R
echo '  #sum reads per allele per locus for heterozygous samples' >> HDplots.R
echo '  A_reads<-apply(alleleReads_1_het,1,sum)' >> HDplots.R
echo '  B_reads<-apply(alleleReads_2_het,1,sum)' >> HDplots.R
echo '  totalReads<-A_reads+B_reads' >> HDplots.R
echo '  ratio<-A_reads/totalReads' >> HDplots.R
echo '  std<-sqrt(totalReads*0.5*0.5)' >> HDplots.R
echo '  z<- -(totalReads/2-A_reads)/std' >> HDplots.R
echo '  #get percent heterozygosity for each locus' >> HDplots.R
echo '  numHets<-apply(hetMatrix,1,sum)' >> HDplots.R
echo '  hetPerc<-numHets/dim(hetMatrix)[2]' >> HDplots.R
echo '' >> HDplots.R
echo '  numGenos<-apply(calledGenos,1,sum,na.rm=TRUE)' >> HDplots.R
echo '  H<-numHets/numGenos' >> HDplots.R
echo '' >> HDplots.R
echo '  #assign results to HDplotTable' >> HDplots.R
echo '  HDplotTable$CHROM<-vcfData@fix[,"CHROM"]' >> HDplots.R
echo '  HDplotTable$POS<-vcfData@fix[,"POS"]' >> HDplots.R
echo '  HDplotTable$ID<-vcfData@fix[,"ID"]' >> HDplots.R
echo '  HDplotTable$depth_a<-A_reads' >> HDplots.R
echo '  HDplotTable$depth_b<-B_reads' >> HDplots.R
echo '  HDplotTable$ratio<-ratio' >> HDplots.R
echo '  HDplotTable$num_hets<-numHets' >> HDplots.R
echo '  HDplotTable$num_samples<-dim(hetMatrix)[2]' >> HDplots.R
echo '  HDplotTable$num_called<-numGenos' >> HDplots.R
echo '  HDplotTable$H_all<-hetPerc' >> HDplots.R
echo '  HDplotTable$H<-H' >> HDplots.R
echo '  HDplotTable$std<-std' >> HDplots.R
echo '  HDplotTable$D<-z' >> HDplots.R
echo '' >> HDplots.R
echo '  return(HDplotTable)' >> HDplots.R
echo '}' >> HDplots.R
echo '' >> HDplots.R
echo '#------------------#' >> HDplots.R
echo '# Loading the data #' >> HDplots.R
echo '#------------------#' >> HDplots.R
echo 'args <- commandArgs(trailingOnly = TRUE); names(args) <- c("file", "Hlim", "Dlim"); args' >> HDplots.R
echo 'vcf <- read.vcfR(as.character(args[1]), verbose = FALSE)' >> HDplots.R
echo '' >> HDplots.R
echo '#------------------------------------------------#' >> HDplots.R
echo '# Estimating parameters to assign paralog status #' >> HDplots.R
echo '#------------------------------------------------#' >> HDplots.R
echo 'res <- HDplot(vcfData = vcf)' >> HDplots.R
echo '# NaN --> no heterozygotes across samples (only 0/0 and 1/1 genotypes) or a std of the ratio is 0 (singleton heterozygotes or no variance in ratios).' >> HDplots.R
echo '# NA (for all estimates) --> multiallelic samples (this function only assesses biallelic loci).' >> HDplots.R
echo '' >> HDplots.R
echo '#---------------------#' >> HDplots.R
echo '# Determining filters #' >> HDplots.R
echo '#---------------------#' >> HDplots.R
echo 'pdf(file = "Distribution_H_D_Quantiles.pdf", width = 7, height = 5)' >> HDplots.R
echo 'hist(res$H, main = "", freq = F, xlab = "H"); abline(v = quantile(res$H, probs = c(0.05, 0.50, 0.95), na.rm = TRUE), col = "red", lty = 2, lwd = 1.5)' >> HDplots.R
echo 'hist(res$D, main = "", freq = F, xlab = "D"); abline(v = quantile(res$D, probs = c(0.05, 0.50, 0.95), na.rm = TRUE), col = "red", lty = 2, lwd = 1.5)' >> HDplots.R
echo 'dev.off()' >> HDplots.R
echo 'round(quantile(res$H, probs = c(0.05, 0.50, 0.95), na.rm = TRUE),3); mean(res$H, na.rm = TRUE)' >> HDplots.R
echo 'round(quantile(res$D, probs = c(0.05, 0.50, 0.95), na.rm = TRUE),3); mean(res$D, na.rm = TRUE)' >> HDplots.R
echo '' >> HDplots.R
echo 'Hlim <- as.numeric(args[2])' >> HDplots.R
echo 'Dlim <- as.numeric(args[3])' >> HDplots.R
echo '' >> HDplots.R
echo 'col <- sapply(1:nrow(res), function(i){' >> HDplots.R
echo '  if(is.na(res$H[i]) | is.na(res$D[i])) return(NA)' >> HDplots.R
echo '  if(res$H[i] > Hlim | abs(res$D[i]) > Dlim) return("red") else return("black")' >> HDplots.R
echo '})' >> HDplots.R
echo '' >> HDplots.R
echo '#-------------------------#' >> HDplots.R
echo '# Visualizing H against D #' >> HDplots.R
echo '#-------------------------#' >> HDplots.R
echo 'pdf(file="HD_plot.pdf", width = 7, height = 5)' >> HDplots.R
echo 'ggplot(data = res, aes(x = H, y = D)) + ' >> HDplots.R
echo '  geom_point(col = col, alpha = 0.6) + ' >> HDplots.R
echo '  geom_hline(yintercept = 0, col = "black", linetype = 2, linewidth = 0.5) +' >> HDplots.R
echo '  geom_vline(xintercept = Hlim, col = "red", linetype = 2, linewidth = 0.5) +' >> HDplots.R
echo '  geom_hline(yintercept = c(Dlim, -Dlim), col = "red", linetype = 2, linewidth = 0.5) +' >> HDplots.R
echo '  theme_minimal()+' >> HDplots.R
echo '  labs(x = "\nH\n", y = "\nD\n") +' >> HDplots.R
echo '  theme(' >> HDplots.R
echo '    axis.line = element_line(),' >> HDplots.R
echo '    axis.ticks = element_line(),' >> HDplots.R
echo '    axis.title = element_text(size = 14, face = "bold"),' >> HDplots.R
echo '    axis.text = element_text(size = 12)' >> HDplots.R
echo '  )' >> HDplots.R
echo 'dev.off()' >> HDplots.R
echo '' >> HDplots.R
echo '#---------------------------------------------------------#' >> HDplots.R
echo '# Creating a file recording site to keep based on filters #' >> HDplots.R
echo '#---------------------------------------------------------#' >> HDplots.R
echo 'nrow(res) # Total number of sites assessed.' >> HDplots.R
echo 'nrow(res[which(is.na(res$D) & res$H == 0),]) # No heterozygotes across samples (only 0/0 and 1/1 genotypes).' >> HDplots.R
echo 'nrow(res[which(is.na(res$D) & res$H != 0),]) # Missing depth for allele A and B (set as 0,0 for depth A,B).' >> HDplots.R
echo 'nrow(res[which(is.na(res$depth_a) & is.na(res$depth_b)),]) # Number of multiallelic sites across samples.' >> HDplots.R
echo 'keep <- res[-which(col == "red"),]' >> HDplots.R
echo 'range(keep$H, na.rm = TRUE) # Should range between 0 and Hlim.' >> HDplots.R
echo 'range(keep$D, na.rm = TRUE) # Should range between -Dlim and Dlim.' >> HDplots.R
echo 'keep <- keep[,c("CHROM", "POS")]' >> HDplots.R
echo 'write.table(keep, file = "keep_paralogs.txt", row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")' >> HDplots.R
echo '' >> HDplots.R
echo '#-----------#' >> HDplots.R
echo '# R Version #' >> HDplots.R
echo '#-----------#' >> HDplots.R
echo 'sessionInfo()' >> HDplots.R
echo 'sessioninfo::package_info()' >> HDplots.R

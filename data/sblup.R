library(Matrix)
options     <- commandArgs(trailingOnly = TRUE)
trait       <- options[1]
sumstatFile <- paste0(trait,".fastGWA")
hsq         <- 0.2 #as.numeric( options[2] )

gwas        <- read.table(sumstatFile,h=T,stringsAsFactors = F)
Msnps       <- nrow(gwas)
Nind        <- max(gwas[,"N"],na.rm=TRUE)
gwas$Z      <- (gwas$BETA/gwas$SE)
gwas$b      <- gwas$Z / sqrt(gwas$Z^2 + gwas$N)
rownames(gwas) <- gwas[,"SNP"]

system.time( ld <- read.table(paste0("/data/SISG2022M14/practicals/prac6/plinkLDMat_chrom20.ld.gz"),h=T,stringsAsFactors = F) )
snps <- unique(c(ld$SNP_A,ld$SNP_B))
snps <- intersect(snps,gwas[,"SNP"])

gwas <- gwas[snps,]
ld   <- ld[which(ld[,"SNP_A"]%in%snps & ld[,"SNP_B"]%in%snps),]
M    <- length(snps)

indexes <- 1:M
names(indexes) <- snps
i  <- c( indexes[ld[,"SNP_A"]], indexes[ld[,"SNP_B"]], indexes )
j  <- c( indexes[ld[,"SNP_B"]], indexes[ld[,"SNP_A"]], indexes )
r  <- c(ld[,"R"],ld[,"R"], rep(1,M))
R  <- sparseMatrix(i=i,j=j,x=r)
Im <- sparseMatrix(i=indexes,j=indexes,x=rep(1,M))

lambda  <- Msnps * (1-hsq) / ( Nind * hsq )
D       <- R + lambda * Im
bgwas   <- gwas[,"b"]
system.time( Bsblup <- solve(D,bgwas) )

SHET <- sqrt( 2*gwas[,"AF1"]*(1-gwas[,"AF1"]) )
results <- cbind.data.frame(gwas[,c("SNP","A1")],BSLUP=Bsblup[,1]/SHET)
write.table(results,paste0("SBLUP/",trait,".sblubInR.res"),quote=F,row.names=F,col.names=F)

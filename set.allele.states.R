####################################################################################################
##  Script Purpose: To update allele states according to reference dbSNP
##  
##  Author: M.Sc. Cainã Max Couto-Silva
##  Email: cmcouto.silva@gmail.com / cmcouto.silva@usp.br
##  Institution: University of Sâo Paulo, Brazil
##  
##  Date: Jan 09, 2018
##  Last modified: Jan 12, 2018
####################################################################################################

#### Section: Argument Parser Options
####################################################################################################

suppressPackageStartupMessages(if(!require(argparse)) install.packages("argparse"))
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-d", "--dataset", type="character", metavar="",
                    help="Path to dataset and its name (without extensions)")

parser$add_argument("-r","--ref-dbSNP", type="character", metavar="",
                    help="Path to SNP data base file (.vcf)")

parser$add_argument("-o","--out", type="character", metavar="",
                    help="Identification to output files (without extentions)")


args <- parser$parse_args()

ref_dbSNP <- args$ref_dbSNP
dataset <- args$dataset
out <- args$out

#### Section: Running Script
####################################################################################################

bimfile <- read.table(file = paste0(dataset,".bim"), header = F, sep = "\t", stringsAsFactors = F)
colnames(bimfile) <- c("Chromosome", "SNP_ID", "Genetic_dist", "Physical_pos", "Allele_1", "Allele_2")

writeLines(bimfile$SNP_ID, con = "bimfile.snps.txt")

message("\nSubsetting target-SNPs info (based on .bim file) from dbSNP.
This may take a while...")

system(paste("grep -Fwf bimfile.snps.txt", ref_dbSNP, paste0("> ", out, ".dbSNP.subset.vcf")))
message("Done\n")

message("Matching files according to their SNPs...")
ref_dbSNP <- read.table(file = paste0(out, ".dbSNP.subset.vcf"), header = F, sep = "\t", stringsAsFactors = F)
colnames(ref_dbSNP) <- c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO")

bimfile <- bimfile[bimfile$SNP_ID %in% ref_dbSNP$ID, ]
ref_dbSNP <- ref_dbSNP[order(match(ref_dbSNP$ID, bimfile$SNP_ID)), ]
row.names(ref_dbSNP) <- NULL

df <- data.frame(bimfile$SNP_ID, bimfile$Allele_1, bimfile$Allele_2, ref_dbSNP$REF, ref_dbSNP$ALT)
write.table(x = df, file = paste0(out, "_allele_update.txt"), quote = F, sep = "\t", row.names = F, col.names = F)

# Adjusting alleles using Plink Software

plink.version <- system("plink --version", intern = T)
message(paste("Adjusting alleles using", plink.version, "...\n"))
system(paste('plink --bfile', dataset, '--update-alleles', paste0(out,'_allele_update.txt'), '--make-bed --out', out))

message(paste("\nYou may want to remove", paste0(out, ".dbSNP.subset.vcf file.")))
unlink(c("bimfile.snps.txt", "my_output_allele_update.txt"))

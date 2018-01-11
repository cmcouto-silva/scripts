####################################################################################################
##  Script Purpose: Run Shapeit in check mode, performing all required tasks to run Shapeit 
##  officially
##  
##  Author: M.Sc. Cainã Max Couto-Silva
##  Email: cmcouto.silva@gmail.com / cmcouto.silva@usp.br
##  Institution: University of Sâo Paulo, Brazil
##  
##  Date: Jan 09, 2018
##  Last modified: Jan 10, 2018
####################################################################################################

#### Section: Argument Parser Options
####################################################################################################

suppressPackageStartupMessages(if(!require(argparse)) {install.packages("argparse") })
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-d", "--dataset", type="character", metavar="",
                    help="Path to dataset and its name (without extensions)")

parser$add_argument("-r","--ref", type="character", metavar="",
                    help="Path to genome reference panel")

parser$add_argument("-o","--out", type="character", metavar="",
                    help="Identification to output files (without extentions)")

parser$add_argument( "--plink-maf", default= 0.05, type="double", metavar="",
                     help="Numeric value for plink MAF parameter [default %(default)s]")

parser$add_argument( "--plink-mind", default= 0.1, type="double", metavar="",
                     help="Numeric value for plink mind parameter [default %(default)s]")

parser$add_argument( "--plink-geno", default= 0.1, type="double", metavar="",
                     help="Numeric value for plink geno parameter [default %(default)s]")

args <- parser$parse_args()


# Setting arguments to new variables

# Obrigatory
dataset_id <- args$out
raw_dataset_path <- args$dataset
ref_path <- args$ref
# Optional
plink_maf <- args$plink_maf
plink_mind <- args$plink_mind
plink_geno <- args$plink_geno

#### Section: Setting & Parameters
####################################################################################################

# Creating directories
cat("\n## Creating requested directories\n")

if (!dir.exists(dataset_id)){ dir.create(dataset_id, showWarnings = F) 
  cat(paste0(" Directory '", dataset_id, "' successfully created\n"))
} else {
  message(paste0("\n Directory", " \"", dataset_id, "\" already exists (it will be used)"))
}

directories= c('plink_filter', 'plink_split', 'shapeit_check01', 'shapeit_check02')
invisible(sapply(directories, function(mkdir){
  if (!dir.exists(paste0("./",dataset_id,"/", mkdir))){ dir.create(paste0("./",dataset_id,"/", mkdir), showWarnings = F)
    cat(paste0(" Directory \"", mkdir, "\" successfully created at ", dataset_id, "\n"))
  } else {
    message(paste0(" Directory '", mkdir, "' already exists (it will be used)"))
  }}))

rm(directories)
cat("\n")

# Setting paths
dataset_id_path <- paste0('./',dataset_id,'/')
plink_filter_path <- paste0(dataset_id_path, 'plink_filter/')
plink_split_path <- paste0(dataset_id_path, 'plink_split/')
shapeit_check01_path <- paste0(dataset_id_path, 'shapeit_check01/')
shapeit_check02_path <- paste0(dataset_id_path, 'shapeit_check02/')

# Plink filter
message(paste0("\n", "## Running plink filter parameters\n"))
system(paste("plink ", paste0("--bfile ",raw_dataset_path)," --make-bed --maf", 
             plink_maf," --mind ", plink_mind," --geno", plink_geno," --out ", 
             paste0(plink_filter_path, dataset_id)))

# Splitting chromosomes using Plink
message(paste0("\n", "## Splitting chromosomes\n"))
split_raw_dataset_path <- paste0(plink_filter_path, dataset_id)
split_output <- paste0(plink_split_path, dataset_id, '_chr')

for (i in 1:22) {
  system(paste0("plink --bfile ", split_raw_dataset_path, " --chr ", i, " --make-bed --out ", split_output,i ))
}

# Running Shapeit in checking mode

# Writting shell command for Shapeit as character vector
shapeit_check01 <- paste0("for i in `seq 1 22`; do shapeit -check \\
                          -B ", split_output,"$i \\
                          -M ",ref_path,"/genetic_map_chr$i\\_combined_b37.txt \\
                          --input-ref ",ref_path,"/1000GP_Phase3_chr$i.hap.gz ",ref_path,"/1000GP_Phase3_chr$i.legend.gz ",ref_path,"/1000GP_Phase3.sample \\
                          --output-log ", shapeit_check01_path, dataset_id, '_chr',"$i.aligments; done")

# Running Shapeit first check within terminal
message(paste0("\n", "## Running Shapeit's first check"))
system(shapeit_check01)

directories= c("strand","exclude","logs")
invisible(sapply(directories, function(mkdir){
  if (!dir.exists(paste0(shapeit_check01_path, mkdir))){ 
    dir.create(paste0(shapeit_check01_path, mkdir), showWarnings = F)
    cat(paste0("Directory \"", mkdir, "\" successfully created at ", shapeit_check01_path, "\n"))
  } else { 
    message(paste0(" Directory \"", mkdir, "\" already exists (it will be used)"))
  }
}))

rm(directories)
cat("\n")

check01_strand_path <- paste0(shapeit_check01_path,'strand/')
check01_exclude_path <- paste0(shapeit_check01_path,'exclude/')
check01_logs_path <- paste0(shapeit_check01_path,'logs/')

system(paste0("mv ", shapeit_check01_path,"*.strand ", check01_strand_path))
system(paste0("mv ", shapeit_check01_path,"*.exclude ", check01_exclude_path))
system(paste0("mv ", shapeit_check01_path,"*.log ", check01_logs_path))

# Inverting strands using Plink

message("\n## Strand inversion using plink")
shapeit_strand_chr <- paste0(check01_strand_path,dataset_id, '_chr')
strand_ext <- ".aligments.snp.strand"

for (i in 1:22) {
  file_name <- paste0(shapeit_strand_chr,i,strand_ext)
  data <- read.table(file_name, header = T, row.names = NULL, stringsAsFactors = F)
  rsID <- data[data$type=="Strand", 4]
  RefSnps_to_flip <- unique(rsID)
  write.table(RefSnps_to_flip, paste0(shapeit_strand_chr,i,"_strands_to_flip.txt"), 
              sep="\t", col.names=F, row.names=F, quote=F)
}

if (!dir.exists(paste0(dataset_id_path,"shapeit_plink_flipped_strands"))) { 
  dir.create(paste0(dataset_id_path,"shapeit_plink_flipped_strands"), showWarnings = F)
  cat(paste0(" Directory \"", 'shapeit_plink_flipped_strands', "\" successfully created at ", dataset_id_path, "\n"))
} else {
  message(paste0(" Directory \"", 'shapeit_plink_flipped_strands', "' already exists (it will be used)"))
}

cat("\n")
shapeit_plink_flipped_strands_path <- paste0(dataset_id_path,'shapeit_plink_flipped_strands/')

chrs_raw_dataset_path <- paste0(plink_split_path, dataset_id,'_chr')
chrs_output <- paste0(shapeit_plink_flipped_strands_path, dataset_id,'_chr')
chrs_flip.ref <- shapeit_strand_chr

for (i in 1:22) {
  system(paste0("plink --bfile ", chrs_raw_dataset_path,i, " --flip ", chrs_flip.ref,i,'_strands_to_flip.txt --make-bed --out ', chrs_output,i))
}

# Running Shapeit in checking mode

# Writing shell commands in R character vector
shapeit_check02 <- paste0("for i in `seq 1 22`; do shapeit -check \\
-B ", chrs_output,"$i \\
-M ",ref_path,"/genetic_map_chr$i\\_combined_b37.txt \\
--input-ref ", ref_path,"/1000GP_Phase3_chr$i.hap.gz ", ref_path,"/1000GP_Phase3_chr$i.legend.gz ", ref_path,"/1000GP_Phase3.sample \\
--output-log ", shapeit_check02_path, dataset_id, '_chr',"$i.aligments; done")

# Running Shapeit second check within terminal
message(paste0("\n", "## Running Shapeit's second check"))
system(shapeit_check02)

directories <- c("strand","exclude","logs")
invisible(sapply(directories, function(mkdir){
  if (!dir.exists(paste0(shapeit_check02_path, mkdir))){ 
    dir.create(paste0(shapeit_check02_path, mkdir), showWarnings = F) 
    cat(paste0(" Directory \"", mkdir, "\" successfully created at ", shapeit_check02_path, "\n"))
  } else { 
    message(paste0(" Directory \"", mkdir, "' already exists (it will be used)"))
  }
}))

rm(directories)
cat("\n")

check02_strand_path= paste0(shapeit_check02_path,'strand/')
check02_exclude_path= paste0(shapeit_check02_path,'exclude/')
check02_logs_path= paste0(shapeit_check02_path,'logs/')

system(paste0("mv ", shapeit_check02_path,"*.strand ", check02_strand_path))
system(paste0("mv ", shapeit_check02_path,"*.exclude ", check02_exclude_path))
system(paste0("mv ", shapeit_check02_path,"*.log ", check02_logs_path))


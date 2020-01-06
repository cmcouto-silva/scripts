####################################################################################################
##  Script Purpose: Run Shapeit in check mode, performing all required tasks to run Shapeit
##  officially
##
##  Author: M.Sc. Cainã Max Couto-Silva
##  Email: cmcouto.silva@gmail.com / cmcouto.silva@usp.br
##  Institution: University of Sâo Paulo, Brazil
##
##  Date: Jan 09, 2018
##  Last modified: Jan 23, 2019
####################################################################################################

#### Section: Argument Parser Options
####################################################################################################

suppressPackageStartupMessages(if(!require(argparse)) { install.packages("argparse") })
suppressPackageStartupMessages(library(argparse))

parser <- ArgumentParser()

parser$add_argument("-d", "--dataset", type="character", metavar="",
                    help="[REQUIRED] Path to raw dataset (without extensions)")

parser$add_argument("-r","--ref", type="character", metavar="",
                    help="[REQUIRED] Path to map and genome reference panel")

parser$add_argument("-o","--out", type="character", metavar="",
                    help="[REQUIRED] Name for output files (without extentions)")

parser$add_argument("--chrs", type="character", metavar="", default = "1:22",
                    help="Chromosomes IDs. It accepts IDs comma-delimited (e.g. 1,3,5,10) or colon (e.g. 1:10) as sequence delimiter. Default set from chromosome 1 to 22.")

parser$add_argument( "--plink-maf", type="double", metavar="",
                     help="Numeric value for plink MAF parameter")

parser$add_argument( "--plink-mind", type="double", metavar="",
                     help="Numeric value for plink mind parameter")

parser$add_argument( "--plink-geno", type="double", metavar="",
                     help="Numeric value for plink geno parameter")

args <- parser$parse_args()

## Setting arguments to new variables

# Obrigatory
dataset_id <- args$out
raw_dataset_path <- args$dataset
ref_path <- args$ref
chrs <- args$chrs

# Optional
if(!is.null(args$plink_maf)) {plink_maf <- args$plink_maf}
if(!is.null(args$plink_mind)) {plink_mind <- args$plink_mind}
if(!is.null(args$plink_geno)) {plink_geno <- args$plink_geno}

if(!is.null(chrs)) {
  
  if(grepl(",", chrs)) {
    chrs <- as.numeric(unlist(strsplit(chrs, ",")))
  } else {
    if(grepl(":", chrs)) {
      idx <- unlist( strsplit(chrs, ":") )
      chrs <- seq(idx[1], idx[2])  
    } else {
      chrs <- as.numeric(chrs)
    }
  }
  
}

#### Section: Setting & Parameters
####################################################################################################

# Creating directories
cat("\n## Creating requested directories\n")

if (!dir.exists(dataset_id)){ dir.create(dataset_id, showWarnings = F)
  cat(paste0(" Directory '", dataset_id, "' successfully created\n"))
} else {
  message(paste0("\n Directory", " \"", dataset_id, "\" already exists (it will be used)"))
}

directories <- c('plink_filter', 'plink_split', 'shapeit_check')
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
shapeit_check_path <- paste0(dataset_id_path, 'shapeit_check/')

if (grepl("/", dataset_id)) dataset_id <- gsub(".*/(.+)", "\\1", dataset_id)

# Plink filter
plink_v <- gsub(".*v(\\S..).*", "\\1", system("plink --version", T))
plink <- ifelse(plink_v >= 1.9, "plink", "plink --noweb")

message(paste0("\n", "## Running plink filter parameters\n"))
system(paste(plink, paste0("--bfile ",raw_dataset_path),"--make-bed", 
             if(exists("plink_maf")){paste("--maf", plink_maf)},
             if(exists("plink_mind")){paste("--mind", plink_mind)},
             if(exists("plink_geno")){paste("--geno", plink_geno)},
             "--allow-no-sex --out ", paste0(plink_filter_path, dataset_id)))

# Splitting chromosomes using Plink
message(paste0("\n", "## Splitting chromosomes\n"))
split_raw_dataset_path <- paste0(plink_filter_path, dataset_id)
split_output <- paste0(plink_split_path, dataset_id, '_chr')

for (i in chrs) {
  system(paste0(plink, " --bfile ", split_raw_dataset_path, " --chr ", i, " --allow-no-sex --make-bed --out ", split_output,i ))
}

# Running Shapeit in checking mode

# Writting shell command for Shapeit as character vector
shapeit_check <- paste0("for i in $(echo ", paste(chrs, collapse = " "), "); do shapeit -check \\
                          -B ", split_output,"$i \\
                          -M ", ref_path,"/genetic_map_chr$i\\_combined_b37.txt \\
                          --output-log ", shapeit_check_path, dataset_id, '_chr',"$i; done")

# Running Shapeit first check within terminal
message(paste0("\n", "## Running Shapeit's first check"))
system(shapeit_check)

directories <- c("snp.mm","ind.mm","logs")
invisible(sapply(directories, function(mkdir){
  if (!dir.exists(paste0(shapeit_check_path, mkdir))){
    dir.create(paste0(shapeit_check_path, mkdir), showWarnings = F)
    cat(paste0("Directory \"", mkdir, "\" successfully created at ", shapeit_check_path, "\n"))
  } else {
    message(paste0(" Directory \"", mkdir, "\" already exists (it will be used)"))
  }
}))

rm(directories)
cat("\n")

check_snp.mm_path <- paste0(shapeit_check_path,'snp.mm/')
check_ind.mm_path <- paste0(shapeit_check_path,'ind.mm/')
check_logs_path <- paste0(shapeit_check_path,'logs/')

system(paste0("mv ", shapeit_check_path,"*.snp.mm ", check_snp.mm_path))
system(paste0("mv ", shapeit_check_path,"*.ind.mm ", check_ind.mm_path))
system(paste0("mv ", shapeit_check_path,"*.log ", check_logs_path))

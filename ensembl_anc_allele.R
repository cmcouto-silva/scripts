ensembl_snp_annot <- function(snps, output, grch = 38) {

  ### List Ensembl Main Options
  # listEnsembl()

  # Set Ensembl main option
  ensembl <- biomaRt::useEnsembl("snp", GRCh = grch)

  ### List all data sets
  # listDatasets(ensembl)

  # Set SNP Dataset
  ensembl.dataset <- biomaRt::useDataset(dataset = "hsapiens_snp", ensembl)

  ### List Attributes for selected data set
  # listAttributes(ensembl.dataset)

  # Set SNP Attributes
  ensembl.attributes <- c('chr_name', 'refsnp_id', 'chrom_start', 'allele_1', 'allele')

  # Set SNP filter
  ensembl.filter <- snps

  # Run Annotation
  ensembl.annot <- biomaRt::getBM(attributes = ensembl.attributes, filters = 'snp_filter', values = ensembl.filter, mart = ensembl.dataset)

  # Write annotation file
  write.table(x = ensembl.annot, file = paste0(output, ".txt"), append = F, quote = F, sep = "\t", row.names = F, col.names = T)

}

## USAGE
# snps <- gt::read.bim(bim_file = "~/cmcouto.silva@usp.br/lab_files/all_datasets/Merged_data/merged.bim")[, SNP]
# output <- "~/cmcouto.silva@usp.br/lab_files/all_datasets/Annot/ensembl_grch37.txt"
#
# # Running function for GRCh 37 / hg19:
# ensembl_snp_annot(snps = snps, output = output, grch = 37)


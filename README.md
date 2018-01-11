# Genetic Tools

#### A set of programs for handling genetic data.

### Shapeit-check.R

This program performs all shapeit-check steps and associate steps to run properly [Shapeit software][shapeit.fr] .

It uses as input binary plink files (.bed .bim and .fam), and the [1000 Genomes][1000gen] as reference panel. Multiple folders pretty organized will be created, plotting respective files on each. Plink filter is applied accordingly user preferences, then plink files are splitted by chromossome, and the first check from Shapeit runs. Missaligned strands are recorded and inverted by plink, and second and last running of Shapeit is performed in order to get useful info files to officialy run Shapeit Phasing.

##### Requirements:
- Unix-based OS,
- [Shapeit Software][shapeit.fr],
- [Plink software (version => 1.9)][plink]
- [1000 Genomes' Reference Panel][1000gen.shapeit]

##### Run shapeit-check.R with the following command in Shell:
```
$ Rscript shapeit-check.R [options with arguments]
$ Rscript shapeit-check.R --help # help info
```

##### Here is an example of the --help session
```
usage: shapeit-check.R [-h] [-d] [-r] [-o] [--plink-maf] [--plink-mind] [--plink-geno]

optional arguments:
  -h, --help       show this help message and exit
  -d , --dataset   Path to dataset and its name (without extensions)
  -r , --ref       Path to genome reference panel
  -o , --out       Identification to output files (without extentions)
  --plink-maf      Numeric value for plink MAF parameter [default 0.05]
  --plink-mind     Numeric value for plink mind parameter [default 0.1]
  --plink-geno     Numeric value for plink geno parameter [default 0.1]
```

### set.allele.states.R

It uses as input binary plink files (.bed .bim and .fam), and a reference SNP data base in vcf format[link]. A temporary file with the subset of dbSNP which matches SNPs from .bim file is created, and these both files of same row numbers are equally organized line-by-line. This scripts uses Plink software for updating allele state information, generating new .bed .bim and .fam files named by user in -o (--out) parameter. See usage below.

##### Requirements:
- Unix-based OS,
- [Reference dbSNP][dbSNP],
- [Plink software (version => 1.9)][plink]

##### Run set.allele.states.R with the following command in Shell:
```
$ Rscript set.allele.states.R [options with arguments]
$ Rscript set.allele.states.R --help # help info
```

##### Here is an example of the --help session
```
usage: set.allele.states.R [-h] [-d] [-r] [-o]

optional arguments:
  -h, --help         show this help message and exit
  -d , --dataset     Path to dataset and its name (without extensions)
  -r , --ref-dbSNP   Path to SNP data base file (.vcf)
  -o , --out         Identification to output files (without extentions)
```

### LiftOver Conversion Script

Slightly modified script from [Patrick Deelen][liftover.script] for converting Hg18 (b36) to Hg19 (b37) genomes (positions).
A lot of very good resources/scripts can be found at his GitHub page.

##### Requirements:
- [LiftOver files,][liftover.files]
- [LiftOver Software][liftover.download],
- [Plink Software (version => 1.9)][plink]


[shapeit.fr]: <https://mathgen.stats.ox.ac.uk/genetics_software/shapeit/shapeit.html>
[plink]: <https://www.cog-genomics.org/plink2>

[dbSNP]: <ftp://ftp.ncbi.nlm.nih.gov/snp/organisms/human_9606_b150_GRCh37p13/VCF/>
[1000gen]: <http://www.internationalgenome.org/>
[1000gen.shapeit]: <https://mathgen.stats.ox.ac.uk/impute/1000GP_Phase3.html>

[liftover.script]: <https://github.com/molgenis/Imputation/issues/4>
[liftover.files]: <http://hgdownload.cse.ucsc.edu/downloads.html>
[liftover.download]: <http://hgdownload.soe.ucsc.edu/admin/exe/>

# SNPipeline [![Build Status](https://travis-ci.org/cbkmephisto/SNPipeline.svg?branch=master)](https://travis-ci.org/cbkmephisto/SNPipeline)

- Tool set for converting AB genotype calls from Illumina BeadExpress FinalFormat files or Affymetrix genotype files into a simpler "one row per animal" format, namely **ab-genotype** file.
- Tool set for converting **ab-genotype** files to appropriate format for genomic selection tool GenSel.
- Tool set for converting **ab-genotype** files to appropriate format for imputation/phasing tools Beagle v3/v4, FImpute, and findhap.
- Tool set for converting imputed or phased output files from Beagle v3/v4, FImpute, findhap back to **ab-genotype** files, or combining with map position info to create haplotype files.

- Incorporated into [BOLT software](http://www.thetasolutionsllc.com/bolt-software.html) package since [version 1.1.6](http://104.236.183.143/mediawiki/index.php/BOLT_Reference_Manual#SNPipeline)


## Compile
- Just cd into the folder and execute
```
make
```
or
```
make -j8
```
may use 8 cores to run in parallel

- Notice that ```SNPipeline.abg2bin``` depends on [eigen](http://eigen.tuxfamily.org/index.php), but other tools won't be affected if not installed.
  - for Ubuntu Linux, install ```libeigen3-dev``` before ```make```
  ```
  sudo apt-get update && sudo apt-get install libeigen3-dev
  ```
  - for Mac OS X, I used ```port```:
  ```
  sudo port install eigen3
  ```
  - ```Makefile``` may need to be modified for different systems.

- All executables will be created under the folder
```
bin/
```

## Executables
Run them without parameter to see syntaxes.


#### SNPipeline
> Creates **ab-genotype** file from finalReport/affyInput file(s).


#### SNPipeline.mergingAB
> Merges **ab-genotype** files to a single **ab-genotype** file.


### SNPipeline.trimDown2Map
> Trims the given **ab-genotype** file down on SNP names according to a 3-column (SNPName  chr  pos) map file.


### SNPipeline.pooledAB2VCF41
> Converts merged/unmerged **ab-genotype** files to VCF4.1 format, for imputation/phasing in beagle v4.


### SNPipeline.VCF2AB
> Converts phased/unphased, but must be imputed VCF4.1 file to merged/unmerged **ab-genotype** format.


### SNPipeline.abgExcl
> Pipes a **ab-genotype** file to a new one, excluding animalIDs in the blackLisk file.


### SNPipeline.abg2bin
> Converts the **ab-genotype** file to a -10/0/10 coded genotype file, replace missing data with column mean


### finalReportReformer
> To re-format ill-formatted finalReport files.


### gReplace
> To group-replace according to a given cross-reference list file.


### ggrep
> Group grep - greps lines from a target file acording to a given list file.


### abg2FImpute
> Proceed **ab-genotype** file to be a FImpute genotype file format. A **snp_info.txt** file is needed.


### mapUniter
> Unite 3-column ```SNP_ID chr pos``` map files to get a **snp_info.txt** genetic map file for FImpute.


### abg2findhap
> Converts **ab-genotype** files with pedigree info and genetic map info to be imputable by findhap V3.


### fhout2hapview
> Converts findhap output to HapView input to make LD plots.


### fhout2fiout
> Converts findhap output to the format of FImpute (because SNPipeline toolset contains other tools to convert from FImpute)


### fout2abg
> Converts FImpute output file to **ab-genotype** file.


### fout2haplotype
> Converts FImpute output file to haplotype file.


### abg2M
> Converts **ab-genotype** file to recoded (defined by a rule file) dense genotype matrix file.


### bout2genotype
> Reads beagle v3 output files to generate files for doing SNPRelate PCA (matrix: 0-1-2 & SNP-in-row genotype files).


### bout2haplotype
> Reads beagle v3 output files to generate haplotype files (matrix: 0-1 & SNP-in-col & window-grouped haplotype files).


### vcf2haplotype
> Reads beagle v4 phased output vcf format, converting to haplotype files.

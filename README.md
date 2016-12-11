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
will use 8 cores to run in parallel

- Notice that ```SNPipeline.abg2bin``` depends on [eigen](http://eigen.tuxfamily.org/index.php), but other tools won't be affected if not installed.
- All executables will be created under the folder
```
bin/
```

## Docs
### SNPipeline
```
####################################################################
#######
###                 SNPipeline
####
############
##                              by hailins@iastate.edu
#####                           created 2013-10-01
####                            updated 2015-08-12
#######
####################################################################

    SNPipeline - step 1. Creates AB-genotype file from finalReport/affyInput file(s). An optional space/tab-delimited >=2-col file named 'xref' is considered: [0]->[last].

[Warnings]
    - ALWAYS MAKE THE INPUT FILES INTO DIFFERENT FOLDERS!! If 2 or more finalReport/affyInput files with the same [CONTENT] and numberSNPs exist in the SAME DIRectory, the resulting ab-genotype file for these files would have the same fileName, and only the last processed file will be in the folder-related block of the resulting ab-genotype file. It's BLOCK-REPLACING, neither APPENDING nor MERGING!

[xref]
    - A space/tab-delimited->=2column-per-line sampleID cross-reference file in the same directory with finalReport/affyInput file named 'xref' will be used if exists.
    - This program will replace all the spaces in the sampleIDs in the input finalReport/affyInput files by underscore char('_'), thus no white space char(' ') is expected in either column of the IDs in the xref file. Actually ' ' is a delimit char, while tab('	', '\t') is the other.

[Progress % Legend]
    C for  pre-checking
    P for post-checking
    L for loading genotype
    T for creating temp file
    W for writing


usage:
        SNPipeline 'Red Angus 50K 29feb2012_FinalReport.txt' (*.txt)
```
### SNPipeline.mergingAB
```
####################################################################
#######
###                 SNPipeline.mergingAB
####
############
##                              by hailins@iastate.edu
#####                           created 2013-11-22
####                            updated 2014-08-07
#######
####################################################################

    SNPipeline - setp 2. Creates merged-AB-genotype-file from AB-genotype file(s).

usage:
        1) SNPipeline.mergingAB -x SNPxref -o outputFile -l list_file [-v] [-r 0.95]
        2) SNPipeline.mergingAB -x SNPxref -o outputFile -i AB-genotypeFile1 [AB-genotypeFile2 *] [-v] [-r 0.95]
        3) SNPipeline.mergingAB -x SNPxref -o outputFile -i AB-genotypeFile1 [AB-genotypeFile2 *] -l list_file [-v] [-r 0.95]

options:
        -x SNPxref file name
        -o output file name
        -i input file name (list)
        -l input list file containing input files, seperated by new line
        -v verbose output
        -r 0~1, calling rate limitation: only write out samples with calling rate greater than/o/equal to the given value.
           Not specifying this value is the same as specifying 0.

The order of [-r], [-v], -x block, -i block, -l block and -o block doesn't matter.
```
### SNPipeline.trimDown2Map
```
####################################################################
#######
###                 TrimDown2Map
####
############
##                              by hailins@iastate.edu
#####                           created 2014-01-07
####                            updated N/A
#######
####################################################################

    TrimDown2Map - trim the given ab-genotype file down on SNP names given in a 3-column (SNPName  chr  pos) map file.

usage:
        SNPipeline.trimDown2Map mapFile genotypeFile
```
### SNPipeline.pooledAB2VCF41
```
####################################################################
#######
###                 SNPipeline.pooledAB2VCF41
####
############
##                              by hailins@iastate.edu
#####                           created 2014-02-22
####                            updated N/A
#######
####################################################################

    SNPipeline - setp 4. Convert merged/unmerged AB genotype files to VCF4.1 format.

usage:
SNPipeline.pooledAB2VCF41 mapFile AB-GenotypeFile

mapFile
  3-column file, SNPName chr pos, whitespace/tab delimited
  NO HEADER LINE, or headerline starts with '#'. ANY '#' line would be ignored while making the map.

AB-GenotypeFile
  merged/unmerged ab-genotype file, 1st line containing 'SNP_IDs' and SNPNames delimited by #. Any other lines starting with # would be discarded.

Discarding SNPs if
  - SNP not in map file
  - SNP genotype is all missing in the current file
  - SNPs sharing the same position, only 1 would be kept
```
### SNPipeline.VCF2AB
```
####################################################################
#######
###                 SNPipeline.VCF2AB
####
############
##                              by hailins@iastate.edu
#####                           created 2014-02-26
####                            updated N/A
#######
####################################################################

    SNPipeline - setp ?. Convert phased/unphased, but must be imputed VCF4.1 file to merged/unmerged AB genotype format.

usage:
    SNPipeline.VCF2AB VCF4.1_File
```
### SNPipeline.abgExcl
```
SNPipeline.abgExcl    Version<0.01>    Mar 2014    hailins@iastate.edu

    Pipe the ab-genotype file to a new one excluding animalIDs in the blackLisk file.

Syntax:
    SNPipeline.abgExcl <blackList_file_name> <ab-genotype_file_name> -o <abg-output_file_name>
```
### SNPipeline.abg2bin
```
abg2bin    Version<0.02>    Apr 2015    hailins@iastate.edu

Syntax: SNPipeline.abg2bin [-a2g | -g2b] genotype_file_name
  * option -a2g converts the ab-genotype file to a -10/0/10 coded genotype file, replace missing data with column mean
    - non-segregating SNP locus will be coded as -77
  * option -g2b converts the -10/0/10 genotype file to a binary file
  * if neither -a2g nor -g2b was specified, abg2bin will convert ab-genotype file to binary format coding missing value as "1" in a -10/0/10 coding system.
```
### finalReportReformer
```
####################################################################
#######
###                 finalReportReformer
####
############
##                              by hailins@iastate.edu
#####                           created 2013-10-30
####                            updated N/A
#######
####################################################################

    Re-format ill-formatted finalReport files.

usage:
        finalReportReformer -h header -r SNP_Ref -o outputFile -i inputFile

options:
        -r SNP_Ref file name: the expected output SNP_Name template, one SNP name per line, output with the same order
        -o output  file name
        -i input   file name: full ill-formatted finalReport file
        -h header  file name: the expected/corrected header part, including [Header] section, [Data] line and anotation line
```
### gReplace
```
gReplace: Group Replace by hailins@iastate.edu Thr Apr 10 10:58 AM 2014

Usage:

gReplace 2col_xref_file target_file > newFile

This program
 - replaces the words in the col1 of 2col_xref_file to col2 if the col1 words appeared in the target_file
 - print out on screen, delimetered by whitespace.
```
### ggrep
```
ggrep: Group Grep by hailins@iastate.edu Wed Mar 26 2:18PM 2014

Usage:

ggrep lst_file target_file > newFile
ggrep -v blacklst_file target_file  > newFile

This program
 - greps the lines acording to that 1st col of target_file listed in lst_file, OR
 - excludes the lines of target_file for all the 1st col listed in blacklst_file (-v).
 - lines started with # will always be grepped and printed out.
 ```
### abg2FImpute
```
####################################################################
#######
###                 abg2FImpute
####
############
##                              by hailins@iastate.edu
#####                           created 2014-03-27
####                            updated 2014-04-23
#######
####################################################################

    Proceed ab-genotype file to be a FImpute format. A snp_info.txt file is needed.

 [+] 20:57:04 | Starting 2016-12-10, with PID = 5134

usage:
        abg2FImpute snp_info.txt genotype1 (*)
```
### mapUniter
```
####################################################################
#######
###                 mapUniter
####
############
##                              by hailins@iastate.edu
#####                           created 2014-03-24
####                            updated N/A
#######
####################################################################

    Unite map files to get a snp_info file for FImpute.

 [+] 20:57:46 | Starting 2016-12-10, with PID = 5138

usage:
        mapUniter map1 map2 (*)
```
### abg2findhap
```
####################################################################
#######
###                 abg2findhap
####
############
##                              by hailins@iastate.edu
#####                           created 2014-11-10
####                            updated 2014-11-10
#######
####################################################################

    Convert ab-genotype files with pedigree and map info to be imputable by findhap V3.

 [+] 20:58:00 | Starting 2016-12-10, with PID = 5140
 [-] 20:58:00 X ERROR: No ped file    (-p) specified.
 [-] 20:58:00 X ERROR: No map file(s) (-m) specified.
 [-] 20:58:00 X ERROR: No gen file(s) (-g) specified.

        ======================================================================
        ========================        ReadMe        ========================
        ======================================================================

        * usage
        abg2findhap [-s] -p PEDIGREE -m MAP_FILE(s) -g GEN_FILE(s)

        * THE OPTIONAL -s OPTION
        - generate output files by chromosome
        - is off by default

        * THE PEDIGREE
        - reads pedRefiner output .csv format(coma-separated-values) 3-col pedigree file, convert to findhap V3 format
        and make xref for IDs
        - no checking ability was implemented, so
        - pedRefiner version >= 2014-11-10 is required so that the output is checked and sorted
        - date of birth will be made up
        - unknown gender will be converted to M

        * THE MAP(s)
        - reads 'SNP_NAME CHR POS' 3-col-through, space/tab delimetered file
        - skip 1st line as header

        * THE GEN_FILE(s)
        - ab-genotype file from SNPipeline
```
### fhout2hapview
```
####################################################################
#######
###                 fhout2hapview
####
############
##                              by hailins@iastate.edu
#####                           created 2015-08-10
####                            updated N/A
#######
####################################################################

    Converts findhap output to HapView input to plot LD...

 [+] 21:06:24 | Starting 2016-12-10, with PID = 5207
 [-] 21:06:24 X ERROR: No ped file (-p) specified.
 [-] 21:06:24 X ERROR: No map file (-m) specified.
 [-] 21:06:24 X ERROR: No gen file (-g) specified.
 [-] 21:06:24 X ERROR: No out file (-o) specified.
 [-] 21:06:24 X ERROR: target chr (-c) should be specified.
 [-] 21:06:24 X ERROR: starting win (-w) should be specified.
 [-] 21:06:24 X ERROR: number of windows (-n) should be specified.

usage:
        fhout2hapview [-i inc] -c 14 -w 24 -n 1 -p pedigree.file -m chromosome.data -g haplotypes.txt -o PREFIX
```
### fhout2fiout
```
####################################################################
#######
###                 fhout2fiout
####
############
##                              by hailins@iastate.edu
#####                           created 2014-11-16
####                            updated 2014-11-16
#######
####################################################################

    Converts output of findhap to the format of FImpute (because SNPipeline toolset contains other tools to convert from FImpute)

 [+] 20:58:32 | Starting 2016-12-10, with PID = 5143
 [-] 20:58:32 X ERROR: No ped file (-p) specified.
 [-] 20:58:32 X ERROR: No map file (-m) specified.
 [-] 20:58:32 X ERROR: No gen file (-g) specified.
 [-] 20:58:32 X ERROR: No out file (-o) specified.

usage:
        fhout2fiout -p pedigree.file -m chromosome.data -g haplotypes.txt -o PREFIX
```

### fout2abg
```
####################################################################
#######
###                 fout2abg
####
############
##                              by hailins@iastate.edu
#####                           created 2014-03-31
####                            updated 2014-04-18
#######
####################################################################

    Converts FImpute output file to ab-genotype file.

 [+] 21:07:06 | Starting 2016-12-10, with PID = 5211

usage:
        fout2abg snp_info.txt genotypes_imp.txt -o ab-genotype.output.abg
```
### fout2haplotype
```
####################################################################
#######
###                 fout2haplotype
####
############
##                              by hailins@iastate.edu
#####                           created 2014-04-21
####                            updated 2016-04-25
#######
####################################################################

    Converts FImpute output file to haplotype file.
	+/- in the output SNP loci means (heterozygous + unable to phase)
	 X  in the output SNP loci means unable to impute. Will read a file windowSize for window size if exists. Default 1M.

 [+] 21:00:12 | Starting 2016-12-10, with PID = 5148

usage:
    fout2haplotype snp_info.txt genotypes_imp.txt -o haplotype.from.fout [ chr win ]
```
### abg2M
```
abg2M    Version<0.01>    Dec 2015    hailins@iastate.edu

Syntax: abg2M <rule_file> <ab-genotype_file_name> <output_prefix>

 - <rule_file> is the file defining conversion. For a -10, 0, 10 coding system, use sample below:
**** file content below
AA -10
AB 0
BB 10
**** file content above
 - <output_prefix> defines the prefix name of the output files. 3 output files will be generated:
 - <output_prefix>.anm_ID, storing animal IDs, one per line
 - <output_prefix>.SNP_ID, storing SNP IDs, one per line
 - <output_prefix>.mtx_ds, storing the dense coded matrix
```
### bout2genotype
```
####################################################################
#######
###                 bout2genotype
####
############
##                              by hailins@iastate.edu
#####                           created 2013-11-04
####                            updated 2014-06-30
#######
####################################################################

    read beagle v3 output files to generate files for doing SNPRelate PCA.
	matrix: 0-1-2 & SNP-in-row genotype files

 [+] 21:08:05 | Starting 2016-12-10, with PID = 5219

usage:
        bout2genotype incFile mapFile beagle_(*.phased)
```
### bout2haplotype
```
####################################################################
#######
###                 bout2haplotype
####
############
##                              by hailins@iastate.edu
#####                           created 2013-11-19 10:59 CDT
####                            updated 2015-06-03 16:59 CDT
#######
####################################################################

    read beagle output files to generate haplotype files.
	matrix: 0-1 & SNP-in-col & window-grouped haplotype files

 [+] 21:10:28 | Starting 2016-12-10, with PID = 5221

usage:
        bout2haplotype 1.0 incFile mapFile beagle_(*.phased)
```
### vcf2haplotype
```
####################################################################
#######
###                 vcf2haplotype
####
############
##                              by hailins@iastate.edu
#####                           created 2014-03-13
####                            updated N/A
#######
####################################################################

    Read beagle4 phased output vcf format, converting to haplotype info.

 [+] 21:11:05 | Starting 2016-12-10, with PID = 5224

usage:
        vcf2haplotype imputed/phased.*.vcf
```

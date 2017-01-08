# Tutorial for SNPipeline

## Files in the example
In a \*nix system terminal, the command `cat FILENAME` would display file content of `FILENAME`. For example,
```
Wed Jan 04 12:05:30 cbkmephisto@cervine $ cat example.map
#SNP_IDs chr pos
ARS-BFGL-BAC-10172 1 1
ARS-BFGL-BAC-1020 2 2
ARS-BFGL-BAC-10245 3 3
ARS-BFGL-BAC-10345 4 4
ARS-BFGL-BAC-10365 4 7
ARS-BFGL-BAC-10591 5 8
ARS-BFGL-BAC-10793 6 9
```
Below are the descriptions for files in the example.
```
- example.SNPipeline.finalReport
  Illumina finalReport file, as an input file for SNPipeline.

- example.SNPipeline.xref
  A simple two or more columns, whitespace-delimitered, first column to last column
  crossreferencing file. Optional for SNPipeline. If it is required to correct animal ID(s)
  while converting finalReport file to ab-genotype file, make a file named "xref" in the
  same folder with finalReport file. In this example folder, it is renamed to make this
  option disabled.

- example.map
  A "SNP_ID chromosome bp_position" genetic map file. Used by some tools.
```

## Converting *finalReport* to *ab-genotype*
Assuming all the executable binary files were put in $PATH so that we don't need to specify the fullpath of SNPipeline
```Bash
cp example.SNPipeline.finalReport testFR # make a copy with simpler file name
SNPipeline testFR # run SNPipeline
```
Outputs on screen should look like this
```
+ 12:11:52 | { testFR }
+ 12:11:52 | Duplicated sampleID(s) found:                           
           10873027
           108373027
+ 12:11:52 | Try continuing with duplicated sampleIDs

X 12:11:52 | WARNING: duplicated sampleID(s) found after xref:
           10873027
           108373027
+ 12:11:52 |   - AB-genotype file [ ab-genotyped-7-test.bpm ] not found, proceed to create A NNNNNEEEEEWWWWW ONE.
+ 12:11:52 |   - Finally, 5 new sample(s) WRITEN into [ ab-genotyped-7-test.bpm ] within 0 seconds. √
```

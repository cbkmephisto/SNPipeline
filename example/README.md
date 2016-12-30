# Tutorial for SNPipeline

## Files in the example
```
- example.SNPipeline.finalReport
  Illumina finalReport file, as an input file for SNPipeline.

- example.SNPipeline.xref
  A simple two or more columns, whitespace-delimitered, first column to last column crossreferencing file. Optional for SNPipeline. If it is required to correct animal ID(s) while converting finalReport file to ab-genotype file, make a file named "xref" in the same folder with finalReport file. In this example, it is renamed to make this option disabled.

- example.map
  A "SNP_ID chromosome bp_position" genetic map file. Used by some tools.
```

## Converting *finalReport* to *ab-genotype*
Assuming all the executable binary files were put in $PATH
```Bash
cp example.SNPipeline.finalReport testFR
SNPipeline testFR
```

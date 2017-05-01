#OVERVIEW

Kourami is a graph-guided assembler for HLA haplotypes covering typing exons (exons 2 and 3 for Class I and exon 3 for Class II) 
using high-coverage whole genome sequencing data. Kourami constructs highly accurate haplotype sequences at 1-bp resolution by 
encoding currently available HLA allelic sequences from IPD-IMGT/HLA Database ( http://www.ebi.ac.uk/ipd/imgt/hla/ ) 
as a partial-ordered graph. Each database allele is naturally encoded as a path through the graph and any detectable genetic 
variations (SNPs or indels) not captured by the known sequences are added to the graph by graph-modification based on read alignment. 
Unlike previously available WGS-based HLA typing methods (database-matching techniques), Kourami direclty assembles the sequence.


# Release

The latest release, including both jar and source code can be downloaded from [here](https://github.com/Kingsford-Group/kourami/releases/tag/v0.9.0).


# Installation

To install Kourami, you must have following installed on your system:

- JDK 1.7+ 

- Apache Maven (3.3+) or Apache Ant (1.9+) is required (we recommend Maven for each dependency downloads)
  - To use Ant, you must have dependencies downloaded.

To compile and generate a jar file run the following command from the kourami directory where pom.xml is located.
```
jvn install
```

This will create a "target" directory and place a packaged jar file in it.

# Usage
```
java -jar <PATH_TO>/Kourami.jar [options] <bam-1> ... <bam-n>
```

Option Tag | Description
------------ | ----------
-d,--msaDirectory <path> | build HLAGraph from gen and nuc MSAs provided by IMGT/HLA DB from given directory (required)
-h,--help | print this message
-o,--outfilePrefix <outfile> | use given outfile prefix for all output files (required)

#Output

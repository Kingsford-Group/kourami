# Overview

Kourami is a graph-guided assembler for HLA haplotypes covering typing exons (exons 2 and 3 for Class I and exon 3 for Class II) 
using high-coverage whole genome sequencing data. Kourami constructs highly accurate haplotype sequences at 1-bp resolution by 
first encoding currently available HLA allelic sequences from IPD-IMGT/HLA Database ( http://www.ebi.ac.uk/ipd/imgt/hla/ ) 
as partial-ordered graphs. Each database allele is naturally encoded as a path through the graph and any detectable genetic 
variations (SNPs or indels) not captured by the known sequences are added to the graph by graph-modification based on read alignment 
to capture differences novel alleles have compared to known sequences. Unlike previously available WGS-based HLA typing methods 
(database-matching techniques), Kourami direclty assembles both haplotypes for each HLA gene (HLA-A, -B, -C, -DQA1, -DQB1, -DRB1). 
It also provides the typing result (6-digit 'G' resolution) by outputing the best matching alleles among the known sequences.


# Release

The latest release, including both jar and source code can be downloaded from [here](https://github.com/Kingsford-Group/kourami/releases/tag/v0.9).


# Installation

To install Kourami, you must have following installed on your system:

- JDK 1.7+ 

- Apache Maven (3.3+) or Apache Ant (1.9+) is required (we recommend Maven for each dependency downloads)
  - OR you must have dependencies downloaded and added to your CLASSPATH and compile using javac.
  - To use Ant, you must have dependencies downloaded.

[MAVEN USERS] To compile and generate a jar file run the following command from the kourami directory where pom.xml is located.
```
mvn install
```

[ANT USERS] To compile and generate a jar file run the following command from the kourami directory where build.xml is located.
```
ant compile jar
```

This will create a "target" directory and place a packaged jar file in it.

# Usage
```
java -jar <PATH_TO>/Kourami.jar [options] <bam-1> ... <bam-n>
```

Option Tag | Description
----------------------- | -----------------------------
-h,--help | print this message
-d,--msaDirectory <path> | build HLAGraph from gen and nuc MSAs provided by IMGT/HLA DB from given directory (required)
-o,--outfilePrefix <outfile> | use given outfile prefix for all output files (required)

# Output


#Dependencies
Dependecies can be easily downloaded by using Maven install command.

In each release, the pre-compiled jar is distributed with all necessary jars for dependencies, and they are:

- JGraphT 0.9.1 ( http://jgrapht.org/ )
- Apache Commons CLI 1.4 ( https://commons.apache.org/proper/commons-cli/ )
- fastutil 7.0.13 : Fast & compact type-specific collections for Java ( http://fastutil.di.unimi.it/ )

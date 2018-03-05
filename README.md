<pre>
-hhy+.                o o       o o       o o o o       o o
.`           -syss:---.`        o     o o o     o o o         o o o     o o o
:+:`     .:/o+++++///ommy+`    o       _  __                               _
`yhs/..:osssooooo++++dmNNNdo`   o     | |/ /___  _   _ _ __ __ _ _ __ ___ (_)
 /syy///++++ooooooooodNMdNdmh: o      | ' // _ \| | | | '__/ _` | '_ ` _ \| |
 -do/` .://++++++++oodmmmmmmd-        | . \ (_) | |_| | | | (_| | | | | | | |
 .+:     `.://///+///ommmmdy-         |_|\_\___/ \__,_|_|  \__,_|_| |_| |_|_|
  .          -syo----..``          
            +y+.                
</pre>

# Overview

Kourami is a graph-guided assembler for HLA haplotypes covering typing exons (exons 2 and 3 for Class I and exon 3 for Class II) 
using high-coverage whole genome sequencing data. Kourami constructs highly accurate haplotype sequences at 1-bp resolution by 
first encoding currently available HLA allelic sequences from IPD-IMGT/HLA Database ( http://www.ebi.ac.uk/ipd/imgt/hla/ ) 
as partial-ordered graphs. Each database allele is naturally encoded as a path through the graph and any detectable genetic 
variations (SNPs or indels) not captured by the known sequences are added to the graph by graph-modification based on read alignment 
to capture differences novel alleles have compared to known sequences. Unlike previously available WGS-based HLA typing methods 
(database-matching techniques), Kourami direclty assembles both haplotypes for each HLA gene (HLA-A, -B, -C, -DQA1, -DQB1, -DRB1).
From version 0.9.4 or later, Kourami supports additional HLA loci.  It also provides the typing result (6-digit 'G' resolution) by
outputing the best matching alleles among the known sequences whenever 'G' grouping information is available.


# Release

The latest release, including both jar and source code can be downloaded from [here](https://github.com/Kingsford-Group/kourami/releases/latest).



# Installation

To install Kourami, you must have following installed on your system:

- JDK 1.8+

- Apache Maven (3.3+) or Apache Ant (1.9+) is required (we **recommend Maven** for easy dependency downloads)
  - OR you must have dependencies downloaded and added to your CLASSPATH. Then you can compile using javac.
  - To use Ant, you must have dependencies downloaded and place jars under 'exjars' directory. 'exjars' directory must be created.

-A copy of the preformatted IMGT-HLA database (Kourami panel) can be obtained using a script. The panel sequence file needs to be bwa indexed before using and this is NOW done by the script when it downloads the database. The script will download and install the database under ```db``` directory under the Kourami installation directory. The download and index script can be run from the kourami installation directory:

```
scripts/download_panel.sh
```

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
NOTE: kourami jar takes a **bam aligned to Kourami reference panel built from IMGT/HLA db** (included in the preformatted IMGT-HLA database). 
Detailed notes on how to generate input bam consisting of HLA loci reads aligned to known alleles is explained in [How to prepare input bam and HLA panel for Kourami](https://github.com/Kingsford-Group/kourami/blob/master/preprocessing.md).

Option Tag | Description
----------------------- | -----------------------------
-h,--help | print this message
-d,--msaDirectory \<path> | build HLAGraph from gen and nuc MSAs provided by IMGT/HLA DB from given directory (required). Can be downloaded by running ```scripts/download_panel.sh```.
-o,--outfilePrefix \<outfile> | use given outfile prefix for all output files (required)
-a,--additionalLoci           | type additional loci (optional)
# Output

\<outfileprefix>.result contains the typing result and the columns are:  
1: Allele  
2: #BasesMatched  
3: Identity (#BasesMatched/MaxLen(query, db_allele))  
4: Length of the assembled allele  
5: Length of the matched allele from IMGT/HLA DB  
6: Combined bottleneck weights of both paths at a position. This is not necessarily same as the sum of column 7 and 8.
7: Weight of the bottleneck edge in path 1
8: Weight of the bottleneck edge in path 2

Note: Given a path, a bottleneck edge is an edge with the minimal weight. For an allele, there are always two entries (lines) reported in the result file. Path 1 is reported first, and path 2 is reported in the following line. The columns 6 7 8 are going to be redundant (same) for both lines. 
   
\<outfileprefix> contiains program log  
  
Assembled allele sequences are outputed in files ending with .typed.fa (multi-FASTA format)


# Dependencies
Dependecies can be easily downloaded by using Maven install command.

In each release, the pre-compiled jar is distributed with all necessary jars for dependencies, and they are:

- JGraphT 0.9.1 ( http://jgrapht.org/ )
- Apache Commons CLI 1.4 ( https://commons.apache.org/proper/commons-cli/ )
- fastutil 7.0.13 : Fast & compact type-specific collections for Java ( http://fastutil.di.unimi.it/ )

# How to cite Kourami
Please cite our [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-018-1388-2) available on Genome Biology: 

Lee, H., & Kingsford, C. Kourami: graph-guided assembly for novel human leukocyte antigen allele discovery. *Genome Biology* 19(16), 2018

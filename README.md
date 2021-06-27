# DinucISG
Dinucleotide calculations for Shaw et al ISG Zap paper

The ISG_CG.java program is a simple script to calculate a range of nucleotide and dinucleotide features from a FASTA sequence file. It expects sequences to be on a single line. If the sequence name is in the pipe delimited form >Gene|Transcript it will extract the Gene and Transcript names for outputting in the output file, otherwise if no pipe is detected, both the gene and transcript name will simply equal the whole sequence name. Usage is:

```
java -jar ISG_CG.jar sequences.fasta
```

This will create a file called sequences_dat.txt with the sequence name and nucleotide/dinucleotide measures (both raw percentage and observed/expected). The DinucData folder contains the results of running ISG_CG.jar across a range of species cDNA and CDS sequences as downloaded from [Ensembl](https://www.ensembl.org/index.html) release 92. The two mouse data files (mouse_cds.. and mouse_cdna...) are additional files with just 100 genes analysed and are included here for completeness. The V6 subfolder contains earlier versions of the same data files, these files lack the AA% (column 31) to UT% (column 44); in addition, in the V6 files the column UP% represents UA% - the percentage of all dinucleotides that are UA (TA in DNA) in the sequence.

The CPB_Machine_ISG.jar program is a slightly stripped out version of the CPB_Machine.jar program available from the [VirusFeatrues repo](https://github.com/rjorton/VirusFeatures) to run more easily on non-virus sequences (taking the whole sequence name as the name, not required an virus accession lookup file, etc). This file produces many more features than ISG_CG.jar including dinucleotides at the bridge and non-bridge positions and codon-pair-bias scores for the sequence. IMPORTANT - this program should only be run on coding sequences, as it is expected the sequence to be in frame. Usage is:

```
java -jar CPB_Machine_ISG.jar seqeunces.fasta
```

This will create a file called sequences_cpb_dat.txt with all the sequence features in it.

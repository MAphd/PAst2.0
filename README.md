# PAST2.0 
Pseudomonas aeruginosa serotyper 2

Blast based tool which enables in silico serotyping of Pseudomonas aeruginosa genomes. Successor to PAst [1] with following modifications:
* Included O15 and O17 clusters in database for typing.
* Calculates coverage better for fragmented query genomes 

Call from commandline using 

```
Rscript PAST.R 
```

Avaliable commands:

```
-i Input directory path (Must be a file compatible with blastn with file extensions: .fasta, .fsa, .fas, .fna, .fa. or .seq). Defaults to ./input/ if not specified
-o Output directory path (defaults to ./output/ if not specified)
-b (If you already previously performed the blast alignment for debugging purposes, defaults to FALSE) 
-p (To enable installation of missing packages required to run script)
-h or -help Display this message
```

Required software:
* Rscript
* blastn 

Required R packages:
* readr

Has been tested using Rscript 3.6.3 for Ubuntu 20.04 with blastn 2.9.0

Citation Anbo M, Jelsbak L. 2023. A bittersweet fate: detection of serotype switching in Pseudomonas aeruginosa. Microb Genomics 9:000919. https://doi.org/10.1099/mgen.0.000919


```
[1] Thrane SW, Taylor VL, Lund O, Lam JS, Jelsbak L. Application of Whole-Genome Sequencing Data for O-Specific Antigen Analysis and In Silico Serotyping of Pseudomonas aeruginosa Isolates. J Clin Microbiol. 2016 Jul;54(7):1782-1788. doi: 10.1128/JCM.00349-16. Epub 2016 Apr 20. PMID: 27098958; PMCID: PMC4922119.
```


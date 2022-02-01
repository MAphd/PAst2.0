
library(readr)

#Load required functions

Usedir <- dirname(rstudioapi::getSourceEditorContext()$path)
#Attempt to define directory automatically



source(paste0(Usedir,"/Functions/Past2.0.R"))
source(paste0(Usedir,"/Functions/WeightedAlignmentRouter.R"))

#Define variables for running PAst:
# output: output directory
# input: input directory
# osadb: osadb.fasta directory  
# blasted: TRUE/FALSE - omits blasting step if set to TRUE

#Define these if not using original file structure
Inputdir <- paste0(Usedir,"/Input/")
Outputdir <- paste0(Usedir,"/Output/")
OsaDBdir <- paste0(Usedir,"/Functions/OSAdb.fasta")

PAst(Outputdir,Inputdir, OsaDBdir, FALSE)

#Output file will be saved to Outputdir with the suffix "_Serotyping.txt".
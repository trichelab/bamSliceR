library(BSgenomeForge)

# source of sequence: gencode.v36.transcripts.fa.gz 
# from https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_36/gencode.v36.transcripts.fa.gz
# trim header to "ENSTXXXXXX.X"
# download faToTwoBit utility from https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/
# convert fasta to 2bit format.

library(BSgenomeForge)
setwd("/varidata/research/projects/triche/Peter/BamSlicing/BSgenome.gencode.v36.transcripts")
forgeBSgenomeDataPkg(x = "./BSgenome.Hsapiens.GENCODE.V36-seed",destdir = "./BSgenome.gencode.v36.transcripts/")

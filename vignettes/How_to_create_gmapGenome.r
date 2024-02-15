#GmapGenome objects can be created from FASTA files or BSgenome objects. Below is 
#an example of how to create human hg38 gmapGenome from BSgenome.

library(gmapR)
library(BSgenome.Hsapiens.UCSC.hg38)

gmapGenomePath <- file.path("PATH to SAVE gmapGenome", "hg38")
gmapGenomeDirectory <- GmapGenomeDirectory(gmapGenomePath, create = TRUE)

gmapGenome <- GmapGenome(genome=Hsapiens, 
                         directory=gmapGenomeDirectory,
                         name="hg38", create=TRUE) 


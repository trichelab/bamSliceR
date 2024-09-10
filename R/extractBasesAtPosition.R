# Define the function to parse the CIGAR string
cigarToOps <- function(cigar) {
  ops <- unlist(strsplit(cigar, "(?<=\\D)", perl=TRUE))
  lengths <- as.integer(gsub("\\D", "", ops))
  types <- gsub("\\d", "", ops)
  return(list(lengths = lengths, types = types))
}

# parse the bamData from scanBAM() for readname, position, sequences, cigars of reads.
parseBamData = function( bamData)
{
  readNames = bamData$qname
  positions = bamData$pos
  sequences = bamData$seq
  cigars = bamData$cigar
  flag = bamData$flag
  return(list("readNames" = readNames,
	   "positions" = positions,
	   "sequences" = sequences,
           "cigars"   = cigars,
	   "flag" = flag) )  
}

# Function to calculate the relative position in the query sequence based on the CIGAR string
relativePosInQuery <- function(genomicPos, readPos, cigarOps) {
  queryPos <- 1
  refPos <- readPos
  for (i in seq_along(cigarOps$types)) {
    opLength <- cigarOps$lengths[i]
    opType <- cigarOps$types[i]
    if (opType %in% c("M", "X", "EQ")) {  # Alignment match (can be a mismatch or match)
      if (genomicPos <= refPos + opLength - 1) {
        queryPos <- queryPos + (genomicPos - refPos)
        queryPosEnd = queryPos
        if (genomicPos == refPos + opLength - 1)
        {
          # check if next cigar operation is "I" for insertion events.
          next_opLength <- cigarOps$lengths[i+1]
          next_opType <- cigarOps$types[i+1]
          if(!is.na(next_opType) ) {
	   if (next_opType == "I")
           {
             queryPosEnd = queryPosEnd + next_opLength
           }
           }
        }
        return(c(queryPos, queryPosEnd) )
      }
      queryPos <- queryPos + opLength
      refPos <- refPos + opLength
    } else if (opType == "I") {  # Insertion
      queryPos <- queryPos + opLength
    } else if (opType == "D") {  # Deletion
      refPos <- refPos + opLength
    } else if (opType == "S") {  # Soft clipping
      queryPos <- queryPos + opLength
    } else if (opType == "H") {  # Hard clipping
      # No change in queryPos or refPos for hard clipping
    } else if (opType == "N") {  # Skipped region
      refPos <- refPos + opLength
    } else if (opType == "P") {  # Padding
      # No change in queryPos or refPos for padding
    }
  }
  return(NA)
}

# Parse the CIGAR string
cigarToOps <- function(cigar) {
  ops <- unlist(strsplit(cigar, "(?<=\\D)", perl=TRUE))
  lengths <- as.integer(gsub("\\D", "", ops))
  types <- gsub("\\d", "", ops)
  return(list(lengths = lengths, types = types))
}

# Extract the base at the specified position for each read
# which here is a 1 row GRanges
#getBasesAtQueryPos = function(pos, cigar, seq) {
#  #pos = readBamData$pos
#  #cigar = readBamData$cigar
#  #seq = readBamData$seq
#  cigarOps <- cigarToOps(cigar)
#  queryPos <- relativePosInQuery(start(which), pos, cigarOps)
#  if (!is.na(queryPos[1]) && queryPos[1] > 0 && queryPos[1] <= nchar(seq)) {
#    return(substr(seq, queryPos[1], queryPos[2]))
#  } else {
#    return(NA)
#  }
#}

# For each Transcript, iterate through all the reads containd the variants, and extract the bases at positions
scanAllReads = function(parsedBamData, which)
{
 mapply( function(pos, cigar, seq, readName) {
  cigarOps <- cigarToOps(cigar)
  queryPos <- relativePosInQuery(start(which), pos, cigarOps)
  if (!is.na(queryPos[1]) && queryPos[1] > 0 && queryPos[1] <= nchar(as.character(seq))) {
    bs = as.character(substr(seq, queryPos[1], queryPos[2]))
    names(bs) = readName
    return (bs)
  } else {
    return(NA)
  } 
  } ,  parsedBamData$positions, parsedBamData$cigars, parsedBamData$sequences, parsedBamData$readNames ) -> bases
 bases_char = lapply(bases, as.character) %>% unlist()
 bases_char_df = data.frame(readn = names(bases_char), 
                            baseAtLocus = bases_char %>% unname() ) 
 return (bases_char_df)
}

#' given variants, extract the reads from BAMs
#'
#' @param bamFile BAM file name
#' @param which GRanges of coordiantes (Txs or Genomic)  of variants.
#'
#' @return vector of readNames for each alleles
#' @import Rsamtools
#' @export
#'
extractBasesAtPosition <- function(bamFile = "", which )
{
 gp <- ScanBamParam(which=which, 
		    what=c("qname", "pos", "seq", "cigar","flag"), flag = scanBamFlag(isDuplicate = FALSE) )
 bamData_list = scanBam(bamFile, param = gp)
 txs_reads_list = lapply(bamData_list, parseBamData) 
 gr_list = split(which)
 Map(scanAllReads, txs_reads_list, gr_list) -> results
 return (results)
}



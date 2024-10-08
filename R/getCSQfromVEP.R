# Deprecated ensemblVEP  
# use this for now
parseCSQ = function(x, VCFRowID=character(), ..., info.key="CSQ")
    {
        ## no 'info.key'
        if (!info.key %in% names(info(x)))
            return(rowRanges(x))

        hdr <- info(header(x))[info.key, "Description"]
        nms <- unlist(strsplit(strsplit(hdr, "Format: ")[[1]][2], "\\|"))
        ulst <- unlist(info(x)[[info.key]], use.names=FALSE)
        ## 'info.key' without data
        if (all(is.na(ulst))) {
            gr <- rowRanges(x)
            csq <-
                DataFrame(setNames(replicate(length(nms), character(0)), nms))
        } else {
            ## 'info.key' with data
            elt <- elementNROWS(info(x)[[info.key]])
            raw <- strsplit(ulst, "\\|")
            csq <- matrix(nrow=length(ulst), ncol=length(nms))
            for (i in 1:nrow(csq))
                csq[i, 1:length(raw[[i]])] <- raw[[i]]
            csq[!nzchar(csq)] <- NA
            colnames(csq) <- nms
            csq <- data.frame(csq, stringsAsFactors=FALSE)

            rd <- rowRanges(x)
            gr <- rd[rep(seq_along(rd), elt)]
            if (length(VCFRowID)) {
                if (any(no_match <- !VCFRowID %in% rownames(x)))
                    warning(paste0("rownames not found in 'x' : ",
                            paste(VCFRowID[no_match], collapse=",")))
                VCFRowID <- rep(match(rownames(x), VCFRowID), elt)
                csq <- DataFrame(VCFRowID=VCFRowID, csq)
            }
        }
        mcols(gr) <- csq
        genome(gr) <- genome(x)
        gr
    }

#' Parse the CSQ column in a VCF object returned from the Ensembl Variant effect Predictor (VEP)
#'
#' @param file path of the vcf file to read
#'
#' @return A VRanges object
#' @examples
#' x = 1+1
#' @export

getCSQfromVEP = function(file = "vep.vcf") 
{
 readVcf(file) -> vep_vcf
 parseCSQ(vep_vcf, info.key = "CSQ") -> csq
 return (csq)
}

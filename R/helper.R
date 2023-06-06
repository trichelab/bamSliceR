#' @import biomaRt

getHotSpot = function(gr, symbols, position, snp = TRUE)
{
    gr_hotspot = gr[which(gr$SYMBOL %in% symbols )]
    gr_hotspot = subset( gr_hotspot, POS == position)
    if (snp)
    {
        gr_hotspot = gr_hotspot[which(sapply(gr_hotspot$REFAA, nchar) == 1)]
        gr_hotspot = gr_hotspot[which(sapply(gr_hotspot$VARAA, nchar) == 1)]

    }
    return(gr_hotspot)
}

GDC_SAMPLE_TYPE = c("Recurrent Blood Derived Cancer - Bone Marrow",
"Primary Blood Derived Cancer - Bone Marrow",
"Primary Blood Derived Cancer - Peripheral Blood",
"Recurrent Blood Derived Cancer - Peripheral Blood",
"Bone Marrow Normal",
"Next Generation Cancer Model",
"Blood Derived Cancer - Bone Marrow, Post-treatment",
"Blood Derived Normal",
"Cell Lines",
"Blood Derived Cancer - Peripheral Blood, Post-treatment")

sampleTypeDecoder = function(df)
{
    results = case_when(
        df$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow" ~ "Relapse",
        df$sample_type == "Primary Blood Derived Cancer - Bone Marrow" ~ "Primary",
        df$sample_type == "Primary Blood Derived Cancer - Peripheral Blood" ~ "Primary",
        df$sample_type == "Recurrent Blood Derived Cancer - Peripheral Blood" ~ "Relapse",
        df$sample_type == "Bone Marrow Normal" ~ "Bone Marrow Normal",
        df$sample_type == "Next Generation Cancer Model" ~ "Next Generation Cancer Model",
        df$sample_type == "Blood Derived Cancer - Bone Marrow, Post-treatment" ~ "Post-treatment",
        df$sample_type == "Blood Derived Normal" ~ "Normal",
        df$sample_type == "Cell Lines" ~ "Cell Lines",
        df$sample_type == "Blood Derived Cancer - Peripheral Blood, Post-treatment" ~ "Post-treatment",
        TRUE ~ "Other"
    )
    return (results)
}

get_proteinPaint = function(gr_snps, disease = "AML")
{
    symbol = gr_snps$SYMBOL
    getRefSeqFromSymbol(symbol)[1,1] -> ref_seq
    time_point = sampleTypeDecoder (gr_snps)
    data.frame( gene = gr_snps$SYMBOL,
                chromosome = seqnames(gr_snps) %>% as.character() ,
                start = start(ranges(gr_snps)),
                aachange = gr_snps$CHANGE,
                class = "missense",
                disease = disease,
                refseq = ref_seq,
                patient = gr_snps$case_id,
                sampletype = gr_snps$sample_type,
                sample = gr_snps$sample,
                "mutant_in_tumor" = gr_snps$count.high.nm,
                "total_in_tumor" = gr_snps$count.high.nm.ref + gr_snps$count.high.nm ,
                VAF = gr_snps$count.high.nm / (gr_snps$count.high.nm.ref + gr_snps$count.high.nm),
                timepoint = time_point) -> df
    return (df)
}


getRefSeqFromSymbol = function(values)
{
    httr::set_config(httr::config(ssl_verifypeer = FALSE))

    # Connect to the Ensembl database through the biomaRt package
    ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

    # Define your gene symbol of interest
    gene_symbol <- values

    # Set the query attributes
    attributes <- c("refseq_mrna")

    # Set the filter for the gene symbol
    filter <- "external_gene_name"
    values <- gene_symbol

    # Execute the query and retrieve the data
    result <- getBM(
        attributes = attributes,
        filters = filter,
        values = values,
        mart = ensembl
    )

    return(result)
}



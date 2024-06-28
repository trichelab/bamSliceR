library(bamSliceR)
library(rtracklayer)

setwd("/varidata/research/projects/triche/Peter/leucegene/BAM/GSE67040/slice/minimap")
readRDS("tr.leu.gencode.v36.minimap2.rds") -> tr.leu.gencode.v36.minimap2
tr.leu.gencode.v36.minimap2.vr = stackSamples(VRangesList(tr.leu.gencode.v36.minimap2))
tr.leu.gencode.v36.minimap2.vr = saveVRinfo(tr.leu.gencode.v36.minimap2.vr)
file_meta = readRDS("/varidata/research/projects/triche/Peter/leucegene/cov/leucegene_transcriptomic_filemeta.rds" )
sample_names = sampleNames(tr.leu.gencode.v36.minimap2.vr) %>% as.character() %>% unique()
names(sample_names) = str_split(sample_names, "\\.", simplify = TRUE)[,1]
file_meta$file_name = sample_names[file_meta$UPC_ID]
file_meta$downloaded_file_name = sample_names[file_meta$UPC_ID]

annotateWithBAMinfo = function(tallied_reads, file_meta, bamfiles_names = NULL )
{
  if (is.null(bamfiles_names))
  {
    bamfiles_names = str_c(file_meta$sample,  "_",
                           file_meta$case_id, "_",
                           file_meta$file_name)
  }
  
  rownames(file_meta) = bamfiles_names
  mcols(tallied_reads) = cbind (mcols(tallied_reads),
                                file_meta[as.character(sampleNames(tallied_reads)), ] )
  return(tallied_reads)
}
tr.leu.gencode.v36.minimap2.vr.noFilter = annotateWithBAMinfo(tr.leu.gencode.v36.minimap2.vr, 
                                                              file_meta, 
                                                              bamfiles_names = file_meta$file_name)
seqlevels(tr.leu.gencode.v36.minimap2.vr.noFilter)  = as.character(seqnames(tr.leu.gencode.v36.minimap2.vr.noFilter)) %>% unique()
gencode.txs.file =  system.file("extdata", "gencode.v36.txs.annotation.subset.gff3", 
                            package = "bamSliceR")
leu.mini.txs = smartFilter(tr.leu.gencode.v36.minimap2.vr.noFilter, VAF_cutoff = 0.05, altDepth_cutoff = 5, 
                           gencode.file = gencode.txs.file)
fixIndelRefCounts(leu.mini.txs, dir = "./",mode = "INDEL", isFlank = FALSE, 
                  totalDepthOnly = FALSE, mc.cores = 30) -> leu.mini.txs.fixIndel
smartFilter(leu.mini.txs.fixIndel, VAF_cutoff = 0.05, altDepth_cutoff = 5, 
            gencode.file =gencode.txs.file) -> leu.mini.txs.fixIndel.filtered

fixMissingTxs(leu.mini.txs.fixIndel.filtered, gencode.file.txs = gencode.txs.file,
              bam.file.dir = "~/triche-secondary/Peter/leucegene/BAM/GSE67040/slice/minimap/") -> leu.mini.txs.fixIndel.filtered.fixtxs

fa = "/varidata/research/projects/triche/Peter/leucegene/GENCODEv36/gencode.v36.transcripts_header.fa"
getVariantAnnotationForTxs(gencode.file.txs = gencode.txs.file, seqSource = fa,
                           query.ranges = leu.mini.txs.fixIndel.filtered.fixtxs) -> gencode.leu.mini.txs.fixIndel.filtered.fixtxs.annotated

files() %>% GenomicDataCommons::select(c(fields,"file_name")  ) %>%
  GenomicDataCommons::filter( cases.project.project_id == projectId ) %>%
  GenomicDataCommons::filter( experimental_strategy == es  ) %>%
  GenomicDataCommons::filter( data_format == "BAM") -> qfiles
if (workflow != "")
{
  qfiles = qfiles %>%
    GenomicDataCommons::filter( analysis.workflow_type == workflow)
}

q = files() %>%
  GenomicDataCommons::select(available_fields('file_name')) %>%
  GenomicDataCommons::filter( cases.project.project_id == '') %>%
  GenomicDataCommons::filter( data_type == 'Aligned Reads') %>% 
  GenomicDataCommons::filter( experimental_strategy == 'RNA-Seq') %>% 
  GenomicDataCommons::filter( data_format == 'BAM') 
file_ids = q %>% results_all() %>% ids()

GenomicDataCommons::slicing("7d4d460b-5ac9-4a6f-83f8-b25f22c5dcad", regions = c("chr20:32358280-32439369", "chr2:25733703-25878537"), token = gdc_token(),
                            overwrite = TRUE, 
                            destination = "./fk.bam")

.gdc_post(endpoint = sprintf("slicing/view/%s", uuid),
          add_headers(`Content-type` = "application/json", Accept = "application/json"), 
          write_disk(destination, overwrite), if (progress) progress() else NULL, body = toJSON(body),
          token = token)
uri = "https://api.gdc.cancer.gov/slicing/view/7d4d460b-5ac9-4a6f-83f8-b25f22c5dcad"
token = "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJzdWIiOiI5ZmIyNTUxMGM2ZTBkYzcxYjJlODRjZWE3NzM0ZTFkZGI4NWE4OGVjNDNlMzg0ZGQwN2RjYTZkNjEwYzIxMzFiIiwiaWF0IjoxNzE5NTIxMDEzLCJleHAiOjE3MjIxMTMwMTMsIm9wZW5zdGFja19tZXRob2RzIjpbInNhbWwyIl0sIm9wZW5zdGFja19hdWRpdF9pZHMiOlsieGRfODgwLTZScVdkMTQ0UFUwaG5NUSJdLCJvcGVuc3RhY2tfZ3JvdXBfaWRzIjpbeyJpZCI6ImFmZDBkZDI0NzdhMzQ4ODdhMjA2OGYwYTgyMWQwNzI4In1dLCJvcGVuc3RhY2tfaWRwX2lkIjoiZXJhX2NvbW1vbiIsIm9wZW5zdGFja19wcm90b2NvbF9pZCI6InNhbWwyIn0.HJzoSrAAQr-s5x5aNZ9kDePwTjwHrnkjHoMctAFJNSKSFHPOuXPAXik6T7IDDDfwsClbnN8nVFL9z65c-MWv8Q"
regions = c("chr20:32358280-32439369","chr2:25733703-25878537")
body = body <- list(regions = regions)
body = toJSON(body)
POST(uri, add_headers(`X-Auth-Token` = token, Accept = "application/json",
                      `Content-Type` = "application/json"), ..., body = body, encode = "json")

"body"   "config" "encode" "handle" "url"

.gdc_post(endpoint = sprintf("slicing/view/%s",
                             uuid),
          write_disk(destination, overwrite), if (progress)
            progress()
          else NULL, body = toJSON(body), token = token)
devtools::install_github("trichelab/bamSliceR")
library(bamSliceR)
file_meta = getGDCBAMs(projectId = "TARGET-AML", es = "RNA-Seq", workflow = "STAR 2-Pass Genome")
head(file_meta)
target_genes_data = system.file("data", "gene_names.rds", package = "bamSliceR")
target_genes = readRDS(target_genes_data)
target_genes

#Get the vector of character() instead.
target_ranges_chars = getGenesCoordinates(target_genes, ret ="DF")
head(target_ranges_chars)

downloadSlicedBAMs(file_df = file_meta, regions = target_ranges_chars, dir = "./")

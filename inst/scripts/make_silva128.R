### Code to generate source files for a MgDb-object with data from the
### Silva v128
library(DECIPHER)
library(Biostrings)
library(metagenomeFeatures) ## Needed for the make_mgdb_sqlite function
library(tidyr)
library(tidyverse)
library(stringr)
## Database URL
db_root_url <- "ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus/"
taxa_url <- paste0(db_root_url, "taxonomy/97_otu_taxonomy.txt")
aligned_seq_url <- paste0(db_root_url, "rep_set_aligned/97_otus.fasta")
seq_url <- paste0(db_root_url, "rep_set/97_otus.fasta")
tree_url <- paste0(db_root_url, "trees/97_otus_unannotated.tree")

## Downloaded files
taxa_file <- tempfile()
seq_file <- tempfile()
tree_file <- "../extdata/gg13.8_97.tre"
rnacentral_file <- "../../rnacentral/greengenes.tsv"

## MD5 check sums from initial download
taxa_md5 <- "56ef15dccf2e931ec173f4f977ed649b"
seq_md5 <- "50b2269712b3738afb41892bed936c29"
tree_md5 <- "e42c0cb6e4eaf641742cd5767c4a0101"

## MgDb database files name
db_file <- "../extdata/silva.sqlite"
metadata_file <- "../extdata/silva_metadata.RData"

### Download database files ####################################################
download_db <- function(url, file_name, md5){
    ## Downloade file and check to make sure MD5 checksum matches checksum for
    ## previously downloaded version

    download.file(url,file_name)
    new_md5 <- tools::md5sum(file_name)
    if (md5 != new_md5) warning("checksum does not match downloaded file.")
}

## Taxa Data
download_db(taxa_url, taxa_file, taxa_md5)

## Seq Data
download_db(seq_url, seq_file, seq_md5)

## Tree Data
download_db(tree_url, tree_file, tree_md5)

##Directly accessing the files locally
taxa_file <- "../../silva_db/taxmap_slv_ssu_ref_128.txt"
seq_file <- "../../silva_db/SILVA_128_SSURef_tax_silva.fasta"
tree_file <- "../../silva_db/tax_slv_ssu_128.tre"
taxmap_file <- "../../silva_db/tax_slv_ssu_128.txt"
rnacentral_file <- "../../silva_db/silva.tsv"

### Create SQLite DB with Taxa and Seq Data ####################################
expand_path <- function(path){
    split_path <- unlist(str_split(path, ";"))
    path_df <- data_frame()
    for (i in 1:(length(split_path) - 1)) {
        path_frag <- paste0(paste(split_path[1:i], collapse = ";"), ";")
        # path_list <- c(path_list, path_frag)
        path_df <- bind_rows(path_df, data_frame(taxa = split_path[i], V1 = path_frag))
    }

    path_df
}

### Parse silva taxonomy
parse_silva <- function(taxa_file, rnacentral_file, taxamap_file){
    taxa <- read.delim(taxa_file, stringsAsFactors = FALSE, header = TRUE)
    taxamap <- read.delim(taxmap_file, stringsAsFactors = FALSE, header = FALSE)
    expanded_taxa <- taxa[1:50,] %>%
        as_tibble() %>%
        select(primaryAccession, path) %>%
        mutate(expanded_path = map(path, expand_path)) %>%
        select(-path) %>%
        unnest()
    mergenewtaxa <- inner_join(x=expanded_taxa, y=taxamap)
    wide_taxa_df <- mergenewtaxa %>%
        select(primaryAccession, taxa, V3) %>%
        spread(V3,taxa,fill = NA)
    taxa_df <- data.frame(wide_taxa_df$primaryAccession,wide_taxa_df$domain,wide_taxa_df$phylum,wide_taxa_df$class, wide_taxa_df$order,wide_taxa_df$family, wide_taxa_df$genus)
    taxa_df <- taxa_df[match(taxa$primaryAccession[1:50], taxa_df$wide_taxa_df.primaryAccession),]
    taxa_df <- cbind(taxa_df, taxa$organism_name[1:50])
    colnames(taxa_df) <- c("Keys","Kingdom","Phylum","Class","Order","Family","Genus","Species")
    ## Add RNAcentral ID and NCBItax ID - The primaryAccession field in taxa_df differs from RNAmap identifier.
    # rnamap <- read.delim(rnacentral_file, stringsAsFactors = FALSE, header = FALSE)
    # rnamap_new <- rnamap[match(taxa$primaryAccession,rnamap$V3),]
    # taxa_df <-cbind(taxa_df,rnamap_new$V1[1:50],rnamap_new$V4[1:50])
    # colnames(taxa_df) <- c("Keys","Kingdom","Phylum","Class","Ord","Family","Genus","Species","RNAcentralID","NCBItaxID")
    
    ## Return as a data.frame
    data.frame(taxa_df)
}

taxa_tbl <- parse_silva(taxa_file, rnacentral_file)

seqs <- Biostrings::readDNAStringSet(seq_file)

metagenomeFeatures::make_mgdb_sqlite(db_name = "silva",
                                     db_file = db_file,
                                     taxa_tbl = taxa_tbl,
                                     seqs = seqs)



### Database Metadata ##########################################################
metadata <- list(ACCESSION_DATE = date(),
                 URL = "ftp://greengenes.microbio.me/greengenes_release/gg_13_8_otus",
                 DB_TYPE_NAME = "GreenGenes",
                 DB_VERSION = "13.8 97% OTUS",
                 DB_TYPE_VALUE = "MgDb",
                 DB_SCHEMA_VERSION = "2.0")

save(metadata, file = metadata_file)

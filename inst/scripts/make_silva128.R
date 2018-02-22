### Code to generate source files for a MgDb-object with data from the
### Silva v128
library(DECIPHER)
library(Biostrings)
library(metagenomeFeatures) ## Needed for the make_mgdb_sqlite function
library(tidyr)
library(tidyverse)
library(stringr)
library(R.utils)
## Database URL
db_root_url <- "https://www.arb-silva.de/fileadmin/silva_databases/release_128/Exports/"
taxa_url <- paste0(db_root_url, "taxonomy/taxmap_slv_ssu_ref_128.txt.gz")
taxmap_url <- paste0(db_root_url, "taxonomy/tax_slv_ssu_128.txt")
aligned_seq_url <- paste0(db_root_url, "SILVA_128_SSURef_tax_silva_full_align_trunc.fasta.gz")
seq_url <- paste0(db_root_url, "SILVA_128_SSURef_tax_silva.fasta.gz")
tree_url <- paste0(db_root_url, "taxonomy/tax_slv_ssu_128.tre")
rnacentral_url <- "ftp://ftp.ebi.ac.uk/pub/databases/RNAcentral/releases/8.0/id_mapping/database_mappings/silva.tsv"

## Downloaded files
taxa_file <- tempfile()
taxagz_file <- tempfile()
seq_file <- tempfile()
seqgz_file <- tempfile()
rnacentral_file <- tempfile()
taxmap_file <- tempfile()
tree_file <- "../extdata/silva128.tre"

## MD5 check sums from initial download
taxa_md5 <- "766e219360430c7117fd6d76a4525b04"
taxmap_md5 <- "fe5ff3ce9681edc8cf523a4706e76243"
seq_md5 <- "e6feea458c5f1194192c0710c0f40938"
tree_md5 <- "81fe22e0f57a478bda5cea67857baade"

## MgDb database files name
db_file <- "../extdata/silva128.sqlite"
metadata_file <- "../extdata/silva128_metadata.RData"

### Download database files ####################################################
download_db <- function(url, file_name, md5){
    ## Downloade file and check to make sure MD5 checksum matches checksum for
    ## previously downloaded version

    download.file(url,file_name)
    new_md5 <- tools::md5sum(file_name)
    if (md5 != new_md5) warning("checksum does not match downloaded file.")
}

## Taxa Data
download_db(taxa_url, taxagz_file, taxa_md5)
gunzip(filename = taxagz_file, destname = taxa_file)

## Taxonomy mapping Data
download_db(taxmap_url, taxmap_file, taxmap_md5)

## Seq Data
download_db(seq_url, seqgz_file, seq_md5)
gunzip(filename = seqgz_file, destname = seq_file)

## Tree Data
download_db(tree_url, tree_file, tree_md5)

##Directly downloading RNAcentral mapping file
download.file(rnacentral_url, rnacentral_file)

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
    colnames(taxa_df) <- c("Keys","Kingdom","Phylum","Class","Ord","Family","Genus","Species")
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
                 URL = "https://www.arb-silva.de/fileadmin/silva_databases/release_128/",
                 DB_TYPE_NAME = "SILVA",
                 DB_VERSION = "128 99% OTUS",
                 DB_TYPE_VALUE = "MgDb",
                 DB_SCHEMA_VERSION = "2.0")

save(metadata, file = metadata_file)

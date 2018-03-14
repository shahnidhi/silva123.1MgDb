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

### Parse silva taxonomy
parse_silva <- function(taxa_file, rnacentral_file, taxamap_file){
    taxa <- read.delim(taxa_file, stringsAsFactors = FALSE, header = TRUE)
    taxamap <- read.delim(taxmap_file, stringsAsFactors = FALSE, header = FALSE)
    tpath <- taxa$path

    leaf <- lapply(taxamap$V1, function(n) {
        last_node <- str_split(n, ";")
        last_node[[1]][length(last_node[[1]]) - 1]
    })

    taxamap$leaf <- leaf

    taxa_sub <- unique(data.frame(node = unlist(leaf),
                                  level = taxamap$V3))

    taxa_list <- taxa_sub$level
    names(taxa_list) <- taxa_sub$node

    hierarchy <- lapply(tpath, function(path) {
        req <- c("domain", "phylum", "class", "order", "family", "genus")

        last_node <- str_split(path, ";")
        nodes <- last_node[[1]][1:length(last_node[[1]]) -1]
        level <- lapply(nodes, function(node) {
            as.character(taxa_list[[node]][1])
        })

        tlist <- nodes
        names(tlist) <- level
        out <- tlist[req]
        names(out) <- req
        out
    })

    df_hierarchy <- map_df(hierarchy, ~bind_rows(.))

    taxa$domain <- df_hierarchy$domain
    taxa$phylum <- df_hierarchy$phylum
    taxa$class <- df_hierarchy$class
    taxa$order <- df_hierarchy$order
    taxa$family <- df_hierarchy$family
    taxa$genus <- df_hierarchy$genus

    taxa$path <- NULL
    taxa$taxid <- NULL
    colnames(taxa) <- c("Accession", "start", "stop","Species", "Kingdom", "Phylum", "Class", "Ord", "Family", "Genus" )
    taxa$Keys <- paste0(taxa_tbl$Accession,".", taxa_tbl$start,".", taxa_tbl$stop)
    taxa_tbl_new <- df[,c(which(colnames(df)=="Keys"),which(colnames(df)!="Keys"))]

    ## Add RNAcentral ID and NCBItax ID - The primaryAccession field in taxa_df differs from RNAmap identifier.
    # rnamap <- read.delim(rnacentral_file, stringsAsFactors = FALSE, header = FALSE)
    # rnamap_new <- rnamap[match(taxa$primaryAccession,rnamap$V3),]
    # taxa_df <-cbind(taxa_df,rnamap_new$V1[1:50],rnamap_new$V4[1:50])
    # colnames(taxa_df) <- c("Keys","Kingdom","Phylum","Class","Ord","Family","Genus","Species","RNAcentralID","NCBItaxID")

    ## Return as a data.frame
    data.frame(taxa_tbl_new)
}

taxa_tbl <- parse_silva(taxa_file, rnacentral_file)

rna_seqs <- Biostrings::readRNAStringSet(seq_file)
seqs <- Biostrings::DNAStringSet(rna_seqs)
seq_names_final <- vapply(strsplit(names(seqs)," ",fixed=TRUE), `[`, 1, FUN.VALUE=character(1))
names(seqs) <- seq_names_final
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

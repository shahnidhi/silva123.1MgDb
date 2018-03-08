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
    # unlist(out)
    out
})

library(tidyverse)

df_hierarchy <- map_df(hierarchy, ~bind_rows(.))

taxa$domain <- df_hierarchy$domain
taxa$phylum <- df_hierarchy$phylum
taxa$class <- df_hierarchy$class
taxa$order <- df_hierarchy$order
taxa$family <- df_hierarchy$family
taxa$genus <- df_hierarchy$genus

write.table(taxa, file="../extdata/taxa_added_map.txt", quote=F, sep=",", row.names=F)

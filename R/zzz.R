###
### Load silva123.1MgDb into namespace
###

.onAttach <- function(libname, pkgname){

    db_file <- system.file("extdata", "silva128.1.sqlite",
                                package = pkgname, lib.loc = libname)

    metadata_file <- system.file("extdata", "silva128.1_metadata.RDS",
                                package = pkgname, lib.loc = libname)

    ## Note no tree for silva123.1

    if (!file.exists(db_file) | !file.exists(metadata_file)) {
        packageStartupMessage("SILVA 128.1 database data not present, use `get_silva128.1.R` In the package inst/scripts directory to download the database into the package inst/extdata/ directory and reinstall the package")
    }
}

.onLoad <- function(libname, pkgname){
    ns <- asNamespace(pkgname)

    db_file <- system.file("extdata", "silva128.1.sqlite",
                                package = pkgname, lib.loc = libname)

    metadata_file <- system.file("extdata", "silva128.1_metadata.RDS",
                                package = pkgname, lib.loc = libname)

    ## Add Tree data

    metadata <- readRDS(metadata_file)

    ## initiate new MgDB object
    slvMgDb <- metagenomeFeatures::newMgDb(db_file = db_file,
                      tree = NULL,
                      metadata = metadata)

    assign("slv128.1MgDb", slvMgDb, envir = ns)
    namespaceExport(ns, "slv128.1MgDb")

}

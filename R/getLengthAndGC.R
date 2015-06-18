############################################################
# 
# author: Ludwig Geistlinger
# date: 2015-03-10 13:32:37
# 
# descr: get.gene.length.and.gc.content
# update: 2015-06-14 exonic sequences 
############################################################

#
# @input:
#   - id: one or more gene IDs (ensembl or entrez)
#   - org: organism three letter code, e.g. 'hsa' for 'Homo sapiens'
#   - mode: 1. biomart (supports all ensembl organisms, but might be a time-consuming)
#           2. org.db (based on BioC annotation, which is much faster but only 
#                    for organisms with a respective TxDb, BSgenome, and OrgDb package)
getGeneLengthAndGCContent <- function(id, org, mode=c("biomart", "org.db"))
{
    id.type <- .autoDetectGeneIdType(id[1])
    if(is.na(id.type))
        stop("Only ENTREZ or ENSEMBL gene IDs are supported.")
    
    mode <- match.arg(mode)
    inp.id <- id

    # (a) based on BioC annotation utilities:
    #       (0) OrgDb: map between identifiers (if necessary)
    #       (1) TxDB: get genomic coordinates of genes
    #       (2) BSgenome: get sequences of genomic coordinates
    #
    if(mode=="org.db")
    {
        # check for TxDb package
        txdb.pkg <- .org2pkg(org, type="TxDb")
        .isAvailable(txdb.pkg, type="TxDb")

        # check for BSgenome package
        bsgen.pkg <- .org2pkg(org, type="BSgenome") 
        .isAvailable(bsgen.pkg, type="BSgenome")

        txdb.spl <- unlist(strsplit(txdb.pkg, "\\."))
        txdb.id.type <- txdb.spl[length(txdb.spl)]
        if(txdb.id.type == "ensGene") {
          txdb.id.type <- "ensembl"
        } else if(txdb.id.type == "knownGene") {
          txdb.id.type <- "entrez"
        } else if(txdb.id.type == "sgdGene") {
          txdb.id.type <- "sgd"
        } else {
          stop(paste("TxDb does not use ENSEMBL or ENTREZ gene IDs"))
        }

        # (0) map ensembl <-> entrez, 
        # if given id.type is entrez, but Txdb uses ensembl (or vice versa)
        if(id.type != txdb.id.type)
        {
            orgdb.pkg <- .org2pkg(org)
            .isAvailable(orgdb.pkg)
            orgdb.pkg <- get(orgdb.pkg)
            id.map <- mapIds(orgdb.pkg, keys=id, 
                column=ifelse(id.type == "entrez", "ENSEMBL", "ENTREZID"), 
                keytype=ifelse(id.type == "entrez", "ENTREZID", "ENSEMBL"))
            id <- id.map[!is.na(id.map)]
        }
        
        # (1) get genomic coordinates
        txdb.pkg <- get(txdb.pkg)
        coords <- exonsBy(txdb.pkg, by="gene")
        id <- id[id %in% names(coords)]
        coords <- lapply(coords[id], reduce)
        len <- sapply(coords, function(x) sum(width(x)))
        
        # (2) get sequences
        bsgen.pkg <- get(bsgen.pkg)
        seqs <- lapply(coords, function(x) getSeq(bsgen.pkg, x))
        gc.cont <- sapply(seqs, function(s)
          mean(sapply(s, function(ss) sum(alphabetFrequency(ss, as.prob=TRUE)[c("C","G")]))))
    }
    # (b) based on BioMart
    #
    #
    else
    {
        id.type <- paste0(id.type, ifelse(id.type=="entrez", "gene", "_gene_id"))
        
        # setting mart
        message("Connecting to BioMart ...")
        ensembl <- useMart("ensembl")
        ds <- listDatasets(ensembl)[,"dataset"]
        ds <- grep(paste0("^", org), ds, value=TRUE)
        if(length(ds) == 0) 
            stop(paste("Mart not found for:", org))
        else if(length(ds) > 1)
        {
            message("Found several marts")
            sapply(ds, function(d)
                message(paste(which(ds==d), d, sep=": ")))
            n <- readline(paste0("Choose mart (1-", length(ds),") : "))
            ds <- ds[as.integer(n)]
        }

        ensembl <- useDataset(ds, mart=ensembl)

        message( paste0( "Downloading sequence", 
            ifelse(length(id) > 1, "s", ""), " ..."))
        if(length(id) > 100) message("This may take a few minutes ...")

        # download sequence
        # (1) get exon coordinates
        attrs <- c(id.type, "ensembl_exon_id", 
            "chromosome_name", "exon_chrom_start", "exon_chrom_end")
        coords <- getBM(filters=id.type, attributes=attrs, values=id, mart=ensembl)
        id <- unique(coords[,id.type])
        coords <- sapply(id, 
            function(i)
            { 
                i.coords <- coords[coords[,1]== i, 3:5]
                g <- GRanges(i.coords[,1], IRanges(i.coords[,2],i.coords[,3]))
                return(g)
            })
        coords <- lapply(coords[id], reduce)
        len <- sapply(coords, function(x) sum(width(x)))
        
        # (2) get genes and sequences
        sel <- c(id.type, "start_position", "end_position")
        gene.pos <- getGene(id=id, type=id.type, mart=ensembl)
        gene.pos <- gene.pos[,sel]
        
        gene.seqs <- getSequence(id=id, 
            type=id.type, seqType="gene_exon_intron", mart=ensembl) 

        # (3) get exonic sequences and correspondig GC content
        gc.cont <- sapply(id, 
            function(i)
            {
                # exon coordinates, gene position & sequence for current id i
                ecoords <- coords[[i]]
                gpos <- gene.pos[gene.pos[,id.type] == i, 
                        c("start_position", "end_position")]
                gseq <- DNAString(
                    gene.seqs[gene.seqs[,id.type] == i, "gene_exon_intron"])
                
                # exon coordinates relative to gene position
                start <- start(ranges(ecoords)) - gpos[1,1] + 1 
                end <- end(ranges(ecoords)) - gpos[1,1] + 1
                eseq <- gseq[IRanges(start, end)]
                gc.cont <- sum(alphabetFrequency(eseq, as.prob=TRUE)[c("C","G")])
                return(gc.cont)
            }
        )
    }

    res <- cbind(len, gc.cont)
    colnames(res) <- c("length", "gc")
    rownames(res) <- id

    # (4) order according to input ids
    if(mode == "org.db")
        if(id.type != txdb.id.type)
            rownames(res) <- names(id)

    not.found <- !(inp.id %in% rownames(res))
    na.col <- rep(NA, sum(not.found))
    rn <- c(rownames(res), inp.id[not.found])
    res <- rbind(res, cbind(na.col, na.col))
    rownames(res) <- rn
    res <- res[inp.id,]
    return(res)
}


.isAvailable <- function(pkg, type="annotation")
{
    if(!(pkg %in% .packages(all.available=TRUE)))
    {   
        message(paste0("Corresponding ", type,  " package not found: ", 
            pkg, "\nMake sure that you have it installed."))
        choice <- readline("Install it now? (y/n): ")
        if(choice == "y")
        {   
            source("http://bioconductor.org/biocLite.R")
            biocLite(pkg)
        }   
        else stop(paste("Package", pkg, "is not available"))
    }   
    require(pkg, character.only = TRUE)
}

.autoDetectGeneIdType <- function(id)
{
    type <- NA
    if(grepl("^[Ee][Nn][Ss][A-Za-z]{0,3}[Gg][0-9]+", id)) type <- "ensembl"
    else if(grepl("^[0-9]+$", id)) type <- "entrez"
    else if(grepl("^[Yy][A-Za-z]{2}[0-9]{3}[A-Za-z]", id)) type <- "sgd"
    else if(grepl("^[Aa][Tt][0-9][A-Za-z][0-9]{5}", id)) type <- "tair"
    return(type)
}

.getOrgIdType <- function(org)
{
    it <- "eg"
    if(org == "At") it <- "tair"
    else if(org == "Pf") it <- "plasmo"
    else if(org == "Sc") it <- "sgd"
    return(it)
}   

.supportedOrganisms <- function() sub(".db0$", "", .availableOrgPkgs())

.availableOrgPkgs <- function(type=c("OrgDb", "TxDb", "BSgenome"), local=TRUE)
{
    if(local) pkgs <- .packages(all.available=TRUE)
    else pkgs <- available.packages(paste0("http://bioconductor.org/",
        "packages/release/data/annotation/src/contrib"))[, "Package"]
    
    type <- match.arg(type)
    org.string <- "^org.[A-z][a-z]+.[a-z]+.db$"
    if(type == "TxDb") 
        org.string <- "^TxDb.[A-Z][a-z]+.UCSC.[a-z]{2}[A-Za-z]*[0-9]{1,3}.[a-z]{3,5}Gene$"
    else if(type == "BSgenome") 
        org.string <- "^BSgenome.[A-Z][a-z]+.UCSC.[a-z]{2}[A-Za-z]*[0-9]{1,3}$"
    org.pkgs <- grep(org.string, pkgs, value=TRUE)
    names(org.pkgs) <- NULL 
    return(org.pkgs)
}

.org2pkg <- function(org, type=c("OrgDb", "TxDb", "BSgenome"))
{
    type <- match.arg(type)

    SPECIES <- rbind(
        c("anopheles", "Anopheles gambiae", "Ag", "aga", "anoGam", "7165"),
        c("arabidopsis", "Arabidopsis thaliana", "At", "ath", NA, "3702"),
        c("bovine", "Bos taurus", "Bt", "bta", "bosTau", "9913"),
        c("canine", "Canis familiaris", "Cf", "cfa", "canFam", "9615"),
        c("chicken", "Gallus gallus", "Gg", "gga", "galGal", "9031"), 
        c("chimp", "Pan troglodytes", "Pt", "ptr", "PanTro", "9598"),
        c("ecoliK12", "Escherichia coli K12", "EcK12", "eco", NA, "562,83333,511145"), 
        c("ecoliSakai", "Escherichia coli Sakai", "EcSakai", "ecs", NA, "83334"),
        c("fly", "Drosophila melanogaster", "Dm", "dme", "dm", "7227"),
        c("human", "Homo sapiens", "Hs", "hsa", "hg", "9606"),
        c("malaria", "Plasmodium falciparum", "Pf", "pfa", NA, "5833"),
        c("mouse", "Mus musculus", "Mm", "mmu", "mm", "10090"),
        c("pig", "Sus scrofa", "Ss", "ssc", "susScr", "9823"),
        c("rat", "Rattus norvegicus", "Rn", "rno", "rn", "10116"), 
        c("rhesus", "Macaca mulatta", "Mmu", "mcc", "rheMac", "9544"),  
        c("worm", "Caenorhabditis elegans", "Ce", "cel", "ce", "6239"),
        c("xenopus", "Xenopus laevis", "Xl", "xla", "NA", "8355"),
        c("yeast", "Saccharomyces cerevisiae", "Sc", "sce", "sacCer", "4932,559292"),
        c("zebrafish", "Danio rerio", "Dr", "dre", "danRer", "7955")
    )
    colnames(SPECIES) <- c("common", "tax", "bioc", "kegg", "ucsc", "ncbi")
    

    # org specification via 
    # (a) 3-letter code, e.g. 'hsa' 
    # (b) genome assembly, e.g. 'hg38'
    is.genome <- sub("[0-9]+$", "", org) %in% SPECIES[,"ucsc"]
    if(is.genome)
    {
        ucsc.id <- org
        i <- grep(sub("[0-9]+$", "", org), SPECIES[,"ucsc"]) 
        bioc.id <- SPECIES[i, "bioc"]
    }
    else
    {
        ind <- apply(SPECIES, 1, function(r) org %in% r)
        if(any(ind)) i <- which(ind)[1]
        else stop(paste0("unrecognized organism ID \'", org, "\'"))
        bioc.id <- SPECIES[i, "bioc"]
        ucsc.id <- SPECIES[i, "ucsc"]
    }

    # TxDB, BSgenome, or OrgDB package?
    if(type %in% c("TxDb", "BSgenome"))
    {
        pkg.string <- paste0("^", type, ".", bioc.id, "[a-z]+.UCSC.", ucsc.id)
        pkg <- grep(pkg.string, .availableOrgPkgs(type), value=TRUE)
        if(length(pkg) == 0)
            pkg <- grep(pkg.string, .availableOrgPkgs(type, local=FALSE), value=TRUE)
        if(length(pkg) == 0)
            stop(paste("No corresponding", type, "package for", org))
        else if(length(pkg) > 1)
        {
            message("Found several genome assemblies")
            sapply(pkg, function(p) 
                message(paste(which(pkg==p), p, sep=": ")))
            n <- readline(paste0("Choose assembly (1-", length(pkg),") : "))
            pkg <- pkg[as.integer(n)]

            #message("Found several genome assemblies")
            #message(paste("Using latest:", pkg))
            #ver <- sapply(pkg, 
            #    function(p)
            #    {
            #        spl <- unlist(strsplit(p, "\\."))
            #        ind <- length(spl)
            #        if(type == "TxDb") ind <- ind - 1
            #        ass <- spl[ind]
            #        ver <- sub("^[a-zA-Z]+", "", ass)
            #        return(as.integer(ver))
            #    })
            #pkg <- pkg[which.max(ver)]
        }
    }
    else
    {
        id.type <- .getOrgIdType(bioc.id)
        pkg <- paste("org", bioc.id, id.type, "db", sep=".")
    }
    return(pkg)
}



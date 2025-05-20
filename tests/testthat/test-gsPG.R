test_that("gsPG runs without error and returns numeric output", {
    kegg <- getKEGG(species="Human")
    gene1 <- unique(as.vector(unlist(kegg[1:5])))
    gene2 <- unique(as.vector(unlist(kegg[6:length(kegg)])))
    gene1 <- rep(gene1,sample(1:10,length(gene1),replace=TRUE))
    gene2 <- rep(gene2,sample(1:10,length(gene2),replace=TRUE))
    p1 <- runif(length(gene1))*(1e-6)
    p2 <- runif(length(gene2))
    geneid <- c(gene1,gene2)
    p <- c(p1,p2)
    Name <- paste0("cg",1:length(p))
    probe.p <- data.frame(Name=Name,p=p)
    GeneProbeTable <- data.frame(Name=Name,entrezid=geneid)
    dat <- matrix(runif(length(p)*100),ncol=100)
    rownames(dat) <- Name
    gsPG(probe.p=probe.p,Data4Cor=dat,GeneProbeTable=GeneProbeTable,
    gSetName="KEGG",species="Human",ncore=1)

    # Check that the output file is created
    expect_true(file.exists("gsPG_KEGG_KEGG.csv"))

    # Optional: load and test contents
    result <- read.csv("gsPG_KEGG_KEGG.csv")
    expect_true("p" %in% colnames(result))
    expect_gt(nrow(result), 100)

    # Clean up
    file.remove("gsPG_KEGG_KEGG.csv")
})

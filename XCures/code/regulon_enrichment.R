##### Cluster Pathway Enrichment Script v 1.4 #####
## Serdar Turkarslan                            ##
## Institute for Systems Biology, 2020          ##
## Last update: 08/20/2020                      ##
## Input set of genes with gene symbol either   ##
## as vector or data frame with Ensembl/Symbol  ##
## Enrichments are performed against different  ##
## MSigSdb signatures                           ##
#################################################

regulon_enrichment <- function(my.genes){
  #library(pathfindR)
  library(clusterProfiler)
  library(msigdbr)
  library(enrichplot)
  library(ReactomePA)
  library(AnnotationDbi)
  library(org.Hs.eg.db)

  ## Input gene list could be data frame or vector of genes
  if(is.vector(my.genes)){
    genes.input <- my.genes
    #entrez.ids <- AnnotationDbi::select(org.Hs.eg.db, keys= my.genes, column="ENTREZID", keytype="SYMBOL", multiVals="first")
  }
  if(is.data.frame(my.genes)){
    genes.input <- my.genes$SYMBOL
    #entrez.ids <- AnnotationDbi::select(org.Hs.eg.db, keys= input.genes$ENSEMBL, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
  }

  #entrez.ids <- AnnotationDbi::select(org.Hs.eg.db, keys= input.genes$ENSEMBL, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
  ## Hallmarks Enrichment (Uses MSIGDB Hallmark signatures)
  h_df = msigdbr(category = "H")
  h_t2g = h_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  edo.h <- enricher(gene = genes.input, TERM2GENE = h_t2g)

  ## Immunologic Signatures (MSIGDb)
  # c7_df = msigdbr(category = "C7")
  # c7_t2g = c7_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  # edo.c7 <- enricher(gene = genes.input, TERM2GENE = c7_t2g)
  #
  #   edo.c7 <- edo.c7 %>%
  #     mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))


  ## Oncogenic Signatures (MSIGDb)
  # c6_df = msigdbr(category = "C6")
  # c6_t2g = c6_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  # edo.c6 <- enricher(gene = genes.input, TERM2GENE = c6_t2g)
  #
  #   edo.c6 <- edo.c6 %>%
  #     mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))


  ## Curated datasets (MSIGDb)
  # c2_df = msigdbr(category = "C2")
  # c2_t2g = c2_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  # edo.c2 <- enricher(gene = genes.input, TERM2GENE = c2_t2g)
  #
  #   edo.c2 <- edo.c2 %>%
  #     mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))


  ## Transcriptin factor target gene set
  # c3_df = msigdbr(category = "C3")
  # c3_t2g = c3_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  # edo.c3 <- enricher(gene = genes.input, TERM2GENE = c3_t2g)
  #
  #   edo.c3 <- edo.c3 %>%
  #    mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))


  ## KEGG  enrichments (MSIGDb)
  # kegg_df = msigdbr(category = "C2", subcategory = "CP:KEGG")
  # kegg_t2g = kegg_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  # edo.kegg <- enricher(gene = genes.input, TERM2GENE = kegg_t2g)
  #

  # ## Curated datasets (MSIGDb)
  # cm_df = msigdbr(category = "C4", subcategory = "CM")
  # cm_t2g = cm_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  # edo.cm <- enricher(gene = genes.input, TERM2GENE = cm_t2g)
  #
  #   edo.cm <- edo.cm %>%
  #     mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))


  ## GO BP (MSIGDb)
  gobp_df = msigdbr(category = "C5", subcategory = "BP")
  gobp_t2g = gobp_df %>% dplyr::select(gs_name, human_gene_symbol) %>% as.data.frame()
  edo.gobp <- enricher(gene = genes.input, TERM2GENE = gobp_t2g)

    edo.gobp <- edo.gobp %>%
      mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))

    #reactome
    my.entrez <- AnnotationDbi::select(org.Hs.eg.db, keys= genes.input, column="ENTREZID", keytype="SYMBOL", multiVals="first")
    my.entrez <- my.entrez$ENTREZID
    #my.entrez <- as.character(sort(as.numeric(my.entrez), decreasing = T))
    reactome.enrich <- enrichPathway(gene = my.entrez, pvalueCutoff = 0.05, pAdjustMethod = "BH")

### Dotplots for each enrichments
  p1 <- dotplot(edo.h,font.size = 8)
  #p2 <- dotplot(edo.c7,font.size = 8)
  #p3 <- dotplot(edo.c6,font.size = 8)
  #p4 <- dotplot(edo.c2,font.size = 8)
  #p5 <- dotplot(edo.cm,font.size = 8)
  p6 <- dotplot(edo.gobp,font.size = 8)
  #p7 <- dotplot(edo.c3,font.size = 8)
  p8 <- dotplot(reactome.enrich,font.size = 8)

## if no enrichment handle empty results
    if(length(data.frame(edo.h)$ID) <= 2){
    q2 = ""
  } else{
    q2 <- try(cnetplot(edo.h, font.size=6))
  }

  # if(length(data.frame(edo.kegg)$ID) <= 2){
  #   q3 = ""
  # } else{
  #   q3 <- try(cnetplot(edo.kegg, font.size=6))
  # }

  ## Add URL links to tables
  try(
    edo.h <- edo.h %>%
      mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))
  )

  # try(
  #   edo.kegg <- edo.kegg %>%
  #     mutate(ID = paste("<a href='https://www.gsea-msigdb.org/gsea/msigdb/cards/", ID, "'>", ID,"</a>",sep = ""))
  # )

  ## Enrichment tables for the above comparisons.
  d1 <- as.data.frame(edo.h)

 # d2 <- as.data.frame(edo.c7)

  #d3 <- as.data.frame(edo.c6)

  #d4 <- as.data.frame(edo.c2)

  #d5 <- as.data.frame(edo.cm)

  d6 <- as.data.frame(edo.gobp)

  #d7 <- as.data.frame(edo.c3)

  d8 <- as.data.frame(reactome.enrich)

  ## Return result list
  return(list(hallmark=p1, reactome=p8, GO.BP=p6, table1=d1, table6=d6, table8=d8))

}

library(tidyverse)
library(httr)
library(stringr)
library(forcats)
library(biomaRt)
library(Biostrings)

human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

lrpFam <- tibble(geneName = as.factor(c("LDLR", "VLDLR", "LRP8", "LRP4", "LRP1", "LRP1B", "LRP2",
                                        "LRPAP1", "A2M", "APOE", "F8", "PDGFB", "TGFB1")),
                 altGeneName = as.factor(c("LDLR", "VLDLR", "APOER2", "LRP4", "LRP1", "LRP1B", "Megalin",
                                           "RAP", "A2M", "APOE", "FVIII", "PDGFB", "TGFB1")),
                 ensembl_gene_id = c("ENSG00000130164",
                                     "ENSG00000147852",
                                     "ENSG00000157193",
                                     "ENSG00000134569",
                                     "ENSG00000123384",
                                     "ENSG00000168702",
                                     "ENSG00000081479",
                                     "ENSG00000163956",
                                     "ENSG00000175899",
                                     "ENSG00000130203",
                                     "ENSG00000185010",
                                     "ENSG00000100311",
                                     "ENSG00000105329"))

attrEndings <- c("_homolog_ensembl_gene",
                 "_homolog_orthology_type",
                 "_homolog_associated_gene_name",
                 "_homolog_perc_id")

species <- c("ptroglodytes",
             "mmulatta",
             "cfamiliaris",
             "btaurus",
             "mmusculus",
             "rnorvegicus",
             "ggallus",
             "xtropicalis",
             "drerio",
             "dmelanogaster",
             "celegans")

queryAttr1 <- unlist(map(species[1:6], paste0, attrEndings))
queryAttr2 <- unlist(map(species[7:11], paste0, attrEndings))

queryAttr1 <- c("ensembl_gene_id", "description", "external_gene_name", queryAttr1)
queryAttr2 <- c("ensembl_gene_id", "description", "external_gene_name", queryAttr2)

query1 <- getBM(attributes = queryAttr1,
                filters = "ensembl_gene_id",
                values = lrpFam$ensembl_gene_id,
                mart = human)

query2 <- getBM(attributes = queryAttr2,
                filters = "ensembl_gene_id",
                values = lrpFam$ensembl_gene_id,
                mart = human)

results <- query1 %>%
    full_join(query2, by = c("ensembl_gene_id", "description", "external_gene_name")) %>%
    gather(attr, value, -ensembl_gene_id, -description, -external_gene_name) %>%
    group_by(ensembl_gene_id) %>%
    mutate(species = str_match(attr, "(^.*?)_")[ , 2],
           category = str_match(attr, "_(.*$)")[ , 2]) %>%
    dplyr::select(-attr) %>%
    spread(category, value)

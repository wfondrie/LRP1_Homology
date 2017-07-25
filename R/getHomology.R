library(tidyverse)
library(httr)
library(stringr)
library(forcats)
library(biomaRt)

# Set Marts
human <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
ptroglodytes <- useEnsembl("ensembl", dataset = "ptroglodytes_gene_ensembl")
mmulatta <- useEnsembl("ensembl", dataset = "mmulatta_gene_ensembl")
clupus <- useEnsembl("ensembl", dataset = "mmulatta_gene_ensembl")
btaurus <- useEnsembl("ensembl", dataset = "btaurus_gene_ensembl")
mmusculus <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
rnorvegicus <- useEnsembl("ensembl", dataset = "rnorvegicus_gene_ensembl")
ggallus <- useEnsembl("ensembl", dataset = "ggallus_gene_ensembl")
xtropicalis <- useEnsembl("ensembl", dataset = "xtropicalis_gene_ensembl")
drerio <- useEnsembl("ensembl", dataset = "drerio_gene_ensembl")
dmelanogaster <- useEnsembl("ensembl", dataset = "dmelanogaster_gene_ensembl")
celegans <- useEnsembl("ensembl", dataset = "celegans_gene_ensembl")


alignmentData <- GET("https://www.ncbi.nlm.nih.gov/gene/?Term=ortholog_gene_4292[group]")

searchData <- GET("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&Term=ortholog_gene_4292[group]")

dat <- content(searchData)

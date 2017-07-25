library(tidyverse)
library(httr)
library(stringr)
library(forcats)

queryBase <- c("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/")

searchHomologene <- function(IDs) {
  IDs <- paste(IDs, sep = "", collapse = ",")
  paste0(queryBase, "efetch.fcgi?db=homologene&id=", IDs, "&rettype=alignmentscores&retmode=text")
}

idList <- tibble(homologeneUID = c("56810",  #LRP1B
                                   "1744",   #LRP1,
                                   "20952",  #LRP2
                                   "55469",  #LDLR
                                   "443",    #VLDLR
                                   "31250",  #apoER2 (LRP8)
                                   "17964",  #LRP4
                                   "37248",  #A2M
                                   "37612",  #LRPAP1
                                   "30951",  #APOE
                                   "49153",   #FVIII
                                   "74303"   #PDGFB
                                   ),
                 gene = as.factor(c("LRP1B", "LRP1", "LRP2", "LDLR", "VLDLR", "LRP8", "LRP4",
                          "A2M", "LRPAP1", "APOE", "F8", "PDGFB")),
                 gene2 = as.factor(c("LRP1B", "LRP1", "LRP2", 
                          "LDLR", "VLDLR", "ApoER2", "LRP4", 
                          "A2M", "RAP", "ApoE", "FVIII", "PDGFB")),
                 group = as.factor(c(rep("LDLR Family", 7), rep("LRP1 Ligands", 5))))

alignmentData <- GET(searchHomologene(idList$homologeneUID))

dat <- content(alignmentData)
write_file(dat, "test.tsv")

b <- read_tsv("test.tsv")

# Formatting result ------------------------------------------------------------
orgs<- as.character(str_split(dat, ".: HomoloGene.+\\\n", simplify = T))
orgs2 <- orgs[orgs != ""]
orgs3 <- str_replace_all(orgs2, "\\\n", "")
orgs4 <- str_replace(orgs3, "^.*DNA    ", "")

orgs5 <- str_replace(orgs4, "H.sapiens.*?vs. ", "")
orgs6 <- str_replace(orgs5, "      [^v].*$", "")
orgs7 <- str_replace_all(orgs6, " Blast      vs. ", ";")
orgs8 <- str_replace(orgs7, " Blast$", "")

tab <- tibble(gene = str_match(orgs4, "H.sapiens (.+?)vs")[, 2],
              dat = orgs8) %>% 
  group_by(gene) %>% 
  do(strlist = as.character(str_split(.$dat, ";", simplify = T))) %>%
  unnest(strlist) %>%
  separate(strlist,into = c("Species", "Symbol", "Protein", "DNA"), sep = " ") %>%
  gather(Level, Identity, Protein, DNA) %>%
  mutate(gene = as.factor(gene),
         Species = as.factor(Species),
         Level = as.factor(Level),
         Identity = as.numeric(Identity)) %>%
  left_join(idList)

fullGrid <- expand.grid(gene = levels(tab$gene), 
                        Species = levels(tab$Species), 
                        Level = levels(tab$Level))
fullGrid <- fullGrid %>% left_join(unique(select(tab, gene, gene2, group)))

tab <- tab %>% full_join(fullGrid)

# heatmap ----------------------------------------------------------------------
speciesOrder <- c("P.troglodytes",
                  "M.mulatta",
                  "C.lupus",
                  "B.taurus",
                  "M.musculus",
                  "R.norvegicus",
                  "G.gallus",
                  "X.tropicalis",
                  "D.rerio",
                  "D.melanogaster",
                  "A.gambiae",
                  "C.elegans")

geneOrder <- c("LDLR",
               "VLDLR",
               "LRP8",
               "LRP4",
               "LRP1",
               "LRP1B",
               "LRP2",
               "LRPAP1", 
               "A2M", 
               "APOE", 
               "F8", 
               "PDGFB")

tab <- tab %>% 
  mutate(Species = fct_relevel(Species, speciesOrder),
         Species2 = as.factor(str_replace(Species, "\\.", ". ")),
         Species2 = fct_reorder(Species2, as.numeric(Species)),
         gene = fct_relevel(gene, geneOrder),
         gene2 = fct_reorder(gene2, as.numeric(gene)))
  
tab %>%
  filter(group == "LDLR Family") %>%
  ggplot(aes(x = gene2, y = fct_rev(Species2), fill = Identity)) +
  geom_tile(color = "white", size = 1) +
  facet_wrap(~ Level, ncol = 2, strip.position = "bottom") +
  scale_fill_gradient(name = "% Identity",low = "#d7191c", high = "#2c7bb6", 
                       na.value = "white", limits = c(0, 100),
                      guide = guide_colorbar()) +
  scale_x_discrete(position = "top") +
  scale_y_discrete() + 
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black"),
        axis.title = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 0, vjust = 0),
        strip.background = element_rect(color = "black"),
        panel.border = element_rect(color = "black", fill = NA),
        legend.key = element_rect(color = "black", fill = NA))

ggsave("results/lrpFamHomologyFigure.png", height = 3, width = 6)

lrp <- tab %>%
  filter(gene == "LRP1") %>%
  group_by(Species, gene, gene2, Level) %>%
  summarize(Identity = mean(Identity))

tab %>% 
  filter(group == "LRP1 Ligands") %>%
  group_by(Species, gene, gene2, Level) %>%
  summarize(Identity = mean(Identity)) %>%
  ggplot(aes(x = Species, y = Identity, color = gene2, group = gene2)) +
  geom_point(data = lrp, color = "black", aes(fill = gene2)) +
  geom_line(data = lrp, color = "black", aes(fill = gene2), show.legend = T) +
  scale_fill_discrete(name = NULL, guide = guide_legend(order = 1)) +
  scale_color_discrete(name = "Ligand", guide = guide_legend(order = 2)) +
  geom_point() +
  geom_line() +
  facet_grid(. ~ Level) +
  ylab("% Identity") +
  theme(panel.border = element_rect(color = "black", fill = NA),
        panel.background = element_blank(),
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.title.x = element_blank(),
        strip.background = element_rect(color = "black"),
        legend.key = element_rect(color = NA, fill = NA)) 

ggsave("results/lrpLigandFigure.png", width = 6, height = 3)

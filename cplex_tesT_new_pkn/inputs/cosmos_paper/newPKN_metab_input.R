library(readr)
library(stringr)
library('org.Hs.eg.db')
library(KEGGREST)
library(webchem)
setwd("~/Dropbox/COSMOS_MSB/")

#### Preparation des input metabolomic

ttop_tumour_vs_healthy <- as.data.frame(
  read_csv("data/metab_ttop_tumour_vs_healthy.csv"))
metab_to_kegg <- as.data.frame(
  read_csv("support/metab_to_kegg.txt"))
meta_network_with_X <- as.data.frame(
  read_csv("~/Dropbox/meta_PKN_redhuman/result/cosmos_PKN.csv"))

mapping <- keggConv("pubchem", source = metab_to_kegg$KEGG, querySize = 77)

sid_to_cid <- list()
for(i in 1:length(mapping))
{
  pubchem <- gsub("pubchem:","",mapping[i])
  print(pubchem)
  pubchem <- as.data.frame(get_cid(pubchem, from = "sid", domain = "substance"))
  sid_to_cid[[i]] <- pubchem
}
sid_to_cid_df <- as.data.frame(do.call(rbind, sid_to_cid))
sid_to_cid_df <- sid_to_cid_df[complete.cases(sid_to_cid_df),]
sid_to_cid_df_vec <- sid_to_cid_df$cid
names(sid_to_cid_df_vec) <- sid_to_cid_df$query

for(i in 1:length(mapping))
{
  pubchem_sid <- gsub("pubchem:","",mapping[i])
  if(pubchem_sid %in% names(sid_to_cid_df_vec))
  {
    mapping[i] <- paste("pubchem:",sid_to_cid_df_vec[pubchem_sid],sep = "")
  } else
  {
    mapping[i] <- paste("pubchem:sid",pubchem_sid,sep = "")
  }
}

kegg_to_pubchem <- data.frame(mapping)
kegg_to_pubchem$X1 <- gsub("cpd:","",row.names(kegg_to_pubchem))
kegg_to_pubchem$X2 <- gsub("pubchem:","",kegg_to_pubchem$mapping)
kegg_to_pubchem <- kegg_to_pubchem[,c(2,3)]

kegg_to_pubchem$X2 <- paste("XMetab__",kegg_to_pubchem$X2, sep = "")

compartment_codes <- unique(c(meta_network_with_X$source,meta_network_with_X$target))
compartment_codes <- compartment_codes[grepl("Metab",compartment_codes)]
compartment_codes <- unique(str_match(compartment_codes,"___.____"))

compartment_mapping <- list()
for(i in 1:length(compartment_codes))
{
  df <- kegg_to_pubchem
  df$X2 <- paste(df$X2, compartment_codes[i], sep = "")
  compartment_mapping[[i]] <- df
}

compartment_mapping <- as.data.frame(do.call(rbind, compartment_mapping))

compartment_mapping <- compartment_mapping[
  compartment_mapping$X2 %in% meta_network_with_X$source |
    compartment_mapping$X2 %in% meta_network_with_X$target,
]

kegg_to_pubchem_with_comp <- compartment_mapping
names(kegg_to_pubchem_with_comp) <- c("KEGG","pubchem")

full_mapping <- merge(metab_to_kegg, kegg_to_pubchem_with_comp, by = "KEGG")

names(ttop_tumour_vs_healthy)[1] <- "metab_name"

ttop_tumour_vs_healthy$metab_name <- tolower(ttop_tumour_vs_healthy$metab_name)
full_mapping$metab_name <- tolower(full_mapping$metab_name)

ttop_tumour_vs_healthy <- merge(ttop_tumour_vs_healthy, full_mapping, by = "metab_name")
ttop_tumour_vs_healthy <- ttop_tumour_vs_healthy[,c(9,2:7)]

ttop_tumour_vs_healthy <- ttop_tumour_vs_healthy[ttop_tumour_vs_healthy$P.Value <= 0.1,]

carnival_metab_input <- as.data.frame(t(ttop_tumour_vs_healthy$t))
names(carnival_metab_input) <- ttop_tumour_vs_healthy$pubchem

write_csv(carnival_metab_input, file = "../cosmos_tests/inputs/cosmos_paper/metab_input_COSMOS_newPKN.csv")                                      

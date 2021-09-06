library(OmnipathR)
library(org.Hs.eg.db)
library(readr)

full_pkn <- as.data.frame(import_all_interactions())

clean_PKN <- full_pkn[full_pkn$consensus_stimulation == 1 | full_pkn$consensus_inhibition == 1,]
clean_PKN <- clean_PKN[clean_PKN$dorothea_level %in% c("A","A;B","A;C","A;D","B","B;D","C","C;D") | is.na(clean_PKN$dorothea_level),]

clean_PKN$sign <- clean_PKN$is_stimulation - clean_PKN$consensus_inhibition

clean_PKN <- clean_PKN[,c(3,4,19)]

clean_PKN_supp <- clean_PKN[clean_PKN$sign == 0,]
clean_PKN_supp$sign <- -1
clean_PKN[clean_PKN$sign == 0,"sign"] <- 1

clean_PKN <- as.data.frame(rbind(clean_PKN, clean_PKN_supp))

prots <- unique(c(clean_PKN$source_genesymbol, clean_PKN$target_genesymbol))

mapping_symbole_to_entrez <- mapIds(org.Hs.eg.db, prots, 'ENTREZID', 'SYMBOL')
head(mapping_symbole_to_entrez)

for(i in 1:2)
{
  clean_PKN[,i] <- sapply(clean_PKN[,i], function(x, mapping_symbole_to_entrez){
    x <- mapping_symbole_to_entrez[x]
  }, mapping_symbole_to_entrez = mapping_symbole_to_entrez)
}

for(i in 1:2)
{
  clean_PKN[,i] <- paste("X",clean_PKN[,i],sep = "")
}

names(clean_PKN) <- c("source","target","sign")
write_csv(x = clean_PKN, file = "Dropbox/meta_PKN_redhuman/result/clean_omnipath_PKN.csv")

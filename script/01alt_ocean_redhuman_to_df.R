library(KEGGREST)
library(webchem)
library(stringr)
library(ocean)
library('org.Hs.eg.db')
library(readr)


all_pathways <- unique(recon2_redhuman$pathway)
sub_network <- model_to_pathway_sif(pathway_to_keep = all_pathways$X1)

sub_network_nocofact <- remove_cofactors(sub_network)

sub_network_nocofact <- compress_transporters(sub_network_nocofact = sub_network_nocofact)

sub_network_nocofact <- split_transaminases(sub_network_nocofact = sub_network_nocofact)

sub_network_nocofact <- sub_network_nocofact$reaction_network

metabs <- unique(c(sub_network_nocofact$source, sub_network_nocofact$target))
metabs <- metabs[grepl("cpd:",metabs)]
metabs <- gsub("cpd:","",metabs)
metabs <- gsub("_.*","",metabs)
metabs <- unique(metabs)

mapping_p1 <- keggConv("pubchem", source = metabs[1:99], querySize = 99)
mapping_p2 <- keggConv("pubchem", source = metabs[100:190], querySize = 90)

mapping <- c(mapping_p1,mapping_p2)

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

row.names(sub_network_nocofact) <- 1:length(sub_network_nocofact[,1])

for(i in 1:length(sub_network_nocofact[,1]))
{
  for(j in 1:length(sub_network_nocofact[1,]))
  {
    if(gsub("_.*","",sub_network_nocofact[i,j]) %in% names(mapping))
    {
      sub_network_nocofact[i,j] <- paste(gsub("pubchem:","XMetab__",mapping[gsub("_.*","",sub_network_nocofact[i,j])]), gsub(".*_","",sub_network_nocofact[i,j]), sep = "___")
      sub_network_nocofact[i,j] <- paste(sub_network_nocofact[i,j],"____",sep = "")
    }
  }
}

sub_network_nocofact$source <- gsub("[.]",">",sub_network_nocofact$source)
sub_network_nocofact$target <- gsub("[.]",">",sub_network_nocofact$target)

enzymes <- unique(c(sub_network_nocofact$source, sub_network_nocofact$target))
enzymes <- enzymes[!grepl("XMetab",enzymes)]

mapping_symbole_to_entrez <- mapIds(org.Hs.eg.db, enzymes, 'ENTREZID', 'SYMBOL')
mapping_symbole_to_entrez <- mapping_symbole_to_entrez[!is.na(mapping_symbole_to_entrez)]

###Build the connectors to omnipath

gene_connectors_list <- list()
mega_complexe <- list()
k <- 1
for(i in 1:length(sub_network_nocofact[,1]))
{
  print(i)
  for(j in 1:length(sub_network_nocofact[1,]))
  {
    element <- sub_network_nocofact[i,j]
    if(grepl("XMetab",element))
    {
      sub_network_nocofact[i,j] <- element
    } else
    {
      if(grepl("_[a-z]$",element) & !(grepl("^EX",element)))
      {
        element <- paste("XMetab",element,sep = "")
        sub_network_nocofact[i,j] <- element
      } else
      {
        if(grepl("_reverse",element))
        {
          element <- gsub("_reverse","",element)
          if(grepl("_",element))
          {
            elements <- str_split(string = element, pattern = "_")[[1]]
            if(length(elements) < 10) #if megacomplexe, we remove cause it's to hard to process
            {
              gene_connector_list <- list()
              for(gene in elements)
              {
                gene_connector_list[[gene]] <- c(paste("X",gene,sep=""),paste("XGene__",paste(element,"_reverse",sep=""),sep=""))
              }
              gene_connectors <- as.data.frame(do.call(rbind,gene_connector_list))
              names(gene_connectors) <- c("source","target")
              gene_connectors_list[[k]] <- gene_connectors
              
              sub_network_nocofact[i,j] <- paste("XGene__",paste(element,"_reverse",sep=""),sep="")
            } else
            {
              mega_complexe[[i]] <- i
            }

          } else
          {
            if(gsub(">.*","",element) %in% names(mapping_symbole_to_entrez))
            {
              if(grepl(">",element))
              {
                element <- gsub(".*>",paste(mapping_symbole_to_entrez[gsub(">.*","",element)],">",sep=""),element)
              } else
              {
                element <- mapping_symbole_to_entrez[element]
              }
            }
            clean_elem <- element
            element <- paste("XGene__",paste(element,"_reverse",sep=""),sep="")
            gene_connectors_list[[k]] <- c(paste("X",gsub(">.*","",clean_elem),sep=""),element)
            
            sub_network_nocofact[i,j] <- element
          }
        } else
        {
          if(grepl("_",element))
          {
            elements <- str_split(string = element, pattern = "_")[[1]]
            if(length(elements) < 10)
            {
              gene_connector_list <- list()
              for(gene in elements)
              {
                print(gene)
                gene_connector_list[[gene]] <- c(paste("X",gene,sep=""),paste("XGene__",element,sep=""))
              }
              gene_connectors <- as.data.frame(do.call(rbind,gene_connector_list))
              names(gene_connectors) <- c("source","target")
              gene_connectors_list[[k]] <- gene_connectors
              
              sub_network_nocofact[i,j] <- paste("XGene__",element,sep="")
            } else
            {
              mega_complexe[[i]] <- i
            }
            
          } else
          {
            if(gsub(">.*","",element) %in% names(mapping_symbole_to_entrez))
            {
              if(grepl(">",element))
              {
                element <- gsub(".*>",paste(mapping_symbole_to_entrez[gsub(">.*","",element)],">",sep=""),element)
              } else
              {
                element <- mapping_symbole_to_entrez[element]
              }
              
            }
            clean_elem <- element
            element <- paste("XGene__",element,sep="")
            gene_connectors_list[[k]] <- c(paste("X",gsub(">.*","",clean_elem),sep=""),element)
            
            sub_network_nocofact[i,j] <- element
          }
        }
      }
    }
  k <- k+1
  }
}
gene_connectors <- as.data.frame(do.call(rbind,gene_connectors_list))
gene_connectors <- unique(gene_connectors)
gene_connectors <- gene_connectors[!grepl("X[A-Za-z]",gene_connectors[,1]),]

mega_complexe <- unlist(mega_complexe)

sub_network_nocofact <- sub_network_nocofact[-mega_complexe,]

sub_network_nocofact <- as.data.frame(rbind(sub_network_nocofact,gene_connectors))
sub_network_nocofact$source <- gsub("cpd:","",sub_network_nocofact$source)
sub_network_nocofact$target <- gsub("cpd:","",sub_network_nocofact$target)

#looks good at this stage
sub_network_nocofact$sign <- 1
row.names(sub_network_nocofact) <- 1:length(sub_network_nocofact$source)

#there are some > left in the connectors part, will remove them manually with a magic number
sub_network_nocofact[2044:2753,1] <- gsub(">","",sub_network_nocofact[2044:2753,1])

sub_network_nocofact$source <- gsub("XGene__","XGene0__",sub_network_nocofact$source)
sub_network_nocofact$target <- gsub("XGene__","XGene0__",sub_network_nocofact$target)

write_csv(sub_network_nocofact, file = "Dropbox/meta_PKN_redhuman/result/redhuman_network.csv")

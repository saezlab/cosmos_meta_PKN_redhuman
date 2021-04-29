library(readr)

setwd("Dropbox/meta_PKN_redhuman/")

clean_omnipath_PKN <- as.data.frame(read_csv("result/clean_omnipath_PKN.csv"))

redhuman_network <- as.data.frame(read_csv("result/redhuman_network.csv"))

connectors_to_keep <- apply(redhuman_network, 1, function(x,clean_omnipath_PKN){
  if(sum(grepl("Metab",x)) == 0)
  {
    test <- x[1] %in% clean_omnipath_PKN$target
    return(test)
  } else
  {
    return(TRUE)
  }
},clean_omnipath_PKN = clean_omnipath_PKN)

redhuman_network <- redhuman_network[connectors_to_keep,]

cosmos_PKN <- as.data.frame(rbind(clean_omnipath_PKN,redhuman_network))

STITCH_900 <- as.data.frame(read_csv("result/STITCH_900.csv"))

STITCH_900 <- STITCH_900[(STITCH_900$target %in% cosmos_PKN$source | STITCH_900$target %in% cosmos_PKN$target),]

# STITCH_900 <- STITCH_900[(STITCH_900$source %in% cosmos_PKN$source | STITCH_900$source %in% cosmos_PKN$target),]

cosmos_PKN <- as.data.frame(rbind(cosmos_PKN,STITCH_900))

cosmos_PKN <- cosmos_PKN[,c(1,3,2)]
names(cosmos_PKN)[2] <- "interaction"

cosmos_PKN <- unique(cosmos_PKN)

write_csv(cosmos_PKN, "result/cosmos_PKN.csv")

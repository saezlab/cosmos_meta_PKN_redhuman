library(cosmosR)
library(readr)

data(metabolite_to_pubchem)
data(omnipath_ptm)

setwd(dir = "~/Dropbox/meta_PKN_redhuman/cplex_tesT_new_pkn/")


#In order to adapt options to users specification we can load them into a variable 
#that will then be passed to preprocess_COSMOS_signaling_to_metabolism CARNIVAL_options parameter
my_options <- default_CARNIVAL_options()

#Here the user should provide a path to its CPLEX executable (only cplex at the moment, other solvers will be documented soon !)
my_options$solverPath <- "~/Documents/cplex" #or cbc solver executable
my_options$solver <- "cplex" #or cbc
my_options$timelimit <- 1800
my_options$mipGAP <- 0.05
my_options$threads <- 6

signaling_input_COSMOS <- as.data.frame(read_csv("inputs/cosmos_paper/signaling_input_COSMOS.csv"))
signaling_input_COSMOS_vec <- as.numeric(signaling_input_COSMOS[1,])
names(signaling_input_COSMOS_vec) <- names(signaling_input_COSMOS)

metab_input_COSMOS <- as.data.frame(read_csv("inputs/cosmos_paper/metab_input_COSMOS_newPKN.csv"))
metab_input_COSMOS_vec <- as.numeric(metab_input_COSMOS[1,])
names(metab_input_COSMOS_vec) <- names(metab_input_COSMOS)

ttop_rna <- as.data.frame(read_csv("inputs/cosmos_paper/RNA_ttop_tumorvshealthy.csv"))
RNA_input <- ttop_rna[,"t"]
names(RNA_input) <- paste0("X",ttop_rna$ID)

meta_network <- as.data.frame(read_csv("inputs/cosmos_paper/cosmos_PKN.csv"))

lost <-  metab_input_COSMOS_vec[!(names(metab_input_COSMOS_vec) %in% meta_network$source |
                                                           names(metab_input_COSMOS_vec) %in% meta_network$target)]


metab_input_COSMOS_vec <- metab_input_COSMOS_vec[names(metab_input_COSMOS_vec) %in% meta_network$source |
                                                   names(metab_input_COSMOS_vec) %in% meta_network$target]

test_for <- preprocess_COSMOS_signaling_to_metabolism(meta_network = meta_network,
                                                      signaling_data = signaling_input_COSMOS_vec,
                                                      metabolic_data = metab_input_COSMOS_vec,
                                                      diff_expression_data = RNA_input,
                                                      maximum_network_depth = 8,
                                                      remove_unexpressed_nodes = T,
                                                      CARNIVAL_options = my_options
                                                      
)

my_options$timelimit <- 7200

test_result_for <- run_COSMOS_signaling_to_metabolism(data = test_for,
                                                      CARNIVAL_options = my_options)
save(test_result_for, file = "cplex_res_formard.RData")

metabolite_to_pubchem$name <- gsub(",","_",metabolite_to_pubchem$name)
metabs <- metabolite_to_pubchem$pubchem
metabolite_to_pubchem <- metabolite_to_pubchem$name
names(metabolite_to_pubchem) <- metabs
metabolite_to_pubchem["65359"] <- "GSSG"
metabolite_to_pubchem["439756"] <- "O-acetylcarnitine"
metabolite_to_pubchem["85"] <- "Carnitine"
metabolite_to_pubchem["440641"] <- "F6P"
metabolite_to_pubchem["165007"] <- "Sedoheptulose_7-phosphate"
metabolite_to_pubchem["34755"] <- "S-adenosyl-L-methionine"

test_result_for <- format_COSMOS_res(test_result_for,
                                     metab_mapping = metabolite_to_pubchem,
                                     measured_nodes = unique(c(names(metab_input_COSMOS_vec),
                                                               names(signaling_input_COSMOS_vec))),
                                     omnipath_ptm = omnipath_ptm)

#View(test_result_for[[1]]) #SIF
#View(test_result_for[[2]]) #ATTRIBUTES
write_csv(test_result_for[[1]], file = "forward_sif.csv")
write_csv(test_result_for[[2]], file = "forward_att.csv")

#### BACKWARD run of COSMOS, to connect metabolism to signaling
my_options$timelimit <- 1800

test_back <- preprocess_COSMOS_metabolism_to_signaling(meta_network = meta_network,
                                                       signaling_data = signaling_input_COSMOS_vec,
                                                       metabolic_data = metab_input_COSMOS_vec,
                                                       diff_expression_data = RNA_input,
                                                       maximum_network_depth = 8,
                                                       remove_unexpressed_nodes = T,
                                                       CARNIVAL_options = my_options
                                                       
)

my_options$timelimit <- 7200

test_result_back <- run_COSMOS_metabolism_to_signaling(data = test_back,
                                                       CARNIVAL_options = my_options)

save(test_result_back, file = "cplex_res_backward.RData")

test_result_back <- format_COSMOS_res(test_result_back,
                                      metab_mapping = metabolite_to_pubchem,
                                      measured_nodes = unique(c(names(metab_input_COSMOS_vec),
                                                                names(signaling_input_COSMOS_vec))),
                                      omnipath_ptm = omnipath_ptm)

#View(test_result_back[[1]]) #SIF
#View(test_result_back[[2]]) #ATTRIBUTES

###Merge forward and backward networks

full_sif <- as.data.frame(rbind(test_result_for[[1]], test_result_back[[1]]))
full_attributes <- as.data.frame(rbind(test_result_for[[2]], test_result_back[[2]]))

full_sif <- unique(full_sif)
full_attributes <- unique(full_attributes)

save.image("cplex_res_translated.RData")

write_csv(full_sif, file = "Dropbox/meta_PKN_redhuman/cplex_tesT_new_pkn/full_sif.csv")
write_csv(full_attributes, file = "Dropbox/meta_PKN_redhuman/cplex_tesT_new_pkn/full_att.csv")

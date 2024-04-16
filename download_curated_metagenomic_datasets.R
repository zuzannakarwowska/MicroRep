library(curatedMetagenomicData)
library(curatedMetagenomicAnalyses)
library(dplyr)
library(doParallel)
library(foreach)

cores <- 4
cl <- makeCluster(cores) 
registerDoParallel(cl)
clusterCall(cl, function(x) .libPaths(x), .libPaths())
pckgs <- c('curatedMetagenomicData', 'curatedMetagenomicAnalyses', 'doParallel', 'foreach', 'dplyr') 


data('sampleMetadata')

out_dir = '/Users/zkarwowska/Desktop/microbiome_gpt/datasets/controls/'

metadata <- sampleMetadata %>%
  dplyr::filter(body_site == 'stool') %>%
  dplyr::filter(is.na(study_condition) == FALSE) %>%
  dplyr::filter(study_condition == 'control')

write.csv(metadata, paste0('/Users/zkarwowska/Desktop/microbiome_gpt/datasets/controls/metadata_healthy.csv'))

studies <- unique(metadata$study_name)

x <- foreach(study=studies, .packages = pckgs, .verbose = TRUE) %dopar% {
  
  dataset_metadata <- metadata %>%
    dplyr::filter(study_name == study)

  metagenomics_object_pathways <- suppressMessages(dataset_metadata %>% returnSamples('pathway_abundance', counts = TRUE))
  metagenomics_object_taxonomy <- suppressMessages(dataset_metadata %>% returnSamples('relative_abundance', counts = TRUE))
  
  # Retrieve the count data
  pathways_countdata <- assays(metagenomics_object_pathways)[[1]] %>%
    t() %>%
    as.data.frame()
  taxonomy_countdata <- assays(metagenomics_object_taxonomy)[[1]] %>%
    t() %>%
    as.data.frame()
  #write to csv

  write.csv(pathways_countdata, paste0(out_dir, study, '_pathways.csv'))
  write.csv(taxonomy_countdata, paste0(out_dir, study, '_taxonomy.csv'))
} 

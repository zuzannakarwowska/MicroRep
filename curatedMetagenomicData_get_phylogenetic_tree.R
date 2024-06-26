library(curatedMetagenomicData)
library('curatedMetagenomicAnalyses')
library('dplyr')
library(phyloseq)
library(TreeSummarizedExperiment)
library(ape)

tse <- curatedMetagenomicData("AsnicarF_2017.relative_abundance", dryrun = FALSE, rownames = "short") |> mergeData()
phylo_tree <- rowTree(tse)
write.tree(phylo_tree, file="/Users/zkarwowska/Desktop/microbiome_gpt/unifrac/phylogenetic_tree.nwk")

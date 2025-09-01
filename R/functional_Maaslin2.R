library(tidyverse)
library(phyloseq)
library(readxl)
library(G4Micro)


MODE = "LOCAL"

if(exists("input")){
  opt <- input
}else{
  if(MODE == "IATA"){
  }else{
    CODEDIR = "/home/carlos/Documentos/CORALS/scripts_PAR/240806scripts/depression_scripts/"
    opt <- list(fc=1,
                pval=0.05,
                num_genes_default=5,
                ptype="adjusted",
                fctype="shrunk",
                mincount= 1,
                minsampleswithcount= 6,
                out= "/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1/functional_food1/", #"/home/carlos/projects/Psoriasis_20231005/results2/",
                input=  "/home/carlos/Escritorio/202311_DEPRESION/mg09_combinempa",
                modules_path=  "/home/carlos/Escritorio/202311_DEPRESION/202311_DEPRESION/funcional1/kegg_modules/",
                cazy="/home/carlos/Escritorio/202311_DEPRESION/202311_DEPRESION/cazy/",
                #r_functions= "/home/carlos/Escritorio/202311_DEPRESION/depression_scripts/metagenomics_core_functions.R", #"/home/carlos/projects/Psoriasis_20231005/Psoriasis/Analysis_code/metagenomics_core_functions.R",
                #r_functions2= "/home/carlos/Escritorio/202311_DEPRESION/depression_scripts/metagenomics_predictive_functions.R",
                #r_functions3= "/home/carlos/Escritorio/202311_DEPRESION/depression_scripts/sankey_functions.R",
                #r_functions4= "/home/carlos/Escritorio/202311_DEPRESION/depression_scripts/functional_auxiliary_functions.R",
                input_funcional = "/home/carlos/Documentos/CORALS/results_cluster_240924/Humann3_analisis_funcional/MERGED_changed_names/",
                #dearesult = "/home/carlos/Escritorio/202311_DEPRESION/results_rstudio_v2_4/DeSEQ2/remove_tanda2/DESEQ2_all_results_remove_tanda2.R",
                #metadata =  "/home/carlos/Escritorio/202311_DEPRESION/metadatos_MC_AL12042023_CMcopy.xlsx",
                phobj_all ="/home/carlos/Documentos/CORALS/results_rstudio/results_Agosto25_1//foodPCA/phyloseq_list_foodPCA_withNMF_withInc.RData"
    )
  }
}

dists <- c("bray", "jaccard", "unifrac", "wunifrac")
dtitles <- c("Bray-Curtis distances", "Jaccard distances", "UniFrac distances", "Weighted UniFrac distances")



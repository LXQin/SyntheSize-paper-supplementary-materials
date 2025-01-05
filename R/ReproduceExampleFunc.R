# Source code containing all functions in DataOrganization.Rmd, Evaluation.Rmd, Visualization.Rmd. 
# Just for ReproduceExample.Rmd. 

###########################################
## Functions from DataOrganization.Rmd ####
###########################################

library(TCGAbiolinks)
library(tidyverse)
library(DANA)
library(SummarizedExperiment)

# clean raw data for any cancer type, ie. LAML and SKCM 
get_count_df <- function(cancer){
  # This function clean the TCGA-cancer.csv file in ../RealData/ folder
  # @param: cancer - character specifying cancer name, eg. "SKCM"
  # @return: list contains df_count - count table with patient (not necessary) and Sample (unique) columns
  #          mat_count - count table 
  
  cancer_ori <- read.csv(paste("../RealData/TCGA-",cancer,".csv",sep=""))
  df_sub <- select(cancer_ori, -contains("million"))
  df_sub <- select(df_sub, -contains("mapped"))
  colnames(df_sub) <-  sub("read_count_", "", colnames(df_sub))
  tryCatch({df_sub <- select(df_sub,-"X")}, 
           error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
  df_count <- as.data.frame(t(select_if(df_sub, is.numeric)))
  colnames(df_count) <- df_sub$miRNA_ID
  # patients may be the same, but related to different sample id.
  df_count$patient  <- sapply(rownames(df_count),function(x){
    paste(strsplit(x,"\\.")[[1]][1:3],collapse="-")})
  df_count$Sample <- rownames(df_count)
  return(list(df_count = df_count, mat_count = Filter(is.numeric, df_count)))
}

# combine two cancer types and add cancer type label groups
get_combine_count_df <- function(df_count_1, df_count_2, group_name){
  # This function combines two cancer count into one and add type labels
  # @param: df_count_1 -  count table with patient and Sample for cancer type 1 
  # @param: df_count_2 -  count table with patient and Sample for cancer type 2 
  # @param: group_name - vector contains the two levels of group
  # @return: count_list - a list including df_count and mat_count for combined two cancer types
  
  features <- colnames(df_count_1)
  df_count <- as.data.frame(rbind(df_count_1, df_count_2[,features]))
  df_count <- Filter(is.numeric, df_count)
  df_count$groups <- rep(group_name, c(nrow(df_count_1), nrow(df_count_2)))
  mat_count <- select_if(df_count, is.numeric)
  return(list(df_count = df_count, mat_count = mat_count))
}

# get positive control datasets for any count table
get_count_positive_df <- function(count_list, thres, poly = F, coords = NULL, plot = T){
  # This function filter positive control genes, positive means mean(log2(count+1)) for this gene is larger than some threshold
  # @param: count_list - a list including df_count and mat_count, result from get_count_df()
  # @param: thres - threshold for positive control genes, in log2 scale
  # @param: poly - whether to use polycistronic clusters in defining the positive control features
  # @param: coords - coordinates data frame based on miRBase miRNA definitions, must provide when poly is TRUE
  # @return: a list containing df_count and mat_count but with positive control genes as columns.
  df_count <- count_list$df_count
  mat_count <- count_list$mat_count
  thres_plot <- thres
  
  gene_means <- apply(mat_count, 2, function(x){mean(log2(x+1))})
  gene_stds <- apply(mat_count, 2, function(x){sd(log2(x+1))})
  
  if(poly){
    clusters <- DANA::defineClusters(genes = coords$name, chr = coords$chr, 
                                     pos = (coords$start + coords$end)/2)
    features_pos <- DANA::defineControls(raw = t(mat_count), 
                                         tZero = 0, tPoor = thres-1, tWell = thres,
                                         clusters = clusters)$posControls
    thres_plot <- log2(thres+1)
  }
  else{
    features_pos <- colnames(mat_count)[gene_means >= thres]
  }
  
  mat_count_pos <- mat_count[, features_pos]
  df_count_pos <- data.frame(mat_count_pos, select_if(df_count, ~!is.numeric(.x)))
  
  if(plot){
    p <- ggplot(data.frame(gene_means, gene_stds))+
      geom_point(aes(x = gene_means, y = gene_stds, color = gene_means), alpha = .4, show.legend = F)+
      scale_color_gradient(low = "#0091ff", high = "#f0650e")+
      geom_vline(xintercept = thres_plot, color = "red")+
      labs(x = "mean(log2(count+1))", y = "sd(log2(count+1))",
           title = paste("# of positive genes = ", length(features_pos)))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    print(p)
  }
  return(list(mat_count = mat_count_pos, df_count = df_count_pos))
}

get_count_positive_df_RNA <- function(cancer, thres_mean = 5, thres_sd = 2, plot = T){
  ## This function get positive control genes from either TCGA-RNAcancer.csv or dataframe
  # @param: cancer - either the name of cancer, e.g RNABRCA, or the dataframe with genes in rows.
  # @param: thres_mean - filtering out genes with log2 mean lower than thres_mean
  # @param: thres_sd - filtering out genes with log2 sd lower than thres_sd
  # @param: plot -  logical, whether plot the log2 mean vs log2 sd
  # @output: a dataframe with positive control genes in columns and samples in rows
  if(!is.character(cancer)){
    mat_count <- t(cancer)
    cat(paste("Input data directly with dimension", ncol(cancer), " and sample size", nrow(cancer)))
  }
  else{
    cancer_ori <- read.csv(paste("../RealData/TCGA-",cancer,".csv",sep=""), row.names = 1)
    mat_count <- t(cancer_ori)
  }
  
  
  gene_means <- apply(mat_count, 2, function(x){mean(log2(x+1))})
  gene_stds <- apply(mat_count, 2, function(x){sd(log2(x+1))})
  features_pos <- colnames(mat_count)[gene_means > thres_mean & gene_stds > thres_sd]
  mat_count_pos <- mat_count[, features_pos]
  mat_count_pos <- as.data.frame(mat_count_pos)
  mat_count_pos$samples <- rownames(mat_count)
  if(plot){
    p <- ggplot(data.frame(gene_means, gene_stds))+
      geom_point(aes(x = gene_means, y = gene_stds, color = gene_means), alpha = .4, show.legend = F)+
      scale_color_gradient(low = "#0091ff", high = "#f0650e")+
      geom_vline(xintercept = thres_mean, color = "red")+
      geom_hline(yintercept = thres_sd, color = "red")
    labs(x = "mean(log2(count+1))", y = "sd(log2(count+1))",
         title = paste("# of positive genes = ", length(features_pos)))+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))
    print(p)
  }
  return(mat_count_pos)
}

# check normalization performance of dataset and get normalized datasets
check_norm <- function(dataname, selected = c("TMM", "DESeq", "PoissonSeq", "RUVr"), save = F){
  # This function checks normalization of dataname and save the normalized counts as new csv files.
  # @param: dataname - characters indicating dataname for csv file
  # @param: selected - normalization names
  # @param: save - whether or not to save the normalized datasets
  # @return: a list of normalized counts
  dat <- read.csv(paste("../RealData/", dataname, ".csv", sep=""))
  count <- t(select_if(dat, is.numeric))
  genes <- unlist(lapply(strsplit(rownames(count),"\\."), function(x){paste(x,collapse = "-")}))
  samples <- dat$samples
  colnames(count) <- samples
  rownames(count) <- genes
  if("groups" %in% colnames(dat)){groups <- dat$groups}
  else{groups <- rep("Single", ncol(count))}
  # coords <- read.csv("../RealData/miRBase_coords.csv") 
  # clusters <- DANA::defineClusters(genes = coords$name, chr = coords$chr, 
  #                                  pos = (coords$start + coords$end)/2)
  # controls <- DANA::defineControls(raw = count, 
  #                                  tZero = 2, tPoor = 50, tWell = 64,
  #                                  clusters = clusters)
  count_normed <- DANA::applyNormalization(data = count, groups = groups, method = selected)
  # assess <- DANA::assessNormalization(raw = count, normalized = count_normed, negControls = controls$negControls, posControls = controls$posControls, clusters = clusters)
  # DANA::plotDANA(assess)
  sel_normed <- count_normed[selected]
  for (dat_ind in 1:length(selected)) {
    norm_name <- selected[dat_ind]
    if(length(unique(groups)) > 1){
      sel_normed[[dat_ind]] <- data.frame(t(sel_normed[[dat_ind]]), samples = samples, groups = groups)
    }
    else{sel_normed[[dat_ind]] <- data.frame(t(sel_normed[[dat_ind]]), samples = samples)}
    if(save){
      write.table(sel_normed[[dat_ind]], 
                  paste("../RealData/", dataname, "_", selected[dat_ind], ".csv", sep = ""), 
                  sep =",", row.names = F, col.names = T)
    }
  }
  return(sel_normed)
}


###########################################
## Functions from Evaluation.Rmd ####
###########################################
library(ggplot2)
library(Rtsne)
library(tidyverse)
library(cowplot)
library(aricode)
library(RColorBrewer)
library(DANA)
library(dgof)
library(precision.seq)
library(epiR)
library(umap)
mycolors <- c(brewer.pal(n = 8, name = "Blues")[c(2,4,6)], brewer.pal(n = 8, name = "Oranges")[c(2,4,6)], brewer.pal(n = 8, name = "Greens")[c(2,4,6)])

listClusters <- function(geneClusters) {
  geneClusters <- as.factor(geneClusters)
  clusters <- list()
  clusters <- lapply(
    levels(geneClusters),
    function(x) names(geneClusters)[geneClusters==x]
  )
  names(clusters) <- levels(geneClusters)
  return(clusters)
}
compute.cc <- function(rawCor, normCor, clusters) {
  # remove genes that are not present in both models
  rawCor <- rawCor[!is.na(match(colnames(rawCor), colnames(normCor))),
                   !is.na(match(colnames(rawCor), colnames(normCor)))]
  normCor <- normCor[!is.na(match(colnames(normCor), colnames(rawCor))),
                     !is.na(match(colnames(normCor), colnames(rawCor)))]
  
  # consider only cluster with multiple genes
  clusters <- clusters[lengths(clusters) > 1]
  clusters.raw <- c()
  clusters.norm <- c()
  for (clust in clusters) {
    # only consider cluster genes that are positive controls
    clust.genes <- clust[stats::na.omit(match(colnames(rawCor), clust))]
    if (length(clust.genes) < 2) {
      next  # disregard clusters with less than 2 genes
    }
    # subset of correlations in the cluster
    rawCor.clust <- rawCor[clust.genes, clust.genes]
    clusters.raw <- c(clusters.raw, rawCor.clust[upper.tri(rawCor.clust)])
    normCor.clust <- normCor[clust.genes, clust.genes]
    clusters.norm <- c(clusters.norm, normCor.clust[upper.tri(normCor.clust)])
  }
  
  ## Compute the concordance correlation between nonzero partial correlations
  # idx <- clusters.raw | clusters.norm
  # if (sum(idx)>0) {
  #   cc <- DescTools::CCC(as.vector(clusters.raw[idx]),
  #                        as.vector(clusters.norm[idx]))$rho.c$est
  # } else {
  #   cc <- 0
  # }
  
  cc <- DescTools::CCC(as.vector(clusters.raw),as.vector(clusters.norm))$rho.c$est
  
  return(cc)
}


get_bind_list <- function(real, generated, real_log = F, generated_log = T){
  
  # This function bind the real samples with the generated samples.
  # @param: real - csv path of real dataset, N * p with header gene names and samples or groups
  # @param: generated -  path of generated dataset, logged, no header no row names
  # @param: real_log - whether real data are already logged, usually F
  # @param: generated_log - whether generated data are already logged, usually T
  # @return: bind_list -lists of 25 replications of combined datasets, each is 2N * p, first N are real samples, second N are generated samples, all samples are logged
  # @return: groups - vector of groups information, if there is no group in real data, groups = NULL, if there is groups, groups is a (1 + 25) * N vector, is the group info for (real, bind_list)
  # @return: samples - sample names of real data
  repli <- 5
  draw <- 5
  
  real <- read.csv(real, header = T)
  sample_names <- real$sample
  if("groups" %in% colnames(real)){groups <- real$groups}
  else{groups <- NULL}
  draw_paths <- paste("../GeneratedData/", generated, "_Draw", 1:draw, ".csv", sep = "")
  generated <- read.csv(draw_paths[1], header = F)
  for(draw_ind in 2:draw) {
    generated_draw <- read.csv(draw_paths[draw_ind], header = F)
    generated <- rbind(generated, generated_draw)
  }
  if("groups" %in% colnames(real)){
    g_default <- groups[1]
    g_other <- unique(groups)[unique(groups) != g_default]
    groups <- c(groups, ifelse(generated[,ncol(generated)] == 0, g_default, g_other))
    generated <- generated[, 1:(ncol(generated)-1)]
  }
  
  real <- select_if(real, is.numeric)
  feature_names <- unlist(lapply(strsplit(colnames(real),"\\."), 
                                 function(x){paste(x, collapse = "-")}))
  
  if(!real_log){real <- log2(real+1)}
  if(!generated_log){generated <- log2(generated+1)}
  stopifnot((nrow(real) * (repli * draw) == nrow(generated)) & (ncol(real) == ncol(generated)))
  
  colnames(generated) <- feature_names
  colnames(real) <- feature_names
  rownames(real) <- sample_names
  
  N <- nrow(real)
  bind_list <- list()
  for (repe in 1:(repli * draw)) {
    bind_list[[repe]] <- rbind(real, generated[((repe-1) * N + 1) : (repe * N), ])
    colnames(bind_list[[repe]]) <- feature_names
  }
  rm(generated)
  return(list(bind_list = bind_list, groups = groups, samples = sample_names))
}

heatmap_eval <-  function(dat_combine, main, log = TRUE){
  
  # This function plot the heatmap of data matrix dat_combine
  # @param: dat_combine - the log2 data matrix with samples in rows, features in columns
  # @param: main - the title for the heatmap
  # @param: log - whether your dat_combine is log2 transformed
  
  if (!log) {
    dat_combine <- log2(dat_combine+1)
  }
  
  dat <- data.frame(column = rep(1:ncol(dat_combine), rep(nrow(dat_combine), ncol(dat_combine))),
                    row = rep(1:nrow(dat_combine), ncol(dat_combine)),
                    value = (c(as.matrix(dat_combine))))
  
  p_heat <- ggplot(data = dat, aes(x = column, y = row, fill = value))+
    geom_tile()+
    theme_bw()+
    labs(title = main, x = "Genes", y = "Samples")
  
  return(p_heat)
}

DEA_eval <- function(dat_combine, groups = NULL, log = TRUE, failure= c("replace", "remove"), plot = FALSE){
  
  # This function compares real data and generated data from the perspective of differential expression analysis
  # @param: dat_combine - data matrix with the upper half real, lower half generated
  # @param: groups - vector containing the group information of samples
  # @param: log - whether the data is log2 transformed
  # @param: failure -  how to deal with failure genes "replace" replace the zero features in generated data with the real ones, "remove"  remove the failure genes in evaluation
  # @param: plot - whether to plot 
  # @return: return the concordance correlation of -log10pvalue and log2FC
  
  
  genes <- colnames(dat_combine)
  
  if(log){
    dat_combine <- 2^(dat_combine) - 1
    dat_combine[dat_combine < 0] <- 0
    dat_combine <- round(dat_combine)
  }
  dat_real <- dat_combine[1:(nrow(dat_combine)/2), ]
  dat_generated <- dat_combine[(nrow(dat_combine)/2+1):nrow(dat_combine), ]
  
  if(failure == "remove"){
    dat_real <- dat_real[, apply(dat_generated, 2, sd) != 0]
    dat_generated <- dat_generated[, apply(dat_generated, 2, sd) != 0]
    dat_combine <- rbind(dat_real, dat_generated)
    genes <- colnames(dat_combine)
  }
  else{
    dat_generated[ ,apply(dat_generated, 2, sd) == 0] <- dat_real[ ,apply(dat_generated, 2, sd) == 0]
    dat_combine <- rbind(dat_real, dat_generated)
  }
  if(nrow(dat_combine) <= 2){res <- NA}
  else{
    if(is.null(groups)){
      groups_real <- rep("A", nrow(dat_real))
      groups_real[sample(1:nrow(dat_real), round(0.5*nrow(dat_real)),replace = F)] <- "B"
      groups_generated <- groups_real
      groups_combine <- c(groups_real, groups_generated)
    }
    else if(length(groups) == 0.5*nrow(dat_combine)){
      groups_combine <- c(groups, groups)
      groups_real <-  groups
      groups_generated <-  groups
    }
    else{
      groups_combine <- groups
      groups_real <- groups_combine[1:(length(groups_combine)/2)]
      groups_generated <- groups_combine[(length(groups_combine)/2+1):length(groups_combine)]
    }
    
    dat_real <- as.matrix(dat_real)
    dat_generated <- as.matrix(dat_generated)
    
    DEA_real <- precision.seq::DE.voom(t(dat_real), groups_real)
    DEA_generated <- precision.seq::DE.voom(t(dat_generated), groups_generated)
    if(is.null(DEA_real$id.list)){DEA_real$id.list = colnames(dat_combine)}
    
    DEA_real$p.val[!is.finite(log10(DEA_real$p.val))] <- min(DEA_real$p.val[is.finite(log10(DEA_real$p.val))])
    
    DEA_generated$p.val[!is.finite(log10(DEA_generated$p.val))] <- min(DEA_generated$p.val[is.finite(log10(DEA_generated$p.val))])
    
    #### generated vs real
    res <- c(ccc_log10pvalue = round(epiR::epi.ccc((-log10(DEA_real$p.val[genes])), (-log10(DEA_generated$p.val[genes])))$rho.c[1], 4),
             ccc_log2FC = round(epiR::epi.ccc(DEA_real$log2.FC[DEA_real$id.list], DEA_generated$log2.FC[DEA_real$id.list])$rho.c[1], 4))
    
    
    if(plot){
      
      layout(matrix(c(1, 2, 3), nrow = 1, byrow = T))
      
      # mean vs sd
      plot(apply(dat_real, 2, mean), apply(dat_real, 2, sd),
           main = "mean vs sd", xlab = "Feature mean", ylab = "Feature std")
      points(apply(dat_generated, 2, mean), apply(dat_generated, 2, sd), col = "red")
      legend("topright", pch = c(1,1), col = c("black","red"), legend = c("Real", "Generated"))
      
      # -log10 pvalue ccc
      plot(-log10((DEA_real$p.val[genes])), -log10((DEA_generated$p.val[genes])),
           xlab = "Real",ylab = "Generated",
           main = paste("-log10 pvalues, ccc=", res$ccc_log10pvalue, sep = ""))
      abline(a = 0, b = 1)
      
      # log2FC of DE or all genes
      plot(DEA_real$log2.FC[DEA_real$id.list], DEA_generated$log2.FC[DEA_real$id.list],
           xlab = "Real",ylab = "Generated", main = paste( "log2FC of real DE genes, ccc =",res$ccc_log2FC, sep=""))
      abline(a = 0, b = 1)
      
      layout(1)
    }
  }
  
  return(res)
}

cluster_eval <- function(dat_combine, groups = NULL, log = TRUE, failure= c("replace", "remove"), plot = FALSE){
  # This function perform the clustering on the combined datasets
  # @param: dat_combine - data matrix with the upper half real, lower half generated
  # @param: groups - vector containing the group information of samples
  # @param: log - whether the data is log2 transformed
  # @param: failure -  how to deal with failure genes "replace" replace the zero features in generated data with the real ones, "remove"  remove the failure genes in evaluation
  # @param: plot - whether to plot 
  # @return: if groups are given, return ARI index, otherwise, return cARI = 1-ARI(datatype, clusters)
  
  if(!log){
    dat_combine <- log2(dat_combine + 1)
  }
  
  dat_real <- dat_combine[1:(nrow(dat_combine)/2), ]
  dat_generated <- dat_combine[(nrow(dat_combine)/2+1):nrow(dat_combine), ]
  
  if(failure == "remove"){
    dat_real <- dat_real[, apply(dat_generated, 2, sd) != 0]
    dat_generated <- dat_generated[, apply(dat_generated, 2, sd) != 0]
    dat_combine <- rbind(dat_real, dat_generated)
  }
  else if(failure == "replace"){
    dat_generated[ ,apply(dat_generated, 2, sd) == 0] <- dat_real[ ,apply(dat_generated, 2, sd) == 0]
    dat_combine <- rbind(dat_real, dat_generated)
  }
  
  if(is.null(groups)){no_clusters <- 2}
  else{
    if(length(groups) == 0.5*nrow(dat_combine)){
      groups_combine <- c(groups, groups)
      groups_real <-  groups
      groups_generated <-  groups
    }
    else{
      groups_combine <- groups
      groups_real <- groups_combine[1:(length(groups_combine)/2)]
      groups_generated <- groups_combine[(length(groups_combine)/2+1):length(groups_combine)]
    }
    no_clusters <- length(unique(groups_combine))
  }
  
  datatype <- rep(c("Real", "Generated"), rep(nrow(dat_combine)/2, 2))
  d <- dist(dat_combine, method = "euclidean") 
  fit <- hclust(d, method = "ward.D2")
  clusters <- cutree(fit, k = no_clusters)
  if(is.null(groups)){
    if(plot){
      rect.hclust(fit, k = no_clusters, border = "red")
      df <- data.frame(dataset = datatype, cluster = as.factor(clusters), sample = rep(1:(nrow(dat_combine)/2), 2)) 
      ggplot(data = df)+
        geom_point(aes(x = sample, y = dataset, color = cluster))+
        labs(x="Sample index", y="Dataset")
    }
    return(1-abs(round(aricode::ARI(datatype, clusters),4)))
  }
  else{
    if(plot){
      rect.hclust(fit, k = no_clusters, border = "red")
      df <- data.frame(groups = groups_combine, dataset = datatype, 
                       cluster = as.factor(clusters), sample = rep(1:(nrow(dat_combine)/2), 2)) 
      ggplot(data = df)+
        geom_point(aes(x = sample, y = groups, color = cluster))+
        facet_wrap(vars(dataset))+
        labs(x="Sample index", y="Groups")
    }
    return(round(aricode::ARI(groups_combine, clusters ),4))
  }
  
  
}

tSNE_eval <- function(dat_combine, groups = NULL, log = TRUE, failure= c("replace", "remove")){
  # This function perform tSNE dimension reduction to real and generated combined dataset 
  # and check whether the main variation lies between real and generated samples
  # @param: dat_combine - data matrix with samples in row, features in col, first N samples are real, second N samples are generated
  # @param: groups - groups information, samples from different groups are represented by different point shape
  # @param: log - whether the data are log2 transformed
  # @param: failure -  how to deal with failure genes "replace" replace the zero features in generated data with the real ones, "remove"  remove the failure genes in evaluation
  if(!log){
    dat_combine <- log2(dat_combine+1)
  }
  
  dat_real <- dat_combine[1:(nrow(dat_combine)/2), ]
  dat_generated <- dat_combine[(nrow(dat_combine)/2+1):nrow(dat_combine), ]
  
  if(failure == "remove"){
    dat_real <- dat_real[, apply(dat_generated, 2, sd) != 0]
    dat_generated <- dat_generated[, apply(dat_generated, 2, sd) != 0]
    dat_combine <- rbind(dat_real, dat_generated)
  }
  else if(failure == "replace"){
    dat_generated[ ,apply(dat_generated, 2, sd) == 0] <- dat_real[ ,apply(dat_generated, 2, sd) == 0]
    dat_combine <- rbind(dat_real, dat_generated)
  }
  if(length(groups) == 0.5*nrow(dat_combine)){
    groups_combine <- c(groups, groups)
    groups_real <-  groups
    groups_generated <-  groups
  }
  else{
    groups_combine <- groups
    groups_real <- groups_combine[1:(length(groups_combine)/2)]
    groups_generated <- groups_combine[(length(groups_combine)/2+1):length(groups_combine)]
  }
  mat <- as.matrix(dat_combine)
  datatype <- rep(c("real","generated"), rep(nrow(mat)/2, 2))
  tSNE_fit <- Rtsne::Rtsne(mat)
  tSNE_df <- tSNE_fit$Y %>% 
    as.data.frame() %>%
    rename(tSNE1 = "V1", tSNE2 = "V2") %>%
    mutate(ID = row_number())
  if(is.null(groups)){
    factor_df <- data.frame(datatype = datatype) %>% mutate(ID=row_number())
    plot_df <- tSNE_df %>% inner_join(factor_df, by="ID")
    plot_df %>% ggplot(aes(x = tSNE1, y = tSNE2, color = datatype))+
      geom_point()+
      theme(legend.position="bottom")+
      labs(title="tSNE plot of combined samples")
  }
  else{
    factor_df <- data.frame(groups = groups_combine, datatype = datatype) %>% mutate(ID=row_number())
    plot_df <- tSNE_df %>% inner_join(factor_df, by="ID")
    plot_df %>% ggplot(aes(x = tSNE1, y = tSNE2, color = datatype, shape = groups))+
      geom_point()+
      theme(legend.position="bottom")+
      labs(title="tSNE plot of combined samples")+
      themebw()
  }
}

UMAP_eval <- function(dat_combine, groups = NULL, log = TRUE, failure= c("replace", "remove")){
  # This function perform UMAP dimension reduction to real and generated combined dataset 
  # and check whether the main variation lies between real and generated samples
  # @param: dat_combine - data matrix with samples in row, features in col, first N samples are real, second N samples are generated
  # @param: groups - groups information, samples from different groups are represented by different point shape
  # @param: log - whether the data are log2 transformed
  # @param: failure -  how to deal with failure genes "replace" replace the zero features in generated data with the real ones, "remove"  remove the failure genes in evaluation
  if(!log){
    dat_combine <- log2(dat_combine+1)
  }
  
  dat_real <- dat_combine[1:(nrow(dat_combine)/2), ]
  dat_generated <- dat_combine[(nrow(dat_combine)/2+1):nrow(dat_combine), ]
  
  if(failure == "remove"){
    dat_real <- dat_real[, apply(dat_generated, 2, sd) != 0]
    dat_generated <- dat_generated[, apply(dat_generated, 2, sd) != 0]
    dat_combine <- rbind(dat_real, dat_generated)
  }
  else if(failure == "replace"){
    dat_generated[ ,apply(dat_generated, 2, sd) == 0] <- dat_real[ ,apply(dat_generated, 2, sd) == 0]
    dat_combine <- rbind(dat_real, dat_generated)
  }
  if(length(groups) == 0.5*nrow(dat_combine)){
    groups_combine <- c(groups, groups)
    groups_real <-  groups
    groups_generated <-  groups
  }
  else{
    groups_combine <- groups
    groups_real <- groups_combine[1:(length(groups_combine)/2)]
    groups_generated <- groups_combine[(length(groups_combine)/2+1):length(groups_combine)]
  }
  mat <- as.matrix(dat_combine)
  datatype <- rep(c("real","generated"), rep(nrow(mat)/2, 2))
  UMAP_fit <- umap::umap(mat)
  UMAP_df <- UMAP_fit$layout %>% 
    as.data.frame() %>%
    rename(UMAP1 = "V1", UMAP2 = "V2") %>%
    mutate(ID = row_number())
  if(is.null(groups)){
    factor_df <- data.frame(datatype = datatype) %>% mutate(ID=row_number())
    plot_df <- UMAP_df %>% inner_join(factor_df, by="ID")
    p_umap <- plot_df %>% ggplot(aes(x = UMAP1, y = UMAP2, 
                                     color = datatype))+
      geom_point()+
      theme(legend.position="bottom")+
      labs(title="UMAP plot of combined samples")
  }
  else{
    factor_df <- data.frame(groups = groups_combine, datatype = datatype) %>% mutate(ID=row_number())
    plot_df <- UMAP_df %>% inner_join(factor_df, by="ID")
    p_umap <- plot_df %>% ggplot(aes(x = UMAP1, y = UMAP2, 
                                     color = datatype, shape = groups))+
      geom_point()+
      theme(legend.position="bottom")+
      labs(title="UMAP plot of combined samples")+
      theme_bw()
  }
  return(list(plot_df = plot_df, p_umap = p_umap))
}

ccpos_eval <- function(dat_combine, log = TRUE, failure = c("replace", "remove"), coords, poly = F, thres = NULL){
  # This function perform DANA comparison of concordance correlation which measures the preservation of biological signals real versus generated
  # @param: dat_combine - data matrix with samples in row, features in col, first N samples are real, second N samples are generated
  # @param: log - whether the data are log2 transformed
  # @param: failure -  how to deal with failure genes "replace" replace the zero features in generated data with the real ones, "remove"  remove the failure genes in evaluation
  # @param: poly - whether the features are all positive controls
  # @param: coords - coordinates data frame based on miRBase miRNA definitions
  # @param: thres - threshold for positive controls features
  if(is.null(coords)){cc <- NA}
  else{
    if(log){
      dat_combine <- 2^(dat_combine) - 1
      dat_combine[dat_combine < 0] <- 0
      dat_combine <- round(dat_combine)
    }
    
    dat_real <- dat_combine[1:(nrow(dat_combine)/2), ]
    dat_generated <- dat_combine[(nrow(dat_combine)/2+1):nrow(dat_combine), ]
    
    if(failure == "remove"){
      dat_real <- dat_real[, apply(dat_generated, 2, sd) != 0]
      dat_generated <- dat_generated[, apply(dat_generated, 2, sd) != 0]
      dat_combine <- rbind(dat_real, dat_generated)
    }
    else if(failure == "replace"){
      dat_generated[ ,apply(dat_generated, 2, sd) == 0] <- dat_real[ ,apply(dat_generated, 2, sd) == 0]
      dat_combine <- rbind(dat_real, dat_generated)
    }
    
    if(nrow(dat_combine) <= 2){cc <- NA}
    else{
      colnames(dat_generated) <- colnames(dat_real)
      rownames(dat_generated) <- rownames(dat_real)
      
      clusters <- DANA::defineClusters(genes = coords$name, chr = coords$chr, 
                                       pos = (coords$start + coords$end)/2)
      
      if(poly){posControls <- colnames(dat_combine)}
      else{
        posControls <- DANA::defineControls(raw = t(dat_real), 
                                            tZero = 2, tPoor = 10, tWell = thres,
                                            clusters = clusters)$posControls
      }
      if(!is.null(posControls)){
        invisible(capture.output(
          corPos_real <- DANA::partialCor(t(dat_real[, posControls]), scale = TRUE)))
        invisible(capture.output(
          corPos_generated <- DANA::partialCor(t(dat_generated[, posControls]), scale = TRUE)))
        cc <- compute.cc(corPos_real, corPos_generated, listClusters(as.factor(clusters)))
      }
      else{cc <- NA}
    }
  }
  return(cc)
}

summary_eval <- function(dat_combine, log = TRUE, failure= c("replace", "remove")){
  # This function compares features in generated and real dataset using mean absolute deviation on log2
  # @param: dat_combine - data matrix with samples in row, features in col, first N samples are real, second N samples are generated
  # @param: log - whether the data are log2 transformed
  # @param: failure -  how to deal with failure genes "replace" replace the zero features in generated data with the real ones, "remove"  remove the failure genes in evaluation
  # @return: a list including median absolute deviation of real-generated of feature mean, var, cv, zero fraction
  if(!log){
    dat_combine <- log2(dat_combine + 1)
  }
  
  dat_real <- dat_combine[1:(nrow(dat_combine)/2), ]
  dat_generated <- dat_combine[(nrow(dat_combine)/2+1):nrow(dat_combine), ]
  
  if(failure == "remove"){
    dat_real <- dat_real[, apply(dat_generated, 2, sd) != 0]
    dat_generated <- dat_generated[, apply(dat_generated, 2, sd) != 0]
    dat_combine <- rbind(dat_real, dat_generated)
  }
  else if(failure == "replace"){
    dat_generated[ ,apply(dat_generated, 2, sd) == 0] <- dat_real[ ,apply(dat_generated, 2, sd) == 0]
    dat_combine <- rbind(dat_real, dat_generated)
  }
  if(nrow(dat_combine) <= 2){return(NA)}
  else{
    mu_generated <- apply(dat_generated, 2, mean)
    sigma_generated <- apply(dat_generated, 2, sd)
    zero_generated <- apply(dat_generated, 2, function(x){mean(x == 0)})
    cv_generated <- sigma_generated/mu_generated
    
    zero_real <- apply(dat_real, 2, function(x){mean(x == 0)})
    mu_real <- apply(dat_real, 2, mean)
    sigma_real <- apply(dat_real, 2, sd)
    cv_real <- sigma_real/mu_real
    
    ks_mean <- mad(mu_generated-mu_real)
    ks_sigma <- mad(sigma_generated-sigma_real)
    ks_zero <-mad(zero_generated-zero_real)
    ks_cv <- mad(cv_generated-cv_real)
    
    return(list(ks_mean = ks_mean, ks_sigma = ks_sigma, ks_zero = ks_zero, ks_cv = ks_cv))
  }
}

fail_features_eval <- function(dat_combine){
  # This function evaluate the generated data with respect to the failue features
  # Define the failure features to be the features with standard deviation 0
  # @param: dat_combine - data matrix with samples in row, features in col, first N samples are real, second N samples are generated
  # @return: proportion of failure features among all the features, specifially, failure genes unique in generated dataset.
  
  dat_real <- dat_combine[1:(nrow(dat_combine)/2), ]
  dat_generated <- dat_combine[(nrow(dat_combine)/2+1):nrow(dat_combine), ]
  
  
  fail_real <- apply(dat_real, 2, sd) == 0
  fail_generated <- apply(dat_generated, 2, sd) == 0
  fail_prop <- mean(apply(data.frame(fail_real = fail_real, fail_generated = fail_generated), 1, function(x){(!x['fail_real']) & x['fail_generated']}))
  return(fail_prop)
}

get_eval <- function(dataname, model_list, pilot_list, log, failure, poly, plot_first = FALSE){
  # This function conduct the DEA, clustering, DANA, summary stat evaluation
  # on the generated samples  under different models and different pilot size
  # @param: dataname - original dataset name
  # @param: model_list - a vector including all the candidates models
  # @param: pilot_list - a vector including all the candidates of pilot size
  # @param: groups - can be groups of N real samples, or groups of 2N real+generated samples
  # @param: log - whether the generated data are log2 transformed, usually TRUE
  # @param: failure -  how to deal with failure genes "replace" replace the zero features in generated data with the real ones, "remove"  remove the failure genes in evaluation
  # @param: poly - whether the features are all positive controls, if TRUE, DANA will use all genes as positive controls
  # @param: plot_first - whether or not plot the heatmap and UMAP of the first draw
  # @return: a data frame containing evaluation results of all model and pilot size combinations
  
  # set up
  real <- paste("../RealData/", dataname, ".csv", sep = "")
  pil <- length(pilot_list)
  repli <- 5
  draw <- 5
  mod <- length(model_list)
  df_res <- data.frame(model = rep(model_list, rep(pil * repli * draw, mod)),
                       pilot_size = rep(pilot_list, rep(repli * draw, pil)),
                       replication = rep(1:(repli * draw), pil * mod),
                       ccc_log10pvalue = rep(NA, pil * repli * draw * mod),
                       ccc_log2FC = rep(NA, pil * repli * draw * mod),
                       ARI = rep(NA, pil * repli * draw * mod),
                       fail_prop = rep(NA, pil * repli * draw * mod),
                       ccc_pos = rep(NA, pil * repli * draw * mod),
                       ks_mean = rep(NA, pil * repli * draw * mod),
                       ks_sd = rep(NA, pil * repli * draw * mod),
                       ks_zero = rep(NA, pil * repli * draw * mod))
  df_res$generated_list <- paste(dataname, df_res$model, df_res$pilot_size, sep = "_")
  if(substr(dataname, 1,3) == "RNA"){coords <- NULL}
  else{coords <- read.csv("../RealData/miRBase_coords.csv")}
  # evaluation
  for ( generated in unique(df_res$generated_list)) {
    cat(generated)
    generated_info <- get_bind_list(real, generated)
    dat_combine_list <- generated_info$bind_list
    groups_all <- generated_info$groups
    samples <- generated_info$samples
    df_index_store <- data.frame(ccc_log10pvalue = rep(NA, repli * draw),
                                 ccc_log2FC = rep(NA, repli * draw),
                                 ARI = rep(NA, repli * draw),
                                 fail_prop = rep(NA, repli * draw),
                                 ccc_pos = rep(NA, repli * draw),
                                 ks_mean = rep(NA, repli * draw),
                                 ks_sd = rep(NA, repli * draw),
                                 ks_zero = rep(NA, repli * draw),
                                 ks_cv =rep(NA, repli * draw))
    for (repe in 1:(repli*draw)) {
      dat_combine <- dat_combine_list[[repe]]
      N <- nrow(dat_combine)/2
      if(is.null(groups_all)){groups <- NULL}
      else{groups <- groups_all[c(1:N, (repe * N + 1) : ((repe+1) * N))]}
      if(repe == 1 & plot_first == TRUE){
        p_heat <- heatmap_eval(dat_combine = dat_combine, 
                               main = "Real in lower panel ", log = log)
        tryCatch({ 
          p_umap <- suppressWarnings(
            UMAP_eval(dat_combine = dat_combine, groups = groups_all,
                      log = log, failure = failure))
        }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        tryCatch({print(p_umap$p_umap)}, 
                 error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
        print(p_heat)
      }      
      tryCatch({ 
        df_index_store[repe, c("ccc_log10pvalue", "ccc_log2FC")] <- DEA_eval(dat_combine = dat_combine, groups = groups, log = log, failure = failure, plot = FALSE)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      tryCatch({ 
        df_index_store[repe, "ARI"] <- cluster_eval(dat_combine = dat_combine, groups = groups, log = log, failure = failure, plot = FALSE)
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      df_index_store[repe, "fail_prop"] <- fail_features_eval(dat_combine = dat_combine)
      tryCatch({ 
        df_index_store[repe, "ccc_pos"] <- ccpos_eval(dat_combine = dat_combine, failure = failure,
                                                      log = log, coords = coords, poly = poly,
                                                      thres = 64) 
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
      tryCatch({ 
        df_index_store[repe, c("ks_mean", "ks_sd", "ks_zero", "ks_cv")] <-
          unlist(suppressWarnings(summary_eval(dat_combine = dat_combine,  failure = failure, log = log)))
      }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
    }
    df_res[df_res$generated_list == generated, 
           c("ccc_log10pvalue", "ccc_log2FC", "ARI", "fail_prop",
             "ccc_pos", "ks_mean", "ks_sd", "ks_zero", "ks_cv")] <- df_index_store 
  }
  return(df_res)
}

visual_eval <- function(df_res,  model_list = NULL){
  # This function help visualization of the evaluation results
  # @params: df_res - data frame containing at least pilot_size, model, 
  #                   ccc_log10pvalue, ccc_log2FC, ARI, fail_prop
  
  df_res$pilot_size <- factor(df_res$pilot_size, levels = pilot_list)
  df_res$model <- factor(df_res$model, levels = model_list)
  p_ccc_log10pvalue <- ggplot(data = df_res)+
    geom_boxplot(aes(x = pilot_size, y = ccc_log10pvalue, fill = model), 
                 show.legend = F, width = 0.6)+
    scale_fill_brewer(palette = "PiYG")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(x="Pilot Size", y = "CCC of -log10(p-values)")
  p_ccc_log2FC <- ggplot(data = df_res)+
    geom_boxplot(aes(x = pilot_size, y = ccc_log2FC, fill = model), 
                 show.legend = F, width = 0.6)+
    scale_fill_brewer(palette = "PiYG")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(x="Pilot Size", y = "CCC of log2FC")
  p_ARI <- ggplot(data = df_res)+
    geom_boxplot(aes(x = pilot_size, y = ARI, fill = model), 
                 show.legend = F, width = 0.6)+
    theme_bw()+ 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_fill_brewer(palette = "PiYG")+
    labs(x="Pilot Size", y = "ARI/cARI of sample groups")
  
  p_succ_prop <- ggplot(data = df_res %>% 
                          group_by(model, pilot_size) %>% 
                          summarise(success =  1-mean(fail_prop)))+
    geom_point(aes(x = pilot_size, y = success, color = model, group = model), 
               show.legend = T)+
    geom_line(aes(x = pilot_size, y = success, color = model, group = model), 
              show.legend = F)+
    theme_bw()+ 
    ylim(0,1)+
    theme(legend.position = c(0.5, 0.2))+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    scale_color_brewer(palette = "PiYG")+
    labs(x="Pilot Size", y = "1-%(all zero genes)")
  p_ccc_pos <- ggplot(data = df_res)+
    geom_boxplot(aes(x = pilot_size, y = ccc_pos, fill = model),
                 show.legend = F, width = 0.6)+
    scale_fill_brewer(palette = "PiYG")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(x="Pilot Size", y = "CCC of partial corr")+
    guides(fill = guide_legend(nrow = 2, byrow = T))
  p_ks_mean <- ggplot(data = df_res)+
    geom_boxplot(aes(x = pilot_size, y = ks_mean, fill = model),
                 show.legend = F, width = 0.6)+
    scale_fill_brewer(palette = "PiYG")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(x="Pilot Size", y = "MAD of feature mean")
  p_ks_sd <- ggplot(data = df_res)+
    geom_boxplot(aes(x = pilot_size, y = ks_sd, fill = model),
                 show.legend = F, width = 0.6)+
    scale_fill_brewer(palette = "PiYG")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(x="Pilot Size", y = "MAD of feature sd")
  p_ks_zero <- ggplot(data = df_res)+
    geom_boxplot(aes(x = pilot_size, y = ks_zero, fill = model),
                 show.legend = F, width = 0.6)+
    scale_fill_brewer(palette = "PiYG")+
    theme_bw()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())+
    labs(x="Pilot Size", y = "MAD feature zero fraction")
  plot_grid(p_ccc_log10pvalue, p_ccc_log2FC, p_ARI, p_succ_prop, 
            p_ccc_pos, p_ks_mean, p_ks_sd, p_ks_zero, nrow = 2, ncol=4 )
}


###########################################
## Functions from Visualization.Rmd ####
###########################################
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggh4x)
library(RColorBrewer)
library(ggsci)
library(beanplot)
library(ggforce)
library(cowplot)

gg_default <- function(p){
  pp <- p +
    theme(legend.position = "bottom")+
    # guides(fill=guide_legend(nrow = 1,byrow = TRUE))+
    theme(strip.text.x = element_text(face = "bold", size = 12),
          strip.text.y = element_text(face = "bold", size = 10),
          legend.title=element_text(size = 10), 
          legend.text=element_text(size = 9),
          axis.text = element_text(size = 10),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank())
  return(pp)
}
UMAP_eval <- function(dat_combine, groups = NULL, log = TRUE, failure= c("replace", "remove")){
  if(!log){
    dat_combine <- log2(dat_combine+1)
  }
  
  dat_real <- dat_combine[1:(nrow(dat_combine)/2), ]
  dat_generated <- dat_combine[(nrow(dat_combine)/2+1):nrow(dat_combine), ]
  
  if(failure == "remove"){
    dat_real <- dat_real[, apply(dat_generated, 2, sd) != 0]
    dat_generated <- dat_generated[, apply(dat_generated, 2, sd) != 0]
    dat_combine <- rbind(dat_real, dat_generated)
  }
  else if(failure == "replace"){
    dat_generated[ ,apply(dat_generated, 2, sd) == 0] <- dat_real[ ,apply(dat_generated, 2, sd) == 0]
    dat_combine <- rbind(dat_real, dat_generated)
  }
  if(length(groups) == 0.5*nrow(dat_combine)){
    groups_combine <- c(groups, groups)
    groups_real <-  groups
    groups_generated <-  groups
  }
  else{
    groups_combine <- groups
    groups_real <- groups_combine[1:(length(groups_combine)/2)]
    groups_generated <- groups_combine[(length(groups_combine)/2+1):length(groups_combine)]
  }
  mat <- as.matrix(dat_combine)
  datatype <- rep(c("real","generated"), rep(nrow(mat)/2, 2))
  UMAP_fit <- umap::umap(mat)
  UMAP_df <- UMAP_fit$layout %>% 
    as.data.frame() %>%
    rename(UMAP1 = "V1", UMAP2 = "V2") %>%
    mutate(ID = row_number())
  if(is.null(groups)){
    factor_df <- data.frame(datatype = datatype) %>% mutate(ID=row_number())
    plot_df <- UMAP_df %>% inner_join(factor_df, by="ID")
    p_umap <- plot_df %>% ggplot(aes(x = UMAP1, y = UMAP2, 
                                     color = datatype))+
      geom_point()+
      theme(legend.position="bottom")+
      labs(title="UMAP plot of combined samples")
  }
  else{
    factor_df <- data.frame(groups = groups_combine, datatype = datatype) %>% mutate(ID=row_number())
    plot_df <- UMAP_df %>% inner_join(factor_df, by="ID")
    p_umap <- plot_df %>% ggplot(aes(x = UMAP1, y = UMAP2, 
                                     color = datatype, shape = groups))+
      geom_point()+
      theme(legend.position="bottom")+
      labs(title="UMAP plot of combined samples")+
      theme_bw()
  }
  return(list(plot_df = plot_df, p_umap = p_umap))
}
vis_all <- function(dataname, nickname, cancer_type = c(1, 2)){
  # This function visualization the results of DGM evaluation
  # dataname = "SKCMPositive_4", nickname = "skcm_pos", cancer_type = 1
  
  
  # Step1: load data ####
  default <- "FinalResAug2023/"
  load(paste(default, dataname, ".RData", sep = ""))
  load(paste(default, dataname, "_AEhead.RData", sep = ""))
  load(paste(default, dataname, "_Normalization.RData", sep = ""))
  load(paste(default, dataname, "_Batch.RData", sep = ""))
  load(paste(default, dataname, "_Epoch.RData", sep = ""))
  overall <- get(nickname)
  offline <- get(paste(nickname, "_AEhead", sep = ""))
  norm <- get(paste(nickname, "_norm", sep = ""))
  batch <-  get(paste(nickname, "_batch", sep = ""))
  epoch <-  get(paste(nickname, "_epoch", sep = ""))
  
  # Step2: generate dataframes for plotting  ####
  if(cancer_type == 1){
    # Evaluation criteria: 1-Prop.successful.genes, 2-CCC.pos.control, 3-cARI, 4-MAD
    mycolors <- c(brewer.pal(n = 8, name = "Blues")[c(2,4,6)], brewer.pal(n = 8, name = "Oranges")[c(2,4,6)], brewer.pal(n = 8, name = "Greens")[c(2,4,6)])
    plot_width <- c(10,10,16,7,7,7)
    plot_height <- c(6,6,6,6,6,6)
    df_proc <- function(res){
      res$succ_prop <- 1 - res$fail_prop
      res$pilot_size <- factor(res$pilot_size, levels = as.character(sort(as.numeric(as.character(unique(res$pilot_size))))))
      if("modelname" %in% colnames(res)){model <- res[,"modelname"]}
      else{model <- res[,"model"]}
      if(length(unique(model)) <= 3){lvs <- c("VAE1-10","MAF","WGANGP")}
      else{lvs <- c("VAE1-1", "VAE1-5", "VAE1-10","RealNVP","GLOW","MAF","GAN", "WGAN", "WGANGP")}
      res$model <- factor(model, levels = lvs)
      res$modeltype <- rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(substr(res$model[i], 1,3)=="VAE"){res$modeltype[i] <- "VAEs"}
        else if(substr(res$model[i], 1,3) %in% c("GAN", "WGA")){res$modeltype[i] <- "GANs"}
        else{res$modeltype[i] <- "FLOWs"}
      }
      res$modeltype <- factor(res$modeltype, levels = c("VAEs","FLOWs","GANs"))
      res$ccc_pos[res$ccc_pos<0 | is.na(res$ccc_pos)] <- 0
      res$ARI[res$ARI<0 | is.na(res$ARI)] <- 0
      res <- tidyr::pivot_longer(res, cols = c("ARI", "ccc_pos", "succ_prop"), names_to = "Metric", values_to = "Rate")
      res$Metric <- factor(res$Metric, levels = c("succ_prop","ARI","ccc_pos"),
                           labels = c("1 - Pct(0-markers)", "cARI", "CCCPCC"))
      if("normalization" %in% colnames(res)){res$normalization <- factor(res$normalization, levels=c('None','TC','TMM','UQ'))}
      return(res)
    }
    heatmap_proc <- function(res){
      if("modelname" %in% colnames(res)){model <- res[,"modelname"]}
      else{model <- res[,"model"]}
      if(length(unique(model)) <= 4){lvs <- c("VAE1-10","MAF","WGANGP")}
      else{lvs <- c("VAE1-1", "VAE1-5", "VAE1-10","RealNVP","GLOW","MAF","GAN", "WGAN", "WGANGP")}
      res$model <- factor(model, levels = rev(lvs))
      res$pilot_size <- factor(res$pilot_size, levels = as.character(sort(as.numeric(as.character(unique(res$pilot_size))))))
      res$modeltype <- rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(substr(res$model[i], 1,3)=="VAE"){res$modeltype[i] <- "VAEs"}
        else if(substr(res$model[i], 1,3) %in% c("GAN", "WGA")){res$modeltype[i] <- "GANs"}
        else{res$modeltype[i] <- "FLOWs"}
      }
      res$modeltype <- factor(res$modeltype, levels = c("VAEs","FLOWs","GANs"))
      res <- res %>% group_by(pilot_size, model, modeltype) %>%
        summarise(
          Mean = mean(ks_mean, na.rm=TRUE),
          Sd = mean(ks_sd, na.rm=TRUE),
          Sparsity = mean(ks_zero, na.rm=TRUE)
        )
      
      res <- tidyr::pivot_longer(res, cols = c("Mean", "Sd", "Sparsity"), names_to = "Metric", values_to = "Rate")
      res$Metric <- factor(res$Metric, levels = (c("Mean", "Sd", "Sparsity")), labels = c("MAD(Mean)","MAD(SD)","MAD(Sparsity)"))
      res$text <- round(res$Rate, 3)
      res$text[res$Rate > 6 | is.na(res$Rate)] <- ">6"
      res$Rate[res$Rate > 6 | is.na(res$Rate)] <- 6
      res$modelorder <-  rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(res$model[i] %in% c("VAE1-1","GAN","RealNVP")){res$modelorder[i] <- "VAE1-1/RealNVP/GAN"}
        else if(res$model[i] %in% c("VAE1-5","WGAN","GLOW")){res$modelorder[i] <- "VAE1-5/GLOW/WGAN"}
        else{res$modelorder[i] <- "VAE1-10/MAF/WGANGP"}
      }
      res$modelorder <- factor(res$modelorder, levels = c("VAE1-10/MAF/WGANGP", "VAE1-5/GLOW/WGAN", "VAE1-1/RealNVP/GAN"))
      if("normalization" %in% colnames(res)){res$normalization <- factor(res$normalization, levels=c('None','TC','TMM','UQ'))}
      return(res)
    }
  }
  else{
    # Evaluation criteria: 1-Prop.successful.genes, 2-CCC.pos.control, 3-ARI, 4-MAD, 5-ccc.log10pvalue, 6-ccc.log2FC
    plot_width <- c(6,4,8,6,6,6)
    plot_height <- c(9.5, 9.5, 9.5, 9.5, 9.5, 9.5)
    mycolors <- c(brewer.pal(n = 8, name = "Blues")[c(2,4,6)],brewer.pal(n = 8, name = "Oranges")[6])
    df_proc <- function(res){
      res$succ_prop <- 1 - res$fail_prop
      res$pilot_size <- factor(res$pilot_size, levels = as.character(sort(as.numeric(as.character(unique(res$pilot_size))))))
      if("modelname" %in% colnames(res)){model <- res[,"modelname"]}
      else{model <- res[,"model"]}
      if(length(unique(model)) <= 2){lvs <- c("CVAE1-10","MAF")}
      else{lvs <- c("CVAE1-1", "CVAE1-5", "CVAE1-10","MAF")}
      res$model <- factor(model, levels = lvs)
      res$modeltype <- rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(substr(res$model[i], 1,3)=="CVA"){res$modeltype[i] <- "VAEs"}
        else{res$modeltype[i] <- "FLOWs"}
      }
      res$modeltype <- factor(res$modeltype, levels = c("VAEs","FLOWs"))
      res$ccc_pos[res$ccc_pos<0 | is.na(res$ccc_pos)] <- 0
      res$ARI[res$ARI<0 | is.na(res$ARI)] <- 0
      res$ccc_log10pvalue[res$ccc_log10pvalue < 0 | is.na(res$ccc_log10pvalue)] <- 0
      res$ccc_log2FC[res$ccc_log2FC < 0 | is.na(res$ccc_log2FC)] <- 0
      if("normalization" %in% colnames(res)){res$normalization <- factor(res$normalization, levels=c('None','TC','TMM','UQ','DESeq'))}
      res <- tidyr::pivot_longer(res, cols = c("ccc_log10pvalue", "ccc_log2FC","ARI", "succ_prop","ccc_pos"), names_to = "Metric", values_to = "Rate")
      res$Metric <- factor(res$Metric, levels = c("succ_prop","ARI","ccc_pos","ccc_log10pvalue", "ccc_log2FC"),
                           labels = c("1 - Pct(0-markers)", "ARI","CCCPCC", "CCC of -log10(pvalue)", "CCC of log2FC"))
      return(res)
    }
    
    heatmap_proc <- function(res){
      if("modelname" %in% colnames(res)){model <- res[,"modelname"]}
      else{model <- res[,"model"]}
      if(length(unique(model)) <= 2){lvs <- c("CVAE1-10","MAF")}
      else{lvs <- c("CVAE1-1", "CVAE1-5", "CVAE1-10","MAF")}
      res$model <- factor(model, levels = lvs)
      res$pilot_size <- factor(res$pilot_size, levels = as.character(sort(as.numeric(as.character(unique(res$pilot_size))))))
      res$modeltype <- rep(NA, nrow(res))
      for (i in 1:nrow(res)) {
        if(substr(res$model[i], 1,3)=="CVA"){res$modeltype[i] <- "VAEs"}
        else{res$modeltype[i] <- "FLOWs"}
      }
      res$modeltype <- factor(res$modeltype, levels = c("VAEs","FLOWs"))
      res <- res %>% group_by(pilot_size, model, modeltype) %>%
        summarise(
          Mean = mean(ks_mean,na.rm=TRUE),
          Sd = mean(ks_sd,na.rm=TRUE),
          Sparsity = mean(ks_zero,na.rm=TRUE)
        )
      
      res <- tidyr::pivot_longer(res, cols = c("Mean", "Sd", "Sparsity"), names_to = "Metric", values_to = "Rate")
      res$Metric <- factor(res$Metric, levels = (c("Mean", "Sd", "Sparsity")), labels = c("MAD(Mean)","MAD(SD)","MAD(Sparsity)"))
      res$text <- round(res$Rate, 3)
      res$text[res$Rate > 6 | is.na(res$Rate)] <- ">6"
      res$Rate[res$Rate > 6 | is.na(res$Rate)] <- 6
      res$modelorder <- factor(res$model, levels=rev(levels(res$model)))
      res$modeltype <- "All models"
      if("normalization" %in% colnames(res)){res$normalization <- factor(res$normalization, levels=c('None','TC','TMM','UQ','DESeq'))}
      return(res)
    }
  }
  # Step 3: make graphs  ####
  gg_default <- function(p){
    pp <- p +
      theme(legend.position = "bottom")+
      # guides(fill=guide_legend(nrow = 1,byrow = TRUE))+
      theme(strip.text.x = element_text(face = "bold", size = 12),
            strip.text.y = element_text(face = "bold", size = 10),
            legend.title=element_text(size = 10), 
            legend.text=element_text(size = 9),
            axis.text = element_text(size = 10),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_blank())
    return(pp)
  }
  #### p11 overall boxplot ####
  res <- df_proc(overall)
  p11 <- ggplot(data = res)+ 
    geom_violin(aes(x = pilot_size, y = Rate, fill = model), position = position_dodge(0.6), width=1)+
    # geom_boxplot(aes(x = pilot_size, y = Rate, fill = model), show.legend = T, width = 0.15,position = position_dodge(0.9))+
    facet_grid(row = vars(Metric), col = vars(modeltype))+
    theme_bw()+
    labs(x = "Pilot size", y = "", switch = "y", fill = "Models")+
    scale_fill_manual(values = mycolors)+
    theme(panel.border = element_rect(colour = "black", fill = NA))+
    ylim(0,1)+
    guides(fill=guide_legend(nrow = 1,byrow = TRUE))
  if(dataname %in% c("SKCM","BRCA","SKCMLAML")){
    p11 <- p11+geom_hline(data = data.frame(Metric = "1 - Pct(0-markers)",y=0.85), aes(yintercept = y), linetype = "dashed")
  }
  else{
    p11 <- p11+geom_hline(data = data.frame(Metric = "1 - Pct(0-markers)",y=1), aes(yintercept = y), linetype = "dashed")
  }
  pdf(file = paste(default, dataname, ".pdf", sep = ""), width = plot_width[1], height = plot_height[1])
  print(gg_default(p11))
  dev.off()
  #### p12 overall heatmap ####
  res <- heatmap_proc(overall)
  p12 <- ggplot(res, aes(x = pilot_size , y = modelorder, fill = Rate)) + 
    geom_tile(show.legend = T)+
    labs(y = "",x = "Pilot size", fill = "MAD")+
    facet_nested(Metric~modeltype)+
    geom_text(aes(pilot_size, modelorder, label = text), color = "black", size = 4) +
    theme_minimal()+
    # scale_y_discrete(position = "left", labels=c( "Mean" = "MAD(Mean)", "Sd" = "MAD(SD)", "Sparsity" = "MAD(Sparsity)"))+
    scale_fill_gradientn(colours = terrain.colors(10))+
    theme(axis.text.y = element_text(angle = 0, hjust = 1, colour = "black", size = 8))+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_blank())
  pdf(file = paste(default, dataname, "_MAD.pdf", sep = ""), width = plot_width[2], height = plot_height[2])
  print(gg_default(p12))
  dev.off()
  
  #### p21 offline violine plot####
  res <- df_proc(offline)
  p21 <- ggplot(data = res, aes(x = pilot_size, y = Rate))+
    geom_violin(aes(color = offline), show.legend = T, position = position_dodge(0.9)) +
    geom_boxplot(aes(color = offline), width = 0.15, position = position_dodge(0.9)) +
    facet_grid(cols=vars(model), rows= vars(Metric), scale = "free_y")+
    theme_bw()+
    labs(x="Pilot size", y="", color = "Offline")+
    scale_y_continuous(position = "left")+
    theme(strip.background = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA),
          axis.line = element_blank())+
    scale_color_npg()+
    ylim(0,1)
  if(dataname %in% c("SKCM","BRCA","SKCMLAML")){
    p21 <- p21+geom_hline(data = data.frame(Metric = "1 - Pct(0-markers)",y=0.85), aes(yintercept = y), linetype = "dashed")
  }
  else{
    p21 <- p21+geom_hline(data = data.frame(Metric = "1 - Pct(0-markers)",y=1), aes(yintercept = y), linetype = "dashed")
  }
  pdf(file = paste(default, dataname, "_Offline.pdf", sep = ""), width = plot_width[3], height = plot_height[3])
  print(gg_default(p21))
  dev.off()
  
  
  #### p3 batch sina plot ####
  res <- df_proc(batch)
  res$batch_size <- paste(res$batch_frac, "of pilot size")
  p3 <- ggplot(res[res$batch_frac!= "50%",])+ 
    # geom_sina(aes(x = pilot_size, y = Rate, color = batch_frac), show.legend = T)+
    geom_violin(aes(x = pilot_size, y = Rate, color = batch_size), show.legend = T, position = position_dodge(0.9)) +
    geom_boxplot(aes(x = pilot_size, y = Rate, color = batch_size), width = 0.15, position = position_dodge(0.9)) +
    facet_nested(Metric~model)+
    theme_bw()+
    labs(x="Pilot size", y="", color = "Batch size")+
    scale_y_continuous(position = "left")+
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.line = element_blank())+
    scale_color_npg()+
    ylim(0,1)
  if(dataname %in% c("SKCM","BRCA","SKCMLAML")){
    p3 <- p3+geom_hline(data = data.frame(Metric = "1 - Pct(0-markers)",y=0.85), aes(yintercept = y), linetype = "dashed")
  }
  else{
    p3 <- p3+geom_hline(data = data.frame(Metric = "1 - Pct(0-markers)",y=1), aes(yintercept = y), linetype = "dashed")
  }
  pdf(file = paste(default, dataname, "_Batch.pdf", sep = ""), width = plot_width[4], height = plot_height[4])
  print(gg_default(p3))
  dev.off()
  
  #### p4 epoch line plot ####
  # if(cancer_type == 1){lvs <- c("VAE1-10","MAF","WGANGP")}
  # else{lvs <- c("CVAE1-10","MAF")}
  # epoch$model <- factor(epoch$modelname, levels = lvs)
  # res <- epoch %>% group_by(model, pilot_size, epoch) %>%
  #   summarise(
  #     sd_ccc_pos = sd(ccc_pos,na.rm = T),
  #     mean_ccc_pos = mean(ccc_pos,na.rm = T),
  #     sd_ARI = sd(ARI,na.rm = T),
  #     mean_ARI = mean(ARI,na.rm = T),
  #   )
  # 
  # res_sd <- tidyr::pivot_longer(res, cols = c("sd_ccc_pos", "sd_ARI"), names_to = "Metric", values_to = "Sd")
  # res_sd <- data.frame(model = res_sd$model,pilot_size = res_sd$pilot_size, epoch = res_sd$epoch,
  #                      Metric = substr(res_sd$Metric, 4, 10), Sd = res_sd$Sd)
  # res_mean <- tidyr::pivot_longer(res, cols = c("mean_ccc_pos", "mean_ARI"), names_to = "Metric", values_to = "Mean")
  # res_mean <- data.frame(model = res_mean$model,pilot_size = res_mean$pilot_size, epoch = res_mean$epoch,
  #                        Metric = substr(res_mean$Metric, 6, 15), Mean = res_mean$Mean)
  # res <- inner_join(res_sd, res_mean, by = c("model","epoch","pilot_size","Metric"))
  # res$Metric <- factor(res$Metric, levels = c("ccc_pos","ARI"),
  #                      labels = c("CCC of PCC", "cARI"))
  # p4 <- ggplot(res, aes(x = pilot_size, y = Mean, group = epoch, color = epoch))+ 
  #   geom_errorbar(aes(ymin = Mean-Sd, ymax = Mean+Sd), width = 4) +
  #   geom_line() + 
  #   geom_point()+
  #   theme_bw()+
  #   facet_nested(Metric~model,scales = "free_y")+
  #   labs(x="Pilot size", y="", color = "Epoch")+
  #   scale_y_continuous(position = "left")+
  #   theme(panel.border = element_rect(colour = "black", fill = NA),
  #         axis.line = element_blank())+
  #   scale_color_npg()
  res <- df_proc(epoch)
  p4 <- ggplot(res)+ 
    geom_violin(aes(x = pilot_size, y = Rate, color = epoch), show.legend = T, position = position_dodge(0.9)) +
    geom_boxplot(aes(x = pilot_size, y = Rate, color = epoch), width = 0.15, position = position_dodge(0.9)) +
    facet_nested(Metric~model)+
    theme_bw()+
    labs(x="Pilot size", y="", color = "Epoch")+
    scale_y_continuous(position = "left")+
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.line = element_blank())+
    scale_color_npg()+
    ylim(0,1)
  if(dataname %in% c("SKCM","BRCA","SKCMLAML")){
    p4 <- p4+geom_hline(data = data.frame(Metric = "1 - Pct(0-markers)",y=0.85), aes(yintercept = y), linetype = "dashed")
  }
  else{
    p4 <- p4+geom_hline(data = data.frame(Metric = "1 - Pct(0-markers)",y=1), aes(yintercept = y), linetype = "dashed")
  }
  pdf(file = paste(default, dataname, "_Epoch.pdf", sep = ""), width = plot_width[5], height = plot_height[5])
  print(gg_default(p4))
  dev.off()
  
  #### p5 normalization violin plot ####
  res <- df_proc(norm)
  p5 <- ggplot(data = res[res$normalization!="DESeq",])+
    geom_violin(aes(x = pilot_size, y = Rate, fill = normalization ), position = position_dodge(0.9))+
    facet_grid(rows = vars(Metric), cols = vars(model),scales = "free_y")+
    theme_bw()+
    labs(x="Pilot size", y="", fill="Normalization")+
    scale_y_continuous(position = "left")+
    theme(panel.border = element_rect(colour = "black", fill = NA),
          axis.line = element_blank())+
    scale_fill_npg()+
    ylim(0,1)
  if(dataname %in% c("SKCM","BRCA","SKCMLAML")){
    p5 <- p5+geom_hline(data = data.frame(Metric = "1 - Pct(0-markers)",y=0.85), aes(yintercept = y), linetype = "dashed")
  }
  else{
    p5 <- p5+geom_hline(data = data.frame(Metric = "1 - Pct(0-markers)",y=1), aes(yintercept = y), linetype = "dashed")
  }
  pdf(file = paste(default, dataname, "_Normalization.pdf", sep = ""), width = plot_width[6], height = plot_height[6])
  print(gg_default(p5))
  dev.off()
  
}
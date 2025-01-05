# SyntheSize-paper-supplementary-materials

This repository stores the data, results, and R scripts to generate these reuslts and figures for the corresponding paper *Optimizing Sample Size for Statistical Learning in Bulk Transcriptome Sequencing: A Learning Curve Approach*.

The method involves two algorithms: SyNG-BTS for augmenting sample sizes using deep generative models and SyntheSize for determing the optimal sample size.

SyNG-BTS is a data augmentation tool synthesizing transcriptomics data with realistic distributions without relying on a predefined formula. Three deep generative models are considered, incluing Variational Auto Encoder (VAE), Generative Adversarial Network (GAN), and flow-based generative model. Those models will be trained on a pilot dataset and then utilized to generate data for any desired number of samples. The workflow of SyNG-BTS is depicted in the following figure:

<p align="center">
  <img src="./pics/syn_bts_workflow.png" width = "1000" alt="method" align=center />
</p>

Original code for SyNG-BTS can be referred to [SyNG-BTS](https://github.com/LXQin/SyNG-BTS). 

SyntheSize is a supervised learning framework designed for determining the optimal sample size by utilizing synthesized data across various sample sizes. This framework employs the inverse power law function (IPLF) to accurately fit augmented data corresponding to different sample sizes and their respective prediction accuracies. SyntheSize also illustrates the generated results through heatmap and UMAP(Uniform Manifold Approximation and Projection).

The workflow of SyntheSize is depicted in the following figure:

<p align="center">
  <img src="./pics/synthesize_workflow.png" width = "1000" alt="method" align=center />
</p>


Original code for SyntheSize can be referred to [SyntheSize](https://github.com/LXQin/SyntheSize). 

## Installation

This repository is *not a package* for SyntheSize. 
It stores the R scripts and data to generate the results and figures in the paper.

## Dependencies

To run the R code, you need to install the following packages:

    ggplot2 3.4.3
    tidyverse 2.0.0 
    DANA 1.1.1
    cowplot 1.1.1 
    ggpubr 0.6.0
    ggsci 3.0.0
    reshape2 1.4.4
    glmnet 4.1-7
    e1071 1.7-13 
    caret 6.0-94
    randomForest 4.7-1.1
    xgboost 1.7.6.1
    ROCR 1.0-11
    class 7.3-22
    
Please make sure to install all dependencies prior to running the code. 
The code presented here was implemented and tested in R version 4.1.1.

## Usage

1. Download this repository.
2. Set your R working directory to the /R/ folder.
3. Check /R/ReproduceExample.Rmd for a step-by-step guide of reproduce the results for one example dataset SKCM with marker filtering.

    This file aims to guide you through the process of reproducing the analyses in the manuscript. In this file, we will use micro RNA sequencing dataset for cancer SKCM with marker filtering as the source/real dataset, and reproduce the evaluation of deep generative models on the pilot datasets.
    
    Before running the code, you should have download this repository. There should be four folders:
    
    * /R/ including this Rmd
    - /RealData/ for saving the downloaded TCGA dataset and store the source datasets for drawing pilot sets.
    - /GeneratedData/ for saving the generated datasets after running SyNG-BTS.
    - /FinalResAug2023/ for saving evaluation results .RData and the visualization figures.
    
    We firstly state using this file:
    
    a. DataOrganization: we firstly download TCGA-SKCM miRNA dataset using TCGAbiolink package. Then filtering out the markers based on mean threshold 4 on log2 scale. There should be a folder ../RealData/ to save the TCGA-SKCM.csv and SKCMPositive_4.csv. This is a part of the code from /R/DataOrganization.Rmd.
    
    b. Using SyNG-BTS PilotExperiments() to run pilot experiments. Go to repository at https://github.com/LXQin/SyNG-BTS. Download the corresponding code, follow the steps in readme. Specify the pilot sample sizes as (20, 60, 100). For each pilot sample size, this function draws 5 random pilot dataset from the SKCMPositive_4 dataset, and for each pilot dataset, this function then runs the specified deep generative models with the specified parameters, and generate new samples with size = 5*sample size of the SKCMPositive_4. Specify the `model` to different deep generative models, `batch_frac` for evaluating the effect of batch size,  `epoch` for different epoch strategies, e.g 1000 for fixed epochs, None for early stopping, `off_aug` for different offline augmentation, "AE_head" for adding auto-encoder before online training. All the generated datasets will be stored in ../GeneratedData/
    
    PilotExperiment(dataname = "SKCMPositive_4", pilot_size = [20, 60, 100],
                    model = "VAE1-10", batch_frac = 0.1, 
                    learning_rate = 0.0005, pre_model = None,
                    epoch = None,  off_aug = None, early_stop_num = 30,
                    AE_head_num = 2, Gaussian_head_num = 9)
                    
                    
    c. Evaluation: with the generated datasets in ../GeneratedData/, we then evaluate the quality of the generated datasets by comparing multiple metrics with the real datasets SKCMPositive_4.csv. This will create .RData file and store it to ../FinalResAug2023/ for visualization. This is a part of the code from /R/Evaluation.Rmd
    
    d. Visualization: with the .RData file, we can visualize the evaluation results. Depending on what generated datasets you have in the ../GeneratedData/ folder, this will generate one of the subfigure in Figure 2. This is a part of the code in /R/Visualization.Rmd

4. If you want to reproduce all results, you can do each step at once

    a. Firstly, go to /R/DataOrganization.Rmd, this will generate all the source datasets for drawing pilot sets, including the raw TCGA datasets, the miRNA and RNA ones, the single-group and two-group ones, the ones with or without marker filtering, and datasets for transfer learning. 
    
    b.Secondly, go do SyNG-BTS, you can use our package or just the code in SyNG-BTS repository. Basically, you will need to run PilotExperiment() multiple times by specifying different parameters to train the deep generative models, and get the generated datasets for further evaluation. Make sure you have enough disk space to save them. For all datasets needed to reproduce the results in our manuscript, make sure you have at least 250GB space to store the generated datasets. 
    
    c. Then, with all the generated datasets in ../GeneratedData/ and source datasets in ../RealData/, you can run /R/Evaluation.Rmd for evaluation results, this will generate all .RData for visualization. 
    
    d. Finally, with all .RData in ../FinalResAug2023/, go to /R/Visualization.Rmd, this will generate all the subfigures for figures in the manuscripts. We also have ../FinalResAug2023/GeneratedFigures/ folder containing these subfigures. You can compare your results with ours. 


5. We also provide the script for the three applications of SyntheSize in /R/CaseStudy.Rmd. All the relevant data files are stored in /Case/
    


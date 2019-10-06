# H3ABioNet 16S rRNA Microbiome Intermediate Bioinformatics Training 2019

## Aim
The H3ABioNet 16S rRNA Microbiome Intermediate Bioinformatics course will provide
training to enable participants to gain the knowledge and skills to perform 16S rRNA
microbiome data analyses using a variety of bioinformatics methods and tools.

## Intended Audience
Anyone who will be working with 16S rRNA microbiome data and would like to learn more
about the bioinformatics analyses.

## Course curriculum
* Module 1: Introduction to the Linux command line / intro to R
* Module 2: Introduction to the microbiome and study design – why 16S?
* Module 3: Sample collection, extraction and library prep for 16S NGS analyses
* Module 4: 16S rRNA gene amplicon sequencing bioinformatics pipeline: the theory
* Module 5: 16S analysis pipeline - QC, ASV picking, taxonomic classification and
alignment
* Module 6: Downstream analysis in R - using the packages phyloseq, NMF, vegan,
metagenomeSeq (among others)
* Potential bonus module: Shotgun sequencing    

## Course objectives
By the end of the course participants should be able to:

* Describe the importance of the microbiome and why it should be studied
* Understand how to design 16S rRNA microbiome studies
* Be able to apply basic syntax and operations in R
* Understand the different NGS data types (e.g. MiSeq reads) produced for a 16S rRNA
microbiome study
* Evaluate the quality of NGS sequence reads and samples
* Understand the various bioinformatics tools used for 16S microbiome studies
* Understand the various 16S rRNA bioinformatics pipelines being used to study the
microbiome.
* Apply the H3ABioNet’s 16S rRNA pipeline and understand how to execute this
* Understand the use of workflow languages (Nextflow) and containerized images
(Singularity) to automate analyses
* Analyze 16S rRNA microbiome data and interpret results

## Learning from this resource
This repository contains the trainer's presentations, tutorial scripts as well as instructions on how to set up Nextflow for automation of 16S rRNA pipeline. The trainers' presentation videos are available on H3ABioNet's Youtube channel on [this playlist](https://www.youtube.com/playlist?list=PLcQ0XMykNhCRG4539nSgLtUJkjVrTwOpg).

This resource assumes that you are running Linux with R and RStudio installed. Additionally, the latest version of the following R packages should be installed:

1. dada2
2. DECIPHER
3. phangorn
4. metagenomeSeq
5. vegan
6. ggplot2
7. NMF
8. gridExtra
9. dplyr
10. phyloseq


This training offers a step-by-step dada2 tutorial on how to analyze 16S rRNA microbiome data (covered in Module 5 and 6) as well as automation of the analysis using Nextflow workflow language and containerized images (Module 6 Session 1), before downstream statistical analysis in R.

To run the [step-by-step tutorial](https://github.com/mbbu/intbt_2019/blob/master/Module_5_16S_analysis_pipeline/dada2_pipeline_tutorial.R), you require R and RStudio and the R packages listed above installed. The dada2 tutorial starts in Module 5 and continues in Module 6 on [downstream analysis](https://github.com/mbbu/intbt_2019/blob/master/Module_6_Nextflow_and_Downstream_analysis_in_R/Session_2_Importing_data_into_R/Dog_microbiome_tutorial_dada2.R). The scripts and required data for the tutorials are also linked in the section below.

To run dada2 nextflow pipeline on a computer cluster, follow the instructions outlined on [this document](https://github.com/mbbu/intbt_2019/blob/master/16S-int-bt-software-setup-and-testing-v1.pdf) for software set up and running the pipeline. It involves setting up Singularity, R and RStudio, and Nextflow. For the training at ICIPE, please follow the instructions outlined [here](https://github.com/mbbu/intbt_2019/blob/master/testing_nextflow_pipeline.md). 

### Links to resources referenced in modules

**Module 1: Introduction to Linux and R**

* [Module resource](https://kviljoen.github.io/H3ABioNet_R/)  
* Rstudio [website](https://rstudio.com/)  
* [RStudio Course Material](https://datacarpentry.org/R-ecology-lesson/00-before-we-start.html#knowing_your_way_around_rstudio)  
* [Debugging code in RStudio](https://resources.rstudio.com/wistia-rstudio-essentials-2/rstudioessentialsprogrammingpart2-2)  
* [Introduction to R](https://datacarpentry.org/R-ecology-lesson/01-intro-to-r.html)  
* [Assignment operators in R](https://renkun.me/2014/01/28/difference-between-assignment-operators-in-r/)  
* [code styling](https://style.tidyverse.org)  
* [R Data types and data structures](https://swcarpentry.github.io/r-novice-inflammation/13-supp-data-structures/)  
* [Good practices in scientific computing](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005510)  

**Module 5: 16S analysis pipeline**

* [dada2 tutorial site](https://iallali.github.io/DADA2_pipeline/16SrRNA_DADA2_pipeline.html)  
* [dada2 tutorial R script](https://github.com/mbbu/intbt_2019/blob/master/Module_5_16S_analysis_pipeline/dada2_pipeline_tutorial.R)  
* Required tutorial data: http://web.cbio.uct.ac.za/~gerrit/downloads/dog_stool_full.tgz  
* RefSeq DB: https://zenodo.org/record/3266798/files/RefSeq-RDP16S_v3_May2018.fa.gz  

**Module 6: Downstream analysis in R**

* [Module resource](https://kviljoen.github.io/H3ABioNet_16S/)  
* [Downstream analysis tutorial R script](https://github.com/mbbu/intbt_2019/blob/master/Module_6_Nextflow_and_Downstream_analysis_in_R/Session_2_Importing_data_into_R/Dog_microbiome_tutorial_dada2.R)  


  
  









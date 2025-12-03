# 16S-Amplicon-Variant-Sequence-Processing-Pipeline
The following a pipeline processing 16S amplicon variant sequence (ASV) data. The DADA2 portion of this script will allow us to process amplicon sequencing data to identify and quantify ASVs. The latter portion of this script uses PHYLOSEQ to visualize the relative abundance of our processed amplicon sequencing data.

## Installation Instructions 
# dada2 

This pipeline uses a plethora of dada2 functions to process and quantify amplicon sequence variants. Installation instructions can be found here: 

https://benjjneb.github.io/dada2/dada-installation.html

The following script was be pasted into the console to download dada2 for this pipeline. 

```{r}
if (!requireNamespace("BiocManager", quietly = TRUE))
install.packages("BiocManager")
BiocManager::install("dada2")
```
Be sure to load your library prior to beginning work. 

```{r}
library(dada2); packageVersion("dada2")
```
# phyloseq

```{r}
# To begin, we will want to install PHYLOSEQ and a set of additional packages (stated below) if they have not already been installed. Do this by adding the below script to the console. 
# install.packages("phyloseq")
# install.packages("Biostrings")
# install.packages("ggplot2")
# install.packages("RColorBrewer")
# install.packages("tidyverse")

# In the even that you receive the warning message: "Package ‘phyloseq’ is not available for this version of R" you can paste the following code into the console to download PHYLOSEQ: 
# if (!require("BiocManager", quietly = TRUE))
    # install.packages("BiocManager")
# BiocManager::install("phyloseq")
```

## Usage 

Examples or demonstrations of how to use the project's features and functionalities.
Code snippets or commands illustrating common use cases.

1. Experimental design!** most important 
2. Quality Control + processing raw reads 
3. Alignment – align reads
4. Identify ASVs (note on OTUs)  - ASVs = # of different types of microbes, you can set an OTU to be 100% and this would make it an ASV – they are different in the sense that an OTU does have 100% percent identity
5. Diversity Indices – how many organisms, what are the relative abundances
6. Classification  ( what is X?) – you don’t need diversity indices for this
7.  Phylogeny ( how related is X to others?)  Not covering in our course
<img width="2105" height="266" alt="image" src="https://github.com/user-attachments/assets/52867775-9374-411e-aa75-2b2edd3e3307" />


# Manual Steps

## Features:
A list or description of the main functionalities and capabilities of the project.

## Technologies Used:
A list of programming languages, frameworks, libraries, and other tools utilized in the project.

## Known Issues or Limitations:
A section detailing any known bugs, limitations, or areas for improvement.

## License:
Information about the project's licensing, specifying how others can use and distribute the code.

## Contact Information and Acknowledgments:
Details on how to contact the project maintainers or authors.
Acknowledgments for any external contributions, resources, or inspirations.

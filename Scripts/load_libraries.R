# Install required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos='http://cran.us.r-project.org' )
BiocManager::install(update = FALSE)

if (!require("stringr", quietly = TRUE))
  install.packages("stringr", repos='http://cran.us.r-project.org')
if (!require("dplyr", quietly = TRUE))
  install.packages("dplyr", repos='http://cran.us.r-project.org')
if (!require("data.table", quietly = TRUE))
  install.packages("data.table", repos='http://cran.us.r-project.org')
if (!require("tidyr", quietly = TRUE))
  install.packages("tidyr", repos='http://cran.us.r-project.org')

# Package names
packages <- c("rtracklayer", "Rsubread", "GenomicRanges")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  BiocManager::install(packages[!installed_packages])
}

# Packages loading
invisible(lapply(packages, library, character.only = TRUE))


# set libraries and standard plotting functions that are used during the analysis

#...................................libraries
packages <- c("ggeconodist", "tidyverse", "here", "ggthemes", "gganimate","writexl",
              "data.table", "ggExtra", "Rsamtools", "GenomicAlignments", "UpSetR",
              "seqTools", "Rsubread", "ape", "DT", "ggpubr", "ggridges", "ggsci","LncFinder",
              "CoverageView", "gghalves", "pryr", "fst", "R.utils", "readxl", "ggseqlogo",
              "GenomicRanges", "magrittr", "dplyr", "fread", "data.table", "stringr", "dplyr", "BiocManager")

invisible(lapply(packages, require, character.only = TRUE))


#### set up parallel processing  #####
#automatic install of packages if they are not installed already
list.of.packages <- c(
  "foreach",
  "doParallel",
  "ranger",
  "doSNOW"
)

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]

if(length(new.packages) > 0){
  install.packages(new.packages, dep=TRUE, repos='http://cran.us.r-project.org')
}

#loading packages
for(package.i in list.of.packages){
  suppressPackageStartupMessages(
    library(
      package.i, 
      character.only = TRUE
    )
  )
}


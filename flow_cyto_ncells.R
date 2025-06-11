# boncat mixtures
# flow cyto 
# n cells
# Jennifer E Harris
# jun 2025

# install packages
BiocManager::install("flowCore")
library(flowCore)

# navigate to place
wd <- "C:/Users/Jenn/The Pennsylvania State University/Burghardt, Liana T - Burghardt Lab Shared Folder/Projects/BONCAT-MicrobialActivity/BONCAT_mixtures/Data/flow_cyto/astrios/JENN HARRIS ALL/2023-6-15 JENN HARRIS"
setwd(wd)


# get list of files
files<-list.files( )
files


###### example ##
## a sample file
fcsFile <- system.file("extdata", "0877408774.B08", package="flowCore")

## read file and linearize values
samp <-  read.FCS(fcsFile, transformation="linearize")
exprs(samp[1:3,])
keyword(samp)[3:6]
class(samp)


## Only read in lines 2 to 5
subset <- read.FCS(fcsFile, which.lines=2:5, transformation="linearize")
exprs(subset)

## Read in a random sample of 100 lines
subset <- read.FCS(fcsFile, which.lines=100, transformation="linearize")
nrow(subset)

#manually supply the alias vs channel options mapping as a data.frame
map <- data.frame(alias = c("A", "B")
                  , channels = c("FL2", "FL4")
)
fr <- read.FCS(fcsFile, channel_alias = map)
fr


######my file##


sample1<-read.FCS(files[1])
sample1@exprs #channels

# not useful params 
sample1@parameters@varMetadata 
sample1@parameters@data #columns
sample1@parameters@.__classVersion__
sample1@parameters@dimLabels


# a just of metadata?
sample1@description





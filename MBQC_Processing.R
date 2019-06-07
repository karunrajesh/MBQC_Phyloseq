set.seed(12341)
rm(list=ls())

#packages for salter data
library(data.table)
library(decontam)
library(phyloseq)
library(reshape2)#;packageVersion("reshape2")
library(xtable)
library(knitr)
options(kableExtra.latex.load_packages = FALSE)
library(kableExtra)
library(DT)
library(Matrix)
library(devtools)
require(ggplot2)
require(sn)
require(fitdistrplus)
require(psych) 
library(PERFect)
library(dirmult)
library(HMP)
library(knitr)
library(png)
library(gridExtra)
library(grid)
library(reshape2)
library(zoo)


#setwd(here::here())
setwd("~/Dropbox/MicrobiomeFiltering/Code/Karun")
#setwd("~/Documents/Katia Lab/MockCommunity/Data")

path.ampli <- "~/Documents/Katia Lab/MockCommunity/Data"
#path.ampli <- "~/Dropbox/MBWG/MockCommunity/Data/"

options(xtable.floating = FALSE)
options(xtable.timestamp = "")
# Read in OTU table
tf.ampli <- read.table(file =  file.path(path.ampli, 'mbqc_integrated_otus.tsv'), sep = '\t', header = TRUE)
view.ampli <- tf.ampli[0:80,0:80]
rownames(view.ampli) <- make.names(view.ampli[,"sample"], unique = TRUE)
# Get the names of all samples
name.ampli <- tf.ampli[0,]
# Get the metadata of the samples, make first row the column names
stat.ampli <- t(tf.ampli[0:71,])
colnames(stat.ampli) = stat.ampli[1,]
stat.ampli = stat.ampli[-1,]

st.ampli <- t(view.ampli[match("Bioinformatics.ID",rownames(view.ampli)):match("thaw_method_sequencing",rownames(view.ampli)),])
colnames(st.ampli) = st.ampli["sample",]
st.ampli = st.ampli[-1,]

# Additional Columns for Statistics
hand.ampli <- read.csv(file =  file.path(path.ampli, 'mbqc_handling_protocols.csv'), header = TRUE, colClasses = "character")
hand.ampli <- hand.ampli[,-1]
other <- colnames(hand.ampli)
sampRow <- rownames(stat.ampli)
x <- matrix(, nrow = length(sampRow), ncol = length(other),dimnames = list(sampRow, other))
labs <- levels(unique(stat.ampli[,'blinded_lab']))
full <- list(length(labs))
for (i in 1:length(labs)){
  inb <- which(stat.ampli[,'blinded_lab']==labs[i])
  full[[i]] <- inb
}
hand.ampli <- as.data.frame(hand.ampli,stringsAsFactors=FALSE)
x <- as.data.frame(x,stringsAsFactors=FALSE)
for(a in 1:length(fullList)){
  for (i in fullList[a]){
    x[i,] <- hand.ampli[a,]
  }
}

stat.ampli <- cbind(data.frame(stat.ampli),x)

# Get the length of the full OTU table
a = dim(tf.ampli)[1]
# Number of taxa in OTU table
begin = match("thaw_method_sequencing",rownames(view.ampli)) + 1
check = (a-begin)+1
# Get all taxa in sample
taxa.ampli <- tf.ampli[begin:a,]
#Dataframe with taxa and sample names
new.ampli <- rbind(name.ampli,taxa.ampli)
rownames(new.ampli) <- new.ampli$sample
new.ampli$sample <- NULL
#stat.ampli$sample <- NULL
otumat <- data.matrix(new.ampli)


oh <- array(1:check)
king <- array(1:check)
phyl <- array(1:check)
cls <- array(1:check)
ord <- array(1:check)
fam <- array(1:check)
gen <- array(1:check)
for (i in 1:check){ 
  oh[i] <- toString(taxa.ampli$sample[i])
  b <- unlist(strsplit(oh[i],"__"))
  king[i] <- substr(b[2],1,nchar(b[2])-2) #kingdom
  phyl[i] <- substr(b[3],1,nchar(b[3])-2) #phylum
  cls[i] <- substr(b[4],1,nchar(b[4])-2)  #class
  ord[i] <- substr(b[5],1,nchar(b[5])-2)  #order
  fam[i] <- substr(b[6],1,nchar(b[6])-2)  #family
  gen[i] <- substr(b[7],1,nchar(b[7])-2)  #genus
}
#Table of taxa with classification
taxmat = matrix(c(king,phyl,cls,ord,fam,gen), nrow = check, ncol = 6)
rownames(taxmat) <- rownames(new.ampli)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")
#OTU Table
OTU = otu_table(otumat, taxa_are_rows = TRUE)
#TAXA Table
TAX = tax_table(taxmat)
write.csv(taxmat, "TaxaTable.csv")
#Sample Metadata
SAMPLE = sample_data(data.frame(stat.ampli))
#Phyloseq Object
physeq = phyloseq(OTU, TAX, SAMPLE)


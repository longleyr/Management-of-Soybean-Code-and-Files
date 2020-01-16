#Title: Crop Management Impacts the Soy (Glycine max) Microbiome
#Authors: Reid Longley, Zachary A. Noel, Gian Maria Niccolo Benucci, Martin I. Chilvers, Frances Trail, Gregory Bonito 
#Journal: Submitted to Frontiers in Microbiology

# Fungal Script - Reid Longley
set.seed(9289)
# load packages ---------------------

library(stringi)
library(rhdf5)
library(zlibbioc)
library(S4Vectors)
library(phyloseq)
library(Biostrings)
library(yaml)
library(colorspace)
library(ggplot2)
library(indicspecies)
library(vegan)
library(phyloseq)


# create a phyloseq object ----------
ITS_otus<- read.delim("otu_table_Total_UPARSE_R1.txt",
                       row.names=1) 
head(ITS_otus)


ITS_otus_phy <-otu_table(ITS_otus,
                          taxa_are_rows = TRUE)
ITS_otus_phy


ITS_metadata <-read.delim("fungi_total_mapping.txt",
                           row.names=1)
ITS_metadata
ITS_metadata_phy <-sample_data(ITS_metadata)

ITS_taxonomy<-read.delim("total_fungal_tax.txt",
                          header=TRUE, 
                          row.names=1)
ITS_taxonomy
ITS_taxonomy_phy <- tax_table(as.matrix(ITS_taxonomy))


ITS_sequences <- readDNAStringSet("otus_total_R1.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
ITS_sequences

physeq_object <- phyloseq(ITS_otus_phy,
                          ITS_metadata_phy,
                          ITS_taxonomy_phy,
                          ITS_sequences)
physeq_object
tax_table(physeq_object)
sample_data(physeq_object)

###formatting taxonomy------------------------
colnames(tax_table(physeq_object)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
tax_table(physeq_object)

tax_table(physeq_object)[, "Kingdom"] <- gsub("d:", "", tax_table(physeq_object)[, "Kingdom"])
tax_table(physeq_object)[, "Phylum"] <- gsub("p:", "", tax_table(physeq_object)[, "Phylum"])
tax_table(physeq_object)[, "Class"] <- gsub("c:", "", tax_table(physeq_object)[, "Class"])
tax_table(physeq_object)[, "Order"] <- gsub("o:", "", tax_table(physeq_object)[, "Order"])
tax_table(physeq_object)[, "Family"] <- gsub("f:", "", tax_table(physeq_object)[, "Family"])
tax_table(physeq_object)[, "Genus"] <- gsub("g:", "", tax_table(physeq_object)[, "Genus"])
tax_table(physeq_object)[, "Species"] <- gsub("s:", "", tax_table(physeq_object)[, "Species"])
tax_table(physeq_object)


#remove chloropast,mitochondria,cyanobacteria
physeq_object <- subset_taxa(physeq_object, Phylum!="Chloroplast")
physeq_object <- subset_taxa(physeq_object, Class!="Chloroplast")
physeq_object <- subset_taxa(physeq_object, Order!="Chloroplast")
physeq_object <- subset_taxa(physeq_object, Family!="Chloroplast")
physeq_object <- subset_taxa(physeq_object, Genus!="Chloroplast")
tax_table(physeq_object)

physeq_object <- subset_taxa(physeq_object, Phylum!="Mitochondria")
physeq_object <- subset_taxa(physeq_object, Class!="Mitochondria")
physeq_object <- subset_taxa(physeq_object, Order!="Mitochondria")
physeq_object <- subset_taxa(physeq_object, Family!="Mitochondria")
physeq_object <- subset_taxa(physeq_object, Genus!="Mitochondria")
tax_table(physeq_object)





###remove contaminants--------------------------------------------------



#removing samples from soil dataset which were noted to not dry properly, and had mold growing in envelope. were analyzed, but only produced a few otus

otu_table(physeq_object) <- subset(otu_table(physeq_object),
                                       select = -c(T4R5AR2S,T4R6BR2S,T4R2CR2S,T4R5CR2S,T2R6BR2S))




# must split by sample source so that contaminants are removed from each individual run
physeq_object_stems <- subset_samples(physeq_object, origin%in%c("stem"))
physeq_object_leaves <- subset_samples(physeq_object,origin%in%c("leaves"))
physeq_object_roots <- subset_samples(physeq_object,origin%in%c("root"))
physeq_object_soil <- subset_samples(physeq_object,origin%in%c("soil"))
sample_data(physeq_object_soil)
library(devtools)
library(processx)
devtools::install_github("benjjneb/decontam", force = TRUE)

# Install decontam package which will be used to remove potential contaminants 
# by checking abundance distributions in negative controls vs regular samples
library(decontam)


#soil
# check library size distribution
write.csv(sample_data(physeq_object_soil), file = "sample_check1_soil.csv")
df_soil <- as.data.frame(sample_data(physeq_object_soil)) # Put sample_data into a ggplot-friendly data.frame
df_soil$LibrarySize_soil <- sample_sums(physeq_object_soil)
df_soil <- df_soil[order(df_soil$LibrarySize_soil),]
df_soil$Index <- seq(nrow(df_soil))
write.csv(df_soil, file = "rank_sums_soil.csv")
ggplot(data=df_soil, aes(x=Index, y=LibrarySize_soil, color=Sample_or_Control)) + geom_point()

# filter by prevelance 
sample_data(physeq_object_soil)$is.neg <- sample_data(physeq_object_soil)$Sample_or_Control == "Control Sample"

contamdf.prev_soil <- isContaminant(physeq_object_soil, method="prevalence", neg="is.neg")
table(contamdf.prev_soil$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_soil <- transform_sample_counts(physeq_object_soil, function(abund) 1*(abund>0))
ps.pa.neg_soil <- prune_samples(sample_data(ps.pa_soil)$Sample_or_Control == "Control Sample", ps.pa_soil)
ps.pa.pos_soil <- prune_samples(sample_data(ps.pa_soil)$Sample_or_Control == "True Sample", ps.pa_soil)
# Make data.frame of prevalence in positive and negative samples
df.pa_soil <- data.frame(pa.pos_soil=taxa_sums(ps.pa.pos_soil), pa.neg_soil=taxa_sums(ps.pa.neg_soil),
                          contaminant=contamdf.prev_soil$contaminant)
ggplot(data=df.pa_soil, aes(x=pa.neg_soil, y=pa.pos_soil, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_soil <- prune_taxa(!contamdf.prev_soil$contaminant, physeq_object_soil)
# with contaminants removed
otu_table(ps.noncontam_soil)



#roots


# check library size distribution
write.csv(sample_data(physeq_object_roots), file = "sample_check1_roots.csv")
df_roots <- as.data.frame(sample_data(physeq_object_roots)) # Put sample_data into a ggplot-friendly data.frame
df_roots$LibrarySize_roots <- sample_sums(physeq_object_roots)
df_roots <- df_roots[order(df_roots$LibrarySize_roots),]
df_roots$Index <- seq(nrow(df_roots))
write.csv(df_roots, file = "rank_sums_roots.csv")
ggplot(data=df_roots, aes(x=Index, y=LibrarySize_roots, color=Sample_or_Control)) + geom_point()

# filter by prevelance 
sample_data(physeq_object_roots)$is.neg <- sample_data(physeq_object_roots)$Sample_or_Control == "Control Sample"
contamdf.prev_roots <- isContaminant(physeq_object_roots, method="prevalence", neg="is.neg")
table(contamdf.prev_roots$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_roots <- transform_sample_counts(physeq_object_roots, function(abund) 1*(abund>0))
ps.pa.neg_roots <- prune_samples(sample_data(ps.pa_roots)$Sample_or_Control == "Control Sample", ps.pa_roots)
ps.pa.pos_roots <- prune_samples(sample_data(ps.pa_roots)$Sample_or_Control == "True Sample", ps.pa_roots)
# Make data.frame of prevalence in positive and negative samples
df.pa_roots <- data.frame(pa.pos_roots=taxa_sums(ps.pa.pos_roots), pa.neg_roots=taxa_sums(ps.pa.neg_roots),
                          contaminant=contamdf.prev_roots$contaminant)
ggplot(data=df.pa_roots, aes(x=pa.neg_roots, y=pa.pos_roots, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_roots <- prune_taxa(!contamdf.prev_roots$contaminant, physeq_object_roots)
# with contaminants removed
otu_table(ps.noncontam_roots)

#stems

# check library size distribution
df_stems <- as.data.frame(sample_data(physeq_object_stems)) # Put sample_data into a ggplot-friendly data.frame
df_stems$LibrarySize_stems <- sample_sums(physeq_object_stems)
df_stems <- df_stems[order(df_stems$LibrarySize_stems),]
df_stems$Index <- seq(nrow(df_stems))
write.csv(df_stems, file = "rank_sums_stems.csv")
ggplot(data=df_stems, aes(x=Index, y=LibrarySize_stems, color=Sample_or_Control)) + geom_point()

# filter by prevelance 
sample_data(physeq_object_stems)$is.neg <- sample_data(physeq_object_stems)$Sample_or_Control == "Control Sample"
contamdf.prev_stems <- isContaminant(physeq_object_stems, method="prevalence", neg="is.neg")
table(contamdf.prev_stems$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_stems <- transform_sample_counts(physeq_object_stems, function(abund) 1*(abund>0))
ps.pa.neg_stems <- prune_samples(sample_data(ps.pa_stems)$Sample_or_Control == "Control Sample", ps.pa_stems)
ps.pa.pos_stems <- prune_samples(sample_data(ps.pa_stems)$Sample_or_Control == "True Sample", ps.pa_stems)
# Make data.frame of prevalence in positive and negative samples
df.pa_stems <- data.frame(pa.pos_stems=taxa_sums(ps.pa.pos_stems), pa.neg_stems=taxa_sums(ps.pa.neg_stems),
                          contaminant=contamdf.prev_stems$contaminant)
ggplot(data=df.pa_stems, aes(x=pa.neg_stems, y=pa.pos_stems, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_stems <- prune_taxa(!contamdf.prev_stems$contaminant, physeq_object_stems)
ps.noncontam_stems# with contaminants removed
otu_table(ps.noncontam_stems)

# leaves

# check library size distribution
df_leaves <- as.data.frame(sample_data(physeq_object_leaves)) # Put sample_data into a ggplot-friendly data.frame
df_leaves$LibrarySize_leaves <- sample_sums(physeq_object_leaves)
df_leaves <- df_leaves[order(df_leaves$LibrarySize_leaves),]
df_leaves$Index <- seq(nrow(df_leaves))
write.csv(df_leaves, file = "rank_sums_leaves.csv")
ggplot(data=df_leaves, aes(x=Index, y=LibrarySize_leaves, color=Sample_or_Control)) + geom_point()

# filter by prevelance 
sample_data(physeq_object_leaves)$is.neg <- sample_data(physeq_object_leaves)$Sample_or_Control == "Control Sample"
contamdf.prev_leaves <- isContaminant(physeq_object_leaves, method="prevalence", neg="is.neg")
table(contamdf.prev_leaves$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_leaves <- transform_sample_counts(physeq_object_leaves, function(abund) 1*(abund>0))
ps.pa.neg_leaves <- prune_samples(sample_data(ps.pa_leaves)$Sample_or_Control == "Control Sample", ps.pa_leaves)
ps.pa.pos_leaves <- prune_samples(sample_data(ps.pa_leaves)$Sample_or_Control == "True Sample", ps.pa_leaves)
# Make data.frame of prevalence in positive and negative samples
df.pa_leaves <- data.frame(pa.pos_leaves=taxa_sums(ps.pa.pos_leaves), pa.neg_leaves=taxa_sums(ps.pa.neg_leaves),
                         contaminant=contamdf.prev_leaves$contaminant)
ggplot(data=df.pa_leaves, aes(x=pa.neg_leaves, y=pa.pos_leaves, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_leaves <- prune_taxa(!contamdf.prev_leaves$contaminant, physeq_object_leaves)
ps.noncontam_leaves# with contaminants removed
otu_table(ps.noncontam_leaves)

# remove negative controls from soil
ps.noncontam_soil <- subset_samples(ps.noncontam_soil, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_soil) <- otu_table(ps.noncontam_soil)[which(rowSums(otu_table(ps.noncontam_soil)) >= 1),]
ps.noncontam_soil
otu_table(ps.noncontam_soil)
sample_data(ps.noncontam_soil)
ps.noncontam_soil

# export soil otu table to check samples
# Following removal of contaminants identified by the decontam package, removing any samples which had 
# less than 1000 reads to avoid biasing beta diversity analyses
# removing samples with less than 1000 reads, including samples that had less than 1000 read following sampling
# using csv files to check
write.csv(otu_table(ps.noncontam_soil), file = "filtering_low_soil.csv")


otu_table(ps.noncontam_soil) <- subset(otu_table(ps.noncontam_soil),
                                       select = -c(T4R1CR6S,T1R1BR6S,T4R2CR6S,T4R6AR2S,T1R1FBR3S,T1R2CCR3S,T1R1AR6S,T1R1FCR3S,T2R2FBR4S,T2R2FCR4S,T1R5CR2S,T1R6CCR4S,T2R5CR2S))


# export roots otu table to check samples

# remove negative controls
ps.noncontam_roots <- subset_samples(ps.noncontam_roots, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_roots) <- otu_table(ps.noncontam_roots)[which(rowSums(otu_table(ps.noncontam_roots)) >= 1),]
ps.noncontam_roots



write.csv(otu_table(ps.noncontam_roots), file = "filtering_low_roots.csv")


otu_table(ps.noncontam_roots) <- subset(otu_table(ps.noncontam_roots),
                                        select = -c(T1R5CBR3R,T2R1CBR6R,T4R5CR2R,T4R2BR2R,T4R1AR2R,T1R6FCR4R,T2R5AR2R,T1R1BR6R,T1R2FCR4R,T1R2AR2R,T1R2FAR4R,T1R6CCR4R,T1R2AR2R,T1R5FAR3R,T1R1FAR3R,T1R6CR2R,T1R2FAR3R,T4R6BR2R))
ps.noncontam_roots


# leaves
# remove negative controls
ps.noncontam_leaves <- subset_samples(ps.noncontam_leaves, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_leaves) <- otu_table(ps.noncontam_leaves)[which(rowSums(otu_table(ps.noncontam_leaves)) >= 1),]
ps.noncontam_leaves


# export leaves otu table to check samples
write.csv(otu_table(ps.noncontam_leaves), file = "filtering_low_leaves.csv")


otu_table(ps.noncontam_leaves) <- subset(otu_table(ps.noncontam_leaves),
                                         select = -c(T2R1FCR6L))
sample_data(ps.noncontam_leaves)
#stems
ps.noncontam_stems <- subset_samples(ps.noncontam_stems, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_stems) <- otu_table(ps.noncontam_stems)[which(rowSums(otu_table(ps.noncontam_stems)) >= 1),]
ps.noncontam_stems



# export stems otu table to check samples
write.csv(otu_table(ps.noncontam_stems), file = "filtering_low_stems.csv")


otu_table(ps.noncontam_stems) <- subset(otu_table(ps.noncontam_stems),
                                        select = -c(T4R1CR2ST,T4R6CR2ST,T1R1AR4ST,T1R6FAR4ST,T2R5FCR3ST,T2R6AR2ST,T2R2CBR6ST,T1R6CR2ST,T1R6FBR6ST,T4R5AR2ST,T1R2AR2ST,T4R6CR6ST,T4R5CR2ST,T1R6FCR3ST,T2R5CR6ST,T2R6FBR3ST,T2R6FCR6ST,T4R5BR6ST,T1R2FBR4ST,T4R1BR6ST,T4R6CR2ST,T4R1AR6ST,T4R6AR6ST,T2R5CBR6ST))
sample_data(ps.noncontam_stems)

###for alpha diversity need data which includes singletons etc. will split the current datasets by objective in order to get read distributions for each objective-----------------------

# checking read distributions to see how sequencing depth is distributed
# soil read distributions
ps.noncontam_soil_obj1 <- subset_samples(ps.noncontam_soil, experiment%in%c("obj_1"))
sample_data(ps.noncontam_soil_obj1)


sums_soil_obj1 <- data.frame(colSums(otu_table(ps.noncontam_soil_obj1)))
colnames(sums_soil_obj1) <- "Sample_totalSeqs_soil"
sums_soil_obj1$Sample <- row.names(sums_soil_obj1)
sums_soil_obj1
sums_roots_obj1


ggplot(sums_soil_obj1, aes(x=Sample_totalSeqs_soil)) +
  geom_histogram(binwidth=500, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sample_totalSeqs_soil, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)

# root read distributions

ps.noncontam_roots_obj1 <- subset_samples(ps.noncontam_roots, experiment%in%c("obj_1"))



sums_roots_obj1 <- data.frame(colSums(otu_table(ps.noncontam_roots_obj1)))
colnames(sums_roots_obj1) <- "Sample_totalSeqs_roots"
sums_roots_obj1$Sample <- row.names(sums_roots_obj1)



ggplot(sums_roots_obj1, aes(x=Sample_totalSeqs_roots)) +
  geom_histogram(binwidth=500, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sample_totalSeqs_roots, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)


# stem read distributions

ps.noncontam_stems_obj1 <- subset_samples(ps.noncontam_stems, experiment%in%c("obj_1"))



sums_stems_obj1 <- data.frame(colSums(otu_table(ps.noncontam_stems_obj1)))
colnames(sums_stems_obj1) <- "Sample_totalSeqs_stems"
sums_stems_obj1$Sample <- row.names(sums_stems_obj1)



ggplot(sums_stems_obj1, aes(x=Sample_totalSeqs_stems)) +
  geom_histogram(binwidth=500, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sample_totalSeqs_stems, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)

# leaves read distributions

ps.noncontam_leaves_obj1 <- subset_samples(ps.noncontam_leaves, experiment%in%c("obj_1"))



sums_leaves_obj1 <- data.frame(colSums(otu_table(ps.noncontam_leaves_obj1)))
colnames(sums_leaves_obj1) <- "Sample_totalSeqs_leaves"
sums_leaves_obj1$Sample <- row.names(sums_leaves_obj1)



ggplot(sums_leaves_obj1, aes(x=Sample_totalSeqs_leaves)) +
  geom_histogram(binwidth=500, colour="black", fill="white") +
  geom_vline(aes(xintercept=mean(Sample_totalSeqs_leaves, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)

# plotting rarefaction curves object 1 
# rarefaction wasn't used, but plotted curves to see how well communities were sampled

# soil
otu_rare_soil <- as.data.frame(otu_table(ps.noncontam_soil_obj1))
metadata_soil_rare <- as.data.frame(sample_data(ps.noncontam_soil_obj1))
rarecurve(t(otu_rare_soil), col = metadata_soil_rare$Management, label = FALSE, 
          step = 50,
          main="Fungi_soil", ylab = "Number of OTUs", xlab = "Number of DNA reads") -> rare_soil_fungi
legend("bottomright", legend=c("No-Till", "organic", "traditional"),
       col=c("black", "green", "red"), lty=1, cex=0.8, box.lty=1) 
# roots
otu_rare_roots <- as.data.frame(otu_table(ps.noncontam_roots_obj1))
metadata_roots_rare <- as.data.frame(sample_data(ps.noncontam_roots_obj1))
rarecurve(t(otu_rare_roots), col = metadata_roots_rare$Management, label = FALSE, 
           step = 50,
          main="Fungi_Roots", ylab = "Number of OTUs", xlab = "Number of DNA reads") -> rare_roots_fungi
legend("bottomright", legend=c("No-Till", "organic", "traditional"),
       col=c("black", "green", "red"), lty=1, cex=0.8, box.lty=1) 

#stems
otu_rare_stems <- as.data.frame(otu_table(ps.noncontam_stems_obj1))
metadata_stems_rare <- as.data.frame(sample_data(ps.noncontam_stems_obj1))
rarecurve(t(otu_rare_stems), col = metadata_stems_rare$Management, label = FALSE, 
          step = 50,
          main="Fungi_stems", ylab = "Number of OTUs", xlab = "Number of DNA reads") -> rare_stems_fungi
legend("bottomright", legend=c("No-Till", "organic", "traditional"),
       col=c("black", "green", "red"), lty=1, cex=0.8, box.lty=1) 


# leaves
otu_rare_leaves <- as.data.frame(otu_table(ps.noncontam_leaves_obj1))
metadata_leaves_rare <- as.data.frame(sample_data(ps.noncontam_leaves_obj1))
rarecurve(t(otu_rare_leaves), col = metadata_leaves_rare$Management, label = FALSE, 
          step = 50,
          main="Fungi_leaves", ylab = "Number of OTUs", xlab = "Number of DNA reads") -> rare_leaves_fungi
legend("bottomright", legend=c("No-Till", "organic", "traditional"),
       col=c("black", "green", "red"), lty=1, cex=0.8, box.lty=1) 

###alpha diversity-----------------------------------
# used the non normalized data as reccomended by phyloseq 
# singletons not removed at this step
#Table 1 fungal alpha diversity
install.packages("BiodiversityR", dependencies = TRUE)
install.packages(pkgs=c("BiodiversityR", "vegan",
                        "Rcmdr", "MASS", "mgcv",
                        "cluster", "RODBC", "rpart", "effects", "multcomp",
                        "ellipse", "maptree", "sp", "splancs", "spatial",
                        "akima", "nnet", "dismo", "raster", "rgdal",
                        "bootstrap", "PresenceAbsence",
                        "maxlike", "gbm", "randomForest", "gam", "earth", "mda",
                        "kernlab", "e1071", "glmnet", "sem", "rgl", "relimp",
                        "lmtest", "leaps", "Hmisc", "colorspace", "aplpack",
                        "abind", "XLConnect", "car", "markdown", "knitr",
                        "geosphere", "maptools", "rgeos", "ENMeval", "red"),
                 dependencies=c("Depends", "Imports"))
install.packages("rgl")
library(rgl)
library(BiodiversityR)

# make data frames of metadata and otu tables for analyzing alpha diversity 
otu_soil_fungi <- as.data.frame(otu_table(ps.noncontam_soil_obj1))
meta_soil_fungi <- as.data.frame(sample_data(ps.noncontam_soil_obj1))
meta_soil_fungi
otu_roots_fungi <- as.data.frame(otu_table(ps.noncontam_roots_obj1))
meta_roots_fungi <- as.data.frame(sample_data(ps.noncontam_roots_obj1))
otu_stems_fungi <- as.data.frame(otu_table(ps.noncontam_stems_obj1))
meta_stems_fungi <- as.data.frame(sample_data(ps.noncontam_stems_obj1))
otu_leaves_fungi <- as.data.frame(otu_table(ps.noncontam_leaves_obj1))
meta_leaves_fungi <- as.data.frame(sample_data(ps.noncontam_leaves_obj1))
meta_roots_fungi
sample_data(ps.noncontam_soil_obj1)
#soil
alpha_div_soil_fungi <- meta_soil_fungi
alpha_div_soil_fungi
alpha_div_soil_fungi$readNO <- sample_sums(ps.noncontam_soil_obj1)
alpha_div_soil_fungi$Observed <- specnumber(otu_soil_fungi, MARGIN = 2)
alpha_div_soil_fungi$Shannon <- diversity(otu_soil_fungi, index="shannon", MARGIN = 2)
jevenness_fungi <- diversityresult(t(otu_soil_fungi), method = "each site", index = "Jevenness")
alpha_div_soil_fungi$Jevenness <- jevenness_fungi$Jevenness
#alpha_div_soil_fungi <- alpha_div_soil_fungi[order(alpha_div_soil_fungi$ReadNO), ]
alpha_div_soil_fungi



# get descriptive stats
install.packages("psych")
library("psych")
# stats were pulled from this table and put into excel
describeBy(alpha_div_soil_fungi, alpha_div_soil_fungi$Indicator_label)
warnings()
# for this purpose, these warnings can be ignored, the data needed to form the 
# alpha diversity table are produced correctly


#roots
alpha_div_roots_fungi <- meta_roots_fungi
alpha_div_roots_fungi
alpha_div_roots_fungi$readNO <- sample_sums(ps.noncontam_roots_obj1)
alpha_div_roots_fungi$Observed <- specnumber(otu_roots_fungi, MARGIN = 2)
alpha_div_roots_fungi$Shannon <- diversity(otu_roots_fungi, index="shannon", MARGIN = 2)
jevenness_fungi <- diversityresult(t(otu_roots_fungi), method = "each site", index = "Jevenness")
alpha_div_roots_fungi$Jevenness <- jevenness_fungi$Jevenness
#alpha_div_roots_fungi <- alpha_div_roots_fungi[order(alpha_div_roots_fungi$ReadNO), ]
alpha_div_roots_fungi

# get descriptive stats
describeBy(alpha_div_roots_fungi, alpha_div_roots_fungi$Indicator_label)

#stems
alpha_div_stems_fungi <- meta_stems_fungi
alpha_div_stems_fungi
alpha_div_stems_fungi$readNO <- sample_sums(ps.noncontam_stems_obj1)
alpha_div_stems_fungi$Observed <- specnumber(otu_stems_fungi, MARGIN = 2)
alpha_div_stems_fungi$Shannon <- diversity(otu_stems_fungi, index="shannon", MARGIN = 2)
jevenness_fungi <- diversityresult(t(otu_stems_fungi), method = "each site", index = "Jevenness")
alpha_div_stems_fungi$Jevenness <- jevenness_fungi$Jevenness
alpha_div_stems_fungi <- alpha_div_stems_fungi[order(alpha_div_stems_fungi$ReadNO), ]
alpha_div_stems_fungi

# get descriptive stats
describeBy(alpha_div_stems_fungi, alpha_div_stems_fungi$Indicator_label)

#leaves
alpha_div_leaves_fungi <- meta_leaves_fungi
alpha_div_leaves_fungi
alpha_div_leaves_fungi$readNO <- sample_sums(ps.noncontam_leaves_obj1)
alpha_div_leaves_fungi$Observed <- specnumber(otu_leaves_fungi, MARGIN = 2)
alpha_div_leaves_fungi$Shannon <- diversity(otu_leaves_fungi, index="shannon", MARGIN = 2)
jevenness_fungi <- diversityresult(t(otu_leaves_fungi), method = "each site", index = "Jevenness")
alpha_div_leaves_fungi$Jevenness <- jevenness_fungi$Jevenness
#alpha_div_leaves_fungi <- alpha_div_leaves_fungi[order(alpha_div_leaves_fungi$ReadNO), ]
alpha_div_leaves_fungi

# get descriptive stats
describeBy(alpha_div_leaves_fungi, alpha_div_leaves_fungi$Indicator_label)
alpha_div_leaves_fungi
#soil - this fixes the dataframes
write.csv(alpha_div_soil_fungi, file = "alpha_soil_check.csv")
alpha_div_soil_fungi <- read.csv("alpha_soil_check.csv")
write.csv(alpha_div_roots_fungi, file = "alpha_roots_check.csv")
alpha_div_roots_fungi <- read.csv("alpha_roots_check.csv")
write.csv(alpha_div_stems_fungi, file = "alpha_stems_check.csv")
alpha_div_stems_fungi <- read.csv("alpha_stems_check.csv")
write.csv(alpha_div_leaves_fungi, file = "alpha_leaves_check.csv")
alpha_div_leaves_fungi<- read.csv("alpha_leaves_check.csv")
alpha


#soil

# get significant differences
install.packages("agricolae")
library("agricolae")
alpha_div_soil_fungi
library(dplyr)
alpha_div_soil_fungi_df
# make data frames split by growth stage
alpha_div_soil_fungi
alpha_div_soil_fungi_df_V2 <- alpha_div_soil_fungi[alpha_div_soil_fungi$Growth_Stage == "V2",]
alpha_div_soil_fungi_df_R2 <- alpha_div_soil_fungi[alpha_div_soil_fungi$Growth_Stage == "R2",]
alpha_div_soil_fungi_df_R2 
alpha_div_soil_fungi_df_R6 <- alpha_div_soil_fungi[alpha_div_soil_fungi$Growth_Stage == "R6",]
alpha_div_soil_fungi_df_R6



#soil V2 rich
aov_fungi_soil_V2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_soil_fungi_df_V2)
summary(aov_fungi_soil_V2_rich)
HSD.test(aov_fungi_soil_V2_rich, "Indicator_label") -> tukeyHSD_fungi_soil_V2_rich
tukeyHSD_fungi_soil_V2_rich
#soil R2 rich
aov_fungi_soil_R2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_soil_fungi_df_R2)
summary(aov_fungi_soil_R2_rich)
HSD.test(aov_fungi_soil_R2_rich, "Indicator_label") -> tukeyHSD_fungi_soil_R2_rich
tukeyHSD_fungi_soil_R2_rich
# soil R6 rich
aov_fungi_soil_R6_rich <- aov(Observed ~ Indicator_label, data=alpha_div_soil_fungi_df_R6)
summary(aov_fungi_soil_R6_rich)
HSD.test(aov_fungi_soil_R6_rich, "Indicator_label") -> tukeyHSD_fungi_soil_R6_rich
tukeyHSD_fungi_soil_R6_rich


#soil V2 shan
aov_fungi_soil_V2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_soil_fungi_df_V2)
summary(aov_fungi_soil_V2_shan)
HSD.test(aov_fungi_soil_V2_shan, "Indicator_label") -> tukeyHSD_fungi_soil_V2_shan
tukeyHSD_fungi_soil_V2_shan
#soil R2 shan
aov_fungi_soil_R2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_soil_fungi_df_R2)
summary(aov_fungi_soil_R2_shan)
HSD.test(aov_fungi_soil_R2_shan, "Indicator_label") -> tukeyHSD_fungi_soil_R2_shan
tukeyHSD_fungi_soil_R2_shan
# soil R6 shan
aov_fungi_soil_R6_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_soil_fungi_df_R6)
summary(aov_fungi_soil_R6_shan)
HSD.test(aov_fungi_soil_R6_shan, "Indicator_label") -> tukeyHSD_fungi_soil_R6_shan
tukeyHSD_fungi_soil_R6_shan

#soil V2 even
aov_fungi_soil_V2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_soil_fungi_df_V2)
summary(aov_fungi_soil_V2_even)
HSD.test(aov_fungi_soil_V2_even, "Indicator_label") -> tukeyHSD_fungi_soil_V2_even
tukeyHSD_fungi_soil_V2_even
#soil R2 even
aov_fungi_soil_R2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_soil_fungi_df_R2)
summary(aov_fungi_soil_R2_even)
HSD.test(aov_fungi_soil_R2_even, "Indicator_label") -> tukeyHSD_fungi_soil_R2_even
tukeyHSD_fungi_soil_R2_even
# soil R6 even
alpha_div_soil_fungi_df_R6 <- alpha_div_soil_fungi_df[alpha_div_soil_fungi_df$Growth_Stage == " R6",]
alpha_div_soil_fungi_df_R6
aov_fungi_soil_R6_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_soil_fungi_df_R6)
summary(aov_fungi_soil_R6_even)
HSD.test(aov_fungi_soil_R6_even, "Indicator_label") -> tukeyHSD_fungi_soil_R6_even
tukeyHSD_fungi_soil_R6_even
#must also split into management systems 
alpha_div_soil_fungi_df_Conventional <- alpha_div_soil_fungi[alpha_div_soil_fungi$Management == "Conventional",]
alpha_div_soil_fungi_df_Conventional
alpha_div_soil_fungi_df_No_Till <- alpha_div_soil_fungi[alpha_div_soil_fungi$Management == "No-Till",]
alpha_div_soil_fungi_df_No_Till
alpha_div_soil_fungi_df_Organic <- alpha_div_soil_fungi[alpha_div_soil_fungi$Management == "Organic",]
alpha_div_soil_fungi_df_Organic

#soil Conventional rich
aov_fungi_soil_Conventional_rich <- aov(Observed ~ Indicator_label, data=alpha_div_soil_fungi_df_Conventional)
aov_fungi_soil_Conventional_rich
summary(aov_fungi_soil_Conventional_rich)
HSD.test(aov_fungi_soil_Conventional_rich, "Indicator_label") -> tukeyHSD_fungi_soil_Conventional_rich
tukeyHSD_fungi_soil_Conventional_rich
#soil No-Till rich
aov_fungi_soil_No_Till_rich <- aov(Observed ~ Indicator_label, data=alpha_div_soil_fungi_df_No_Till)
summary(aov_fungi_soil_No_Till_rich)
HSD.test(aov_fungi_soil_No_Till_rich, "Indicator_label") -> tukeyHSD_fungi_soil_No_Till_rich
tukeyHSD_fungi_soil_No_Till_rich
# soil Organic rich
alpha_div_soil_fungi_df_Organic
aov_fungi_soil_Organic_rich <- aov(Observed ~ Indicator_label, data=alpha_div_soil_fungi_df_Organic)
summary(aov_fungi_soil_Organic_rich)
HSD.test(aov_fungi_soil_Organic_rich, "Indicator_label") -> tukeyHSD_fungi_soil_Organic_rich
tukeyHSD_fungi_soil_Organic_rich


#soil Conventional shan
aov_fungi_soil_Conventional_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_soil_fungi_df_Conventional)
summary(aov_fungi_soil_Conventional_shan)
HSD.test(aov_fungi_soil_Conventional_shan, "Indicator_label") -> tukeyHSD_fungi_soil_Conventional_shan
tukeyHSD_fungi_soil_Conventional_shan
#soil No-Till shan
aov_fungi_soil_No_Till_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_soil_fungi_df_No_Till)
summary(aov_fungi_soil_No_Till_shan)
HSD.test(aov_fungi_soil_No_Till_shan, "Indicator_label") -> tukeyHSD_fungi_soil_No_Till_shan
tukeyHSD_fungi_soil_No_Till_shan
# soil Organic shan
alpha_div_soil_fungi_df_Organic
aov_fungi_soil_Organic_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_soil_fungi_df_Organic)
summary(aov_fungi_soil_Organic_shan)
HSD.test(aov_fungi_soil_Organic_shan, "Indicator_label") -> tukeyHSD_fungi_soil_Organic_shan
tukeyHSD_fungi_soil_Organic_shan

#soil Conventional even
aov_fungi_soil_Conventional_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_soil_fungi_df_Conventional)
summary(aov_fungi_soil_Conventional_even)
HSD.test(aov_fungi_soil_Conventional_even, "Indicator_label") -> tukeyHSD_fungi_soil_Conventional_even
tukeyHSD_fungi_soil_Conventional_even
#soil No-Till even
aov_fungi_soil_No_Till_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_soil_fungi_df_No_Till)
summary(aov_fungi_soil_No_Till_even)
HSD.test(aov_fungi_soil_No_Till_even, "Indicator_label") -> tukeyHSD_fungi_soil_No_Till_even
tukeyHSD_fungi_soil_No_Till_even
# soil Organic even
alpha_div_soil_fungi_df_Organic
aov_fungi_soil_Organic_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_soil_fungi_df_Organic)
summary(aov_fungi_soil_Organic_even)
HSD.test(aov_fungi_soil_Organic_even, "Indicator_label") -> tukeyHSD_fungi_soil_Organic_even
tukeyHSD_fungi_soil_Organic_even


# make data frames split by growth stage
alpha_div_roots_fungi
alpha_div_roots_fungi
alpha_div_roots_fungi_df_V2 <- alpha_div_roots_fungi[alpha_div_roots_fungi$Growth_Stage == "V2",]
alpha_div_roots_fungi_df_V2
alpha_div_roots_fungi_df_R2 <- alpha_div_roots_fungi[alpha_div_roots_fungi$Growth_Stage == "R2",]
alpha_div_roots_fungi_df_R2 
alpha_div_roots_fungi_df_R6 <- alpha_div_roots_fungi[alpha_div_roots_fungi$Growth_Stage == "R6",]
alpha_div_roots_fungi_df_R6



#roots V2 rich
aov_fungi_roots_V2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_roots_fungi_df_V2)
summary(aov_fungi_roots_V2_rich)
HSD.test(aov_fungi_roots_V2_rich, "Indicator_label") -> tukeyHSD_fungi_roots_V2_rich
tukeyHSD_fungi_roots_V2_rich
#roots R2 rich
aov_fungi_roots_R2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_roots_fungi_df_R2)
summary(aov_fungi_roots_R2_rich)
HSD.test(aov_fungi_roots_R2_rich, "Indicator_label") -> tukeyHSD_fungi_roots_R2_rich
tukeyHSD_fungi_roots_R2_rich
# roots R6 rich
aov_fungi_roots_R6_rich <- aov(Observed ~ Indicator_label, data=alpha_div_roots_fungi_df_R6)
summary(aov_fungi_roots_R6_rich)
HSD.test(aov_fungi_roots_R6_rich, "Indicator_label") -> tukeyHSD_fungi_roots_R6_rich
tukeyHSD_fungi_roots_R6_rich


#roots V2 shan
aov_fungi_roots_V2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_roots_fungi_df_V2)
summary(aov_fungi_roots_V2_shan)
HSD.test(aov_fungi_roots_V2_shan, "Indicator_label") -> tukeyHSD_fungi_roots_V2_shan
tukeyHSD_fungi_roots_V2_shan
#roots R2 shan
aov_fungi_roots_R2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_roots_fungi_df_R2)
summary(aov_fungi_roots_R2_shan)
HSD.test(aov_fungi_roots_R2_shan, "Indicator_label") -> tukeyHSD_fungi_roots_R2_shan
tukeyHSD_fungi_roots_R2_shan
# roots R6 shan
aov_fungi_roots_R6_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_roots_fungi_df_R6)
summary(aov_fungi_roots_R6_shan)
HSD.test(aov_fungi_roots_R6_shan, "Indicator_label") -> tukeyHSD_fungi_roots_R6_shan
tukeyHSD_fungi_roots_R6_shan
alpha_div_roots_fungi_df_V2
#roots V2 even
aov_fungi_roots_V2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_roots_fungi_df_V2)
summary(aov_fungi_roots_V2_even)
HSD.test(aov_fungi_roots_V2_even, "Indicator_label") -> tukeyHSD_fungi_roots_V2_even
tukeyHSD_fungi_roots_V2_even
#roots R2 even
aov_fungi_roots_R2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_roots_fungi_df_R2)
summary(aov_fungi_roots_R2_even)
HSD.test(aov_fungi_roots_R2_even, "Indicator_label") -> tukeyHSD_fungi_roots_R2_even
tukeyHSD_fungi_roots_R2_even
# roots R6 even
aov_fungi_roots_R6_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_roots_fungi_df_R6)
summary(aov_fungi_roots_R6_even)
HSD.test(aov_fungi_roots_R6_even, "Indicator_label") -> tukeyHSD_fungi_roots_R6_even
tukeyHSD_fungi_roots_R6_even
#must also split into management systems 
alpha_div_roots_fungi_df_Conventional <- alpha_div_roots_fungi[alpha_div_roots_fungi$Management == "Conventional",]
alpha_div_roots_fungi_df_Conventional
alpha_div_roots_fungi_df_No_Till <- alpha_div_roots_fungi[alpha_div_roots_fungi$Management == "No-Till",]
alpha_div_roots_fungi_df_No_Till
alpha_div_roots_fungi_df_Organic <- alpha_div_roots_fungi[alpha_div_roots_fungi$Management == "Organic",]
alpha_div_roots_fungi_df_Organic

#roots Conventional rich
aov_fungi_roots_Conventional_rich <- aov(Observed ~ Indicator_label, data=alpha_div_roots_fungi_df_Conventional)
aov_fungi_roots_Conventional_rich
summary(aov_fungi_roots_Conventional_rich)
HSD.test(aov_fungi_roots_Conventional_rich, "Indicator_label") -> tukeyHSD_fungi_roots_Conventional_rich
tukeyHSD_fungi_roots_Conventional_rich
#roots No-Till rich
aov_fungi_roots_No_Till_rich <- aov(Observed ~ Indicator_label, data=alpha_div_roots_fungi_df_No_Till)
summary(aov_fungi_roots_No_Till_rich)
HSD.test(aov_fungi_roots_No_Till_rich, "Indicator_label") -> tukeyHSD_fungi_roots_No_Till_rich
tukeyHSD_fungi_roots_No_Till_rich
# roots Organic rich
alpha_div_roots_fungi_df_Organic
aov_fungi_roots_Organic_rich <- aov(Observed ~ Indicator_label, data=alpha_div_roots_fungi_df_Organic)
summary(aov_fungi_roots_Organic_rich)
HSD.test(aov_fungi_roots_Organic_rich, "Indicator_label") -> tukeyHSD_fungi_roots_Organic_rich
tukeyHSD_fungi_roots_Organic_rich


#roots Conventional shan
aov_fungi_roots_Conventional_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_roots_fungi_df_Conventional)
summary(aov_fungi_roots_Conventional_shan)
HSD.test(aov_fungi_roots_Conventional_shan, "Indicator_label") -> tukeyHSD_fungi_roots_Conventional_shan
tukeyHSD_fungi_roots_Conventional_shan
#roots No-Till shan
aov_fungi_roots_No_Till_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_roots_fungi_df_No_Till)
summary(aov_fungi_roots_No_Till_shan)
HSD.test(aov_fungi_roots_No_Till_shan, "Indicator_label") -> tukeyHSD_fungi_roots_No_Till_shan
tukeyHSD_fungi_roots_No_Till_shan
# roots Organic shan
alpha_div_roots_fungi_df_Organic
aov_fungi_roots_Organic_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_roots_fungi_df_Organic)
summary(aov_fungi_roots_Organic_shan)
HSD.test(aov_fungi_roots_Organic_shan, "Indicator_label") -> tukeyHSD_fungi_roots_Organic_shan
tukeyHSD_fungi_roots_Organic_shan

#roots Conventional even
aov_fungi_roots_Conventional_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_roots_fungi_df_Conventional)
summary(aov_fungi_roots_Conventional_even)
HSD.test(aov_fungi_roots_Conventional_even, "Indicator_label") -> tukeyHSD_fungi_roots_Conventional_even
tukeyHSD_fungi_roots_Conventional_even
#roots No-Till even
aov_fungi_roots_No_Till_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_roots_fungi_df_No_Till)
summary(aov_fungi_roots_No_Till_even)
HSD.test(aov_fungi_roots_No_Till_even, "Indicator_label") -> tukeyHSD_fungi_roots_No_Till_even
tukeyHSD_fungi_roots_No_Till_even
# roots Organic even
alpha_div_roots_fungi_df_Organic
aov_fungi_roots_Organic_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_roots_fungi_df_Organic)
summary(aov_fungi_roots_Organic_even)
HSD.test(aov_fungi_roots_Organic_even, "Indicator_label") -> tukeyHSD_fungi_roots_Organic_even
tukeyHSD_fungi_roots_Organic_even

# make data frames split by growth stage
alpha_div_stems_fungi
alpha_div_stems_fungi
alpha_div_stems_fungi_df_V2 <- alpha_div_stems_fungi[alpha_div_stems_fungi$Growth_Stage == "V2",]
alpha_div_stems_fungi_df_V2
alpha_div_stems_fungi_df_R2 <- alpha_div_stems_fungi[alpha_div_stems_fungi$Growth_Stage == "R2",]
alpha_div_stems_fungi_df_R2 
alpha_div_stems_fungi_df_R6 <- alpha_div_stems_fungi[alpha_div_stems_fungi$Growth_Stage == "R6",]
alpha_div_stems_fungi_df_R6



#stems V2 rich
aov_fungi_stems_V2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_stems_fungi_df_V2)
summary(aov_fungi_stems_V2_rich)
HSD.test(aov_fungi_stems_V2_rich, "Indicator_label") -> tukeyHSD_fungi_stems_V2_rich
tukeyHSD_fungi_stems_V2_rich
#stems R2 rich
aov_fungi_stems_R2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_stems_fungi_df_R2)
summary(aov_fungi_stems_R2_rich)
HSD.test(aov_fungi_stems_R2_rich, "Indicator_label") -> tukeyHSD_fungi_stems_R2_rich
tukeyHSD_fungi_stems_R2_rich
# stems R6 rich
aov_fungi_stems_R6_rich <- aov(Observed ~ Indicator_label, data=alpha_div_stems_fungi_df_R6)
summary(aov_fungi_stems_R6_rich)
HSD.test(aov_fungi_stems_R6_rich, "Indicator_label") -> tukeyHSD_fungi_stems_R6_rich
tukeyHSD_fungi_stems_R6_rich


#stems V2 shan
aov_fungi_stems_V2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_stems_fungi_df_V2)
summary(aov_fungi_stems_V2_shan)
HSD.test(aov_fungi_stems_V2_shan, "Indicator_label") -> tukeyHSD_fungi_stems_V2_shan
tukeyHSD_fungi_stems_V2_shan
#stems R2 shan
aov_fungi_stems_R2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_stems_fungi_df_R2)
summary(aov_fungi_stems_R2_shan)
HSD.test(aov_fungi_stems_R2_shan, "Indicator_label") -> tukeyHSD_fungi_stems_R2_shan
tukeyHSD_fungi_stems_R2_shan
# stems R6 shan
aov_fungi_stems_R6_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_stems_fungi_df_R6)
summary(aov_fungi_stems_R6_shan)
HSD.test(aov_fungi_stems_R6_shan, "Indicator_label") -> tukeyHSD_fungi_stems_R6_shan
tukeyHSD_fungi_stems_R6_shan
alpha_div_stems_fungi_df_V2
#stems V2 even
aov_fungi_stems_V2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_stems_fungi_df_V2)
summary(aov_fungi_stems_V2_even)
HSD.test(aov_fungi_stems_V2_even, "Indicator_label") -> tukeyHSD_fungi_stems_V2_even
tukeyHSD_fungi_stems_V2_even
#stems R2 even
aov_fungi_stems_R2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_stems_fungi_df_R2)
summary(aov_fungi_stems_R2_even)
HSD.test(aov_fungi_stems_R2_even, "Indicator_label") -> tukeyHSD_fungi_stems_R2_even
tukeyHSD_fungi_stems_R2_even
# stems R6 even
aov_fungi_stems_R6_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_stems_fungi_df_R6)
summary(aov_fungi_stems_R6_even)
HSD.test(aov_fungi_stems_R6_even, "Indicator_label") -> tukeyHSD_fungi_stems_R6_even
tukeyHSD_fungi_stems_R6_even
#must also split into management systems 
alpha_div_stems_fungi_df_Conventional <- alpha_div_stems_fungi[alpha_div_stems_fungi$Management == "Conventional",]
alpha_div_stems_fungi_df_Conventional
alpha_div_stems_fungi_df_No_Till <- alpha_div_stems_fungi[alpha_div_stems_fungi$Management == "No-Till",]
alpha_div_stems_fungi_df_No_Till
alpha_div_stems_fungi_df_Organic <- alpha_div_stems_fungi[alpha_div_stems_fungi$Management == "Organic",]
alpha_div_stems_fungi_df_Organic

#stems Conventional rich
aov_fungi_stems_Conventional_rich <- aov(Observed ~ Indicator_label, data=alpha_div_stems_fungi_df_Conventional)
aov_fungi_stems_Conventional_rich
summary(aov_fungi_stems_Conventional_rich)
HSD.test(aov_fungi_stems_Conventional_rich, "Indicator_label") -> tukeyHSD_fungi_stems_Conventional_rich
tukeyHSD_fungi_stems_Conventional_rich
#stems No-Till rich
aov_fungi_stems_No_Till_rich <- aov(Observed ~ Indicator_label, data=alpha_div_stems_fungi_df_No_Till)
summary(aov_fungi_stems_No_Till_rich)
HSD.test(aov_fungi_stems_No_Till_rich, "Indicator_label") -> tukeyHSD_fungi_stems_No_Till_rich
tukeyHSD_fungi_stems_No_Till_rich
# stems Organic rich
alpha_div_stems_fungi_df_Organic
aov_fungi_stems_Organic_rich <- aov(Observed ~ Indicator_label, data=alpha_div_stems_fungi_df_Organic)
summary(aov_fungi_stems_Organic_rich)
HSD.test(aov_fungi_stems_Organic_rich, "Indicator_label") -> tukeyHSD_fungi_stems_Organic_rich
tukeyHSD_fungi_stems_Organic_rich


#stems Conventional shan
aov_fungi_stems_Conventional_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_stems_fungi_df_Conventional)
summary(aov_fungi_stems_Conventional_shan)
HSD.test(aov_fungi_stems_Conventional_shan, "Indicator_label") -> tukeyHSD_fungi_stems_Conventional_shan
tukeyHSD_fungi_stems_Conventional_shan
#stems No-Till shan
aov_fungi_stems_No_Till_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_stems_fungi_df_No_Till)
summary(aov_fungi_stems_No_Till_shan)
HSD.test(aov_fungi_stems_No_Till_shan, "Indicator_label") -> tukeyHSD_fungi_stems_No_Till_shan
tukeyHSD_fungi_stems_No_Till_shan
# stems Organic shan
alpha_div_stems_fungi_df_Organic
aov_fungi_stems_Organic_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_stems_fungi_df_Organic)
summary(aov_fungi_stems_Organic_shan)
HSD.test(aov_fungi_stems_Organic_shan, "Indicator_label") -> tukeyHSD_fungi_stems_Organic_shan
tukeyHSD_fungi_stems_Organic_shan

#stems Conventional even
aov_fungi_stems_Conventional_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_stems_fungi_df_Conventional)
summary(aov_fungi_stems_Conventional_even)
HSD.test(aov_fungi_stems_Conventional_even, "Indicator_label") -> tukeyHSD_fungi_stems_Conventional_even
tukeyHSD_fungi_stems_Conventional_even
#stems No-Till even
aov_fungi_stems_No_Till_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_stems_fungi_df_No_Till)
summary(aov_fungi_stems_No_Till_even)
HSD.test(aov_fungi_stems_No_Till_even, "Indicator_label") -> tukeyHSD_fungi_stems_No_Till_even
tukeyHSD_fungi_stems_No_Till_even
# stems Organic even
alpha_div_stems_fungi_df_Organic
aov_fungi_stems_Organic_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_stems_fungi_df_Organic)
summary(aov_fungi_stems_Organic_even)
HSD.test(aov_fungi_stems_Organic_even, "Indicator_label") -> tukeyHSD_fungi_stems_Organic_even
tukeyHSD_fungi_stems_Organic_even

# make data frames split by growth stage
alpha_div_leaves_fungi
alpha_div_leaves_fungi
alpha_div_leaves_fungi_df_V2 <- alpha_div_leaves_fungi[alpha_div_leaves_fungi$Growth_Stage == "V2",]
alpha_div_leaves_fungi_df_V2
alpha_div_leaves_fungi_df_R2 <- alpha_div_leaves_fungi[alpha_div_leaves_fungi$Growth_Stage == "R2",]
alpha_div_leaves_fungi_df_R2 
alpha_div_leaves_fungi_df_R6 <- alpha_div_leaves_fungi[alpha_div_leaves_fungi$Growth_Stage == "R6",]
alpha_div_leaves_fungi_df_R6




#leaves V2 rich
aov_fungi_leaves_V2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_leaves_fungi_df_V2)
summary(aov_fungi_leaves_V2_rich)
HSD.test(aov_fungi_leaves_V2_rich, "Indicator_label") -> tukeyHSD_fungi_leaves_V2_rich
tukeyHSD_fungi_leaves_V2_rich
#leaves R2 rich
aov_fungi_leaves_R2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_leaves_fungi_df_R2)
summary(aov_fungi_leaves_R2_rich)
HSD.test(aov_fungi_leaves_R2_rich, "Indicator_label") -> tukeyHSD_fungi_leaves_R2_rich
tukeyHSD_fungi_leaves_R2_rich
# leaves R6 rich
aov_fungi_leaves_R6_rich <- aov(Observed ~ Indicator_label, data=alpha_div_leaves_fungi_df_R6)
summary(aov_fungi_leaves_R6_rich)
HSD.test(aov_fungi_leaves_R6_rich, "Indicator_label") -> tukeyHSD_fungi_leaves_R6_rich
tukeyHSD_fungi_leaves_R6_rich


#leaves V2 shan
aov_fungi_leaves_V2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_leaves_fungi_df_V2)
summary(aov_fungi_leaves_V2_shan)
HSD.test(aov_fungi_leaves_V2_shan, "Indicator_label") -> tukeyHSD_fungi_leaves_V2_shan
tukeyHSD_fungi_leaves_V2_shan
#leaves R2 shan
aov_fungi_leaves_R2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_leaves_fungi_df_R2)
summary(aov_fungi_leaves_R2_shan)
HSD.test(aov_fungi_leaves_R2_shan, "Indicator_label") -> tukeyHSD_fungi_leaves_R2_shan
tukeyHSD_fungi_leaves_R2_shan
# leaves R6 shan
aov_fungi_leaves_R6_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_leaves_fungi_df_R6)
summary(aov_fungi_leaves_R6_shan)
HSD.test(aov_fungi_leaves_R6_shan, "Indicator_label") -> tukeyHSD_fungi_leaves_R6_shan
tukeyHSD_fungi_leaves_R6_shan
alpha_div_leaves_fungi_df_V2
#leaves V2 even
aov_fungi_leaves_V2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_leaves_fungi_df_V2)
summary(aov_fungi_leaves_V2_even)
HSD.test(aov_fungi_leaves_V2_even, "Indicator_label") -> tukeyHSD_fungi_leaves_V2_even
tukeyHSD_fungi_leaves_V2_even
#leaves R2 even
aov_fungi_leaves_R2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_leaves_fungi_df_R2)
summary(aov_fungi_leaves_R2_even)
HSD.test(aov_fungi_leaves_R2_even, "Indicator_label") -> tukeyHSD_fungi_leaves_R2_even
tukeyHSD_fungi_leaves_R2_even
# leaves R6 even
aov_fungi_leaves_R6_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_leaves_fungi_df_R6)
summary(aov_fungi_leaves_R6_even)
HSD.test(aov_fungi_leaves_R6_even, "Indicator_label") -> tukeyHSD_fungi_leaves_R6_even
tukeyHSD_fungi_leaves_R6_even
#must also split into management systems 
alpha_div_leaves_fungi_df_Conventional <- alpha_div_leaves_fungi[alpha_div_leaves_fungi$Management == "Conventional",]
alpha_div_leaves_fungi_df_Conventional
alpha_div_leaves_fungi_df_No_Till <- alpha_div_leaves_fungi[alpha_div_leaves_fungi$Management == "No-Till",]
alpha_div_leaves_fungi_df_No_Till
alpha_div_leaves_fungi_df_Organic <- alpha_div_leaves_fungi[alpha_div_leaves_fungi$Management == "Organic",]
alpha_div_leaves_fungi_df_Organic

#leaves Conventional rich
aov_fungi_leaves_Conventional_rich <- aov(Observed ~ Indicator_label, data=alpha_div_leaves_fungi_df_Conventional)
aov_fungi_leaves_Conventional_rich
summary(aov_fungi_leaves_Conventional_rich)
HSD.test(aov_fungi_leaves_Conventional_rich, "Indicator_label") -> tukeyHSD_fungi_leaves_Conventional_rich
tukeyHSD_fungi_leaves_Conventional_rich
#leaves No-Till rich
aov_fungi_leaves_No_Till_rich <- aov(Observed ~ Indicator_label, data=alpha_div_leaves_fungi_df_No_Till)
summary(aov_fungi_leaves_No_Till_rich)
HSD.test(aov_fungi_leaves_No_Till_rich, "Indicator_label") -> tukeyHSD_fungi_leaves_No_Till_rich
tukeyHSD_fungi_leaves_No_Till_rich
# leaves Organic rich
alpha_div_leaves_fungi_df_Organic
aov_fungi_leaves_Organic_rich <- aov(Observed ~ Indicator_label, data=alpha_div_leaves_fungi_df_Organic)
summary(aov_fungi_leaves_Organic_rich)
HSD.test(aov_fungi_leaves_Organic_rich, "Indicator_label") -> tukeyHSD_fungi_leaves_Organic_rich
tukeyHSD_fungi_leaves_Organic_rich


#leaves Conventional shan
aov_fungi_leaves_Conventional_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_leaves_fungi_df_Conventional)
summary(aov_fungi_leaves_Conventional_shan)
HSD.test(aov_fungi_leaves_Conventional_shan, "Indicator_label") -> tukeyHSD_fungi_leaves_Conventional_shan
tukeyHSD_fungi_leaves_Conventional_shan
#leaves No-Till shan
aov_fungi_leaves_No_Till_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_leaves_fungi_df_No_Till)
summary(aov_fungi_leaves_No_Till_shan)
HSD.test(aov_fungi_leaves_No_Till_shan, "Indicator_label") -> tukeyHSD_fungi_leaves_No_Till_shan
tukeyHSD_fungi_leaves_No_Till_shan
# leaves Organic shan
alpha_div_leaves_fungi_df_Organic
aov_fungi_leaves_Organic_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_leaves_fungi_df_Organic)
summary(aov_fungi_leaves_Organic_shan)
HSD.test(aov_fungi_leaves_Organic_shan, "Indicator_label") -> tukeyHSD_fungi_leaves_Organic_shan
tukeyHSD_fungi_leaves_Organic_shan

#leaves Conventional even
aov_fungi_leaves_Conventional_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_leaves_fungi_df_Conventional)
summary(aov_fungi_leaves_Conventional_even)
HSD.test(aov_fungi_leaves_Conventional_even, "Indicator_label") -> tukeyHSD_fungi_leaves_Conventional_even
tukeyHSD_fungi_leaves_Conventional_even
#leaves No-Till even
aov_fungi_leaves_No_Till_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_leaves_fungi_df_No_Till)
summary(aov_fungi_leaves_No_Till_even)
HSD.test(aov_fungi_leaves_No_Till_even, "Indicator_label") -> tukeyHSD_fungi_leaves_No_Till_even
tukeyHSD_fungi_leaves_No_Till_even
# leaves Organic even
alpha_div_leaves_fungi_df_Organic
aov_fungi_leaves_Organic_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_leaves_fungi_df_Organic)
summary(aov_fungi_leaves_Organic_even)
HSD.test(aov_fungi_leaves_Organic_even, "Indicator_label") -> tukeyHSD_fungi_leaves_Organic_even
tukeyHSD_fungi_leaves_Organic_even

# note on the table - significance groups before the slash are within individual growth stage
# significance groups after the slash are within a single management


# filtering otus---------------------------------------------------------
# will filter now before creating barplots

# any sample with less than 5 reads for a particular otu will be placed to 0

otu_table(ps.noncontam_soil_obj1)[otu_table(ps.noncontam_soil_obj1) <= 4] <- 0 ### tag switching
ps.noncontam_soil_obj1
sample_data(ps.noncontam_soil_obj1)
otu_table(ps.noncontam_roots_obj1)[otu_table(ps.noncontam_roots_obj1) <= 4] <- 0 
otu_table(ps.noncontam_stems_obj1)[otu_table(ps.noncontam_stems_obj1) <= 4] <- 0 
otu_table(ps.noncontam_leaves_obj1)[otu_table(ps.noncontam_leaves_obj1) <= 4] <- 0 

# removes any OTUs that has less than 10 total reads across all samples
otu_table(ps.noncontam_soil_obj1) <- otu_table(ps.noncontam_soil_obj1)[which(rowSums(otu_table(ps.noncontam_soil_obj1)) >= 10),]### PCR Errors 
sample_data(ps.noncontam_soil_obj1)
otu_table(ps.noncontam_roots_obj1) <- otu_table(ps.noncontam_roots_obj1)[which(rowSums(otu_table(ps.noncontam_roots_obj1)) >= 10),]
otu_table(ps.noncontam_stems_obj1) <- otu_table(ps.noncontam_stems_obj1)[which(rowSums(otu_table(ps.noncontam_stems_obj1)) >= 10),]
otu_table(ps.noncontam_leaves_obj1) <- otu_table(ps.noncontam_leaves_obj1)[which(rowSums(otu_table(ps.noncontam_leaves_obj1)) >= 10),]
otu_table(ps.noncontam_leaves_obj1)

#check for additional samples soil 1000
write.csv(otu_table(ps.noncontam_stems_obj1), file = "stems_check.csv")
write.csv(otu_table(ps.noncontam_roots_obj1), file = "roots_check.csv")
write.csv(otu_table(ps.noncontam_soil_obj1), file = "soil_check.csv")
write.csv(otu_table(ps.noncontam_leaves_obj1), file = "leaves_check.csv")
# remove additional samples that went soil 1000 reads following this processesing
otu_table(ps.noncontam_roots_obj1) <- subset(otu_table(ps.noncontam_roots_obj1),
                                             select = -c(T4R2BV2R, T1R2AV2R))

# leaves
otu_table(ps.noncontam_leaves_obj1) <- subset(otu_table(ps.noncontam_leaves_obj1),
                                             select = -c(T1R5AR2L))


otu_table(ps.noncontam_stems_obj1) <- subset(otu_table(ps.noncontam_stems_obj1),
                                         select = -c(T4R5CV2ST,T4R5AV2ST,T4R6CV2ST,T4R1CV2ST,T2R1BV2ST))
write.csv(otu_table(ps.noncontam_soil_obj1), file = "depth_check_soil.csv")
write.csv(otu_table(ps.noncontam_roots_obj1), file = "depth_check_root.csv")
write.csv(otu_table(ps.noncontam_stems_obj1), file = "depth_check_stems.csv")
write.csv(otu_table(ps.noncontam_leaves_obj1), file = "depth_check_leaves.csv")
###barplots by fungal genus---------------------------
library(data.table)
library(dplyr)
library(ggplot2)


#soil
soil_obj1_barplots <- merge_samples(ps.noncontam_soil_obj1, "bar_label")
sample_data(soil_obj1_barplots)
sample_data(soil_obj1_barplots)$bar_label <- factor(sample_data(soil_obj1_barplots)$bar_label,
                                                 levels = c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))
write.csv(sample_data(soil_obj1_barplots), file ="sample_data_soil_bp.csv")

soil_obj1_barplots <- soil_obj1_barplots %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus
soil_obj1_barplots

soil_obj1_barplots
dat_soil_bp <- data.table(soil_obj1_barplots)
dat_soil_bp
dat_soil_bp[(Abundance <= 0.04), Genus:= "Other"]
dev.off()
# Plot
bar_ITS_soil= ggplot(dat_soil_bp, aes(x = Sample, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))+
  scale_fill_manual(values = c("Trichocladium" ="#652926",
                               "Sporidiobolus" = "steelblue", 
                               "Phoma" = "#C84248",
                               "Phaeosphaeriopsis" = "darksalmon",
                               "Phaeosphaeria" = "green",
                               "Ophiosphaerella" = "#CD9BCD",
                               "Myrothecium" = "#AD6F3B",
                               "Mycosphaerella" = "#673770", 
                               "Lewia" = "#D14285",
                               "Leptosphaerulina" = "#652926",
                               "Other" = "Blue",
                               "Hannaella" ="pink",
                               "Fusarium" ="#673770",
                               "Didymella" = "#AD6F3B",
                               "Diaporthe" ="#CBD588",
                               "Davidiella" = "#5F7FC7", 
                               "Cryptococcus" = "orange",
                               "Coniothyrium" = "#DA5724",
                               "Alternaria" = "#508578",
                               "Massilia" = "#CD9BCD",
                               "Pichia" = "tan",
                               "Dioszegia" = "magenta1",
                               "Edenia" = "gray52",
                               "Filobasidium" = "darkorange4",
                               "Malassezia" = "lightsalmon1",
                               "Microdochium" = "palevioletred3",
                               "Neosetophoma" = "olivedrab1",
                               "Paraphoma" = "cyan1",
                               "Gibellulopsis" = "darkviolet",
                               "Thielaviopsis" = "magenta4",
                               "Rhizophagus" = "greenyellow",
                               "Podospora" = "sienna1",
                               "Periconia" = "deepskyblue1",
                               "Neonectria" = "honeydew2",
                               "Nectria" = "red",
                               "Mycoleptodiscus" = "aquamarine4",
                               "Mortierella" = "gold",
                               "Metacordyceps" = "plum3",
                               "Macrophomina" = "peachpuff3",
                               "Lysurus" = "turquoise3",
                               "Leucoagaricus" = "blueviolet",
                               "Lachnum" = "green4",
                               "Knufia" ="palegreen2",
                               "Herpotrichia" = "hotpink3",
                               "Glomus" = "gray38",
                               "Funneliformis" = "black",
                               "Cylindrocarpon" = "midnightblue",
                               "Cyathus" = "khaki3",
                               "Crocicreas" = "yellow4",
                               "Corynespora" = "coral2",
                               "Ceratobasidium" = "rosybrown3",
                               "Bionectria" = "gray69",
                               "Tetracladium" = "paleturquoise1",
                               "Talaromyces" = "limegreen",
                               "Stropharia" = "lemonchiffon2",
                               "Sistotrema" = "deeppink1",
                               "Rhizophydium" = "plum4",
                               "Phialosimplex" = "darkolivegreen1",
                               "Panaeolus" = "darkorchid1",
                               "Myrmecridium" = "rosybrown2",
                               "Hypocrea" = "cornflowerblue",
                               "Humicola" = "lightgoldenrodyellow",
                               "Geotrichum" = "springgreen1",
                               "Emericellopsis" = "moccasin",
                               "Devriesia" = "lavenderblush",
                               "Coprinellus" = "aquamarine1",
                               "Conocybe" = "navyblue",
                               "Archaeorhizomyces" = "yellow",
                               "Penicillium" = "gray28"
                               
                               
                               
  ))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 0, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("soil")+
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  #theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  
  
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")
plot(bar_ITS_soil)





#roots
root_obj1_barplots <- merge_samples(ps.noncontam_roots_obj1, "bar_label")
sample_data(root_obj1_barplots)
sample_data(root_obj1_barplots)$Sample <- factor(sample_data(root_obj1_barplots)$bar_label,
                                                  levels = c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))
write.csv(sample_data(root_obj1_barplots), file ="sample_data_root_bp.csv")

root_obj1_barplots <- root_obj1_barplots %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus
root_obj1_barplots


dat_root_bp <- data.table(root_obj1_barplots)
dat_root_bp[(Abundance <= 0.04), Genus := "Other"]

# Plot
bar_ITS_root= ggplot(dat_root_bp, aes(x = Sample, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))+
  scale_fill_manual(values = c("Trichocladium" ="#652926",
                               "Sporidiobolus" = "steelblue", 
                               "Phoma" = "#C84248",
                               "Phaeosphaeriopsis" = "darksalmon",
                               "Phaeosphaeria" = "green",
                               "Ophiosphaerella" = "#CD9BCD",
                               "Myrothecium" = "#AD6F3B",
                               "Mycosphaerella" = "#673770", 
                               "Lewia" = "#D14285",
                               "Leptosphaerulina" = "#652926",
                               "Other" = "Blue",
                               "Hannaella" ="pink",
                               "Fusarium" ="#673770",
                               "Didymella" = "#AD6F3B",
                               "Diaporthe" ="#CBD588",
                               "Davidiella" = "#5F7FC7", 
                               "Cryptococcus" = "orange",
                               "Coniothyrium" = "#DA5724",
                               "Alternaria" = "#508578",
                               "Massilia" = "#CD9BCD",
                               "Pichia" = "tan",
                               "Dioszegia" = "magenta1",
                               "Edenia" = "gray52",
                               "Filobasidium" = "darkorange4",
                               "Malassezia" = "lightsalmon1",
                               "Microdochium" = "palevioletred3",
                               "Neosetophoma" = "olivedrab1",
                               "Paraphoma" = "cyan1",
                               "Gibellulopsis" = "darkviolet",
                               "Thielaviopsis" = "lavender",
                               "Rhizophagus" = "greenyellow",
                               "Podospora" = "sienna1",
                               "Periconia" = "deepskyblue1",
                               "Neonectria" = "honeydew2",
                               "Nectria" = "red",
                               "Mycoleptodiscus" = "aquamarine4",
                               "Mortierella" = "gold",
                               "Metacordyceps" = "plum3",
                               "Macrophomina" = "peachpuff3",
                               "Lysurus" = "turquoise3",
                               "Leucoagaricus" = "blueviolet",
                               "Lachnum" = "green4",
                               "Knufia" ="palegreen2",
                               "Herpotrichia" = "hotpink3",
                               "Glomus" = "gray38",
                               "Funneliformis" = "black",
                               "Cylindrocarpon" = "midnightblue",
                               "Cyathus" = "khaki3",
                               "Crocicreas" = "yellow4",
                               "Corynespora" = "cadetblue",
                               "Ceratobasidium" = "rosybrown3",
                               "Bionectria" = "gray69",
                               "Tetracladium" = "paleturquoise1",
                               "Talaromyces" = "limegreen",
                               "Stropharia" = "lemonchiffon2",
                               "Sistotrema" = "deeppink1",
                               "Rhizophydium" = "plum4",
                               "Phialosimplex" = "darkolivegreen1",
                               "Panaeolus" = "darkorchid1",
                               "Myrmecridium" = "rosybrown2",
                               "Hypocrea" = "cornflowerblue",
                               "Humicola" = "lightgoldenrodyellow",
                               "Geotrichum" = "springgreen1",
                               "Emericellopsis" = "moccasin",
                               "Devriesia" = "lavenderblush",
                               "Coprinellus" = "aquamarine1",
                               "Conocybe" = "navyblue",
                               "Archaeorhizomyces" = "yellow"
  ))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme_classic()+
  ggtitle("Roots")+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  #theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")
plot(bar_ITS_root)






#stems
stems_obj1_barplots <- merge_samples(ps.noncontam_stems_obj1, "bar_label")

sample_data(stems_obj1_barplots)$Indicator_label <- factor(sample_data(stems_obj1_barplots)$bar_label,
                                                  levels = c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))
write.csv(sample_data(stems_obj1_barplots), file ="sample_data_stems_bp.csv")

stems_obj1_barplots <- stems_obj1_barplots %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus
stems_obj1_barplots


dat_stems_bp <- data.table(stems_obj1_barplots)
dat_stems_bp
dat_stems_bp[(Abundance <= 0.04), Genus:= "Other"]

# Plot
bar_ITS_stems= ggplot(dat_stems_bp, aes(x = Sample, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))+
  scale_fill_manual(values = c("Trichocladium" ="#652926",
                               "Sporidiobolus" = "steelblue", 
                               "Phoma" = "#C84248",
                               "Phaeosphaeriopsis" = "yellow",
                               "Phaeosphaeria" = "green",
                               "Ophiosphaerella" = "#CD9BCD",
                               "Myrothecium" = "#AD6F3B",
                               "Mycosphaerella" = "#673770", 
                               "Lewia" = "#D14285",
                               "Leptosphaerulina" = "#652926",
                               "Other" = "Blue",
                               "Hannaella" ="pink",
                               "Gibberella" = "#673770",
                               "Fusarium" ="#673770",
                               "Didymella" = "#AD6F3B",
                               "Diaporthe" ="#CBD588",
                               "Davidiella" = "orangered2", 
                               "Cryptococcus" = "orange",
                               "Coniothyrium" = "#DA5724",
                               "Alternaria" = "#508578",
                               "Massilia" = "#CD9BCD",
                               "Pichia" = "tan",
                               "Dioszegia" = "magenta1",
                               "Edenia" = "gray52",
                               "Filobasidium" = "darkorange4",
                               "Malassezia" = "lightsalmon1",
                               "Microdochium" = "palevioletred3",
                               "Neosetophoma" = "olivedrab1",
                               "Paraphoma" = "cyan1",
                               "Gibellulopsis" = "darkviolet"
                            
                               
  ))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("Stems")+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  #theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")
plot(bar_ITS_stems)


#leaves
leaves_obj1_barplots <- merge_samples(ps.noncontam_leaves_obj1, "bar_label")
sample_data(leaves_obj1_barplots)$bar_label <- factor(sample_data(leaves_obj1_barplots)$bar_label,
                                                   levels = c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))
write.csv(sample_data(leaves_obj1_barplots), file ="sample_data_leaves_bp.csv")

leaves_obj1_barplots <- leaves_obj1_barplots %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at Genus level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by Genus
leaves_obj1_barplots


dat_leaves_bp <- data.table(leaves_obj1_barplots)
dat_leaves_bp[(Abundance <= 0.04), Genus := "Other"]

# Plot
bar_ITS_leaves= ggplot(dat_leaves_bp, aes(x = Sample, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))+
  scale_fill_manual(values = c("Trichocladium" ="#652926",
                               "Sporidiobolus" = "burlywood2", 
                               "Phoma" = "#C84248",
                               "Phaeosphaeriopsis" = "yellow",
                               "Phaeosphaeria" = "green",
                               "Ophiosphaerella" = "#CD9BCD",
                               "Myrothecium" = "#AD6F3B",
                               "Mycosphaerella" = "snow4", 
                               "Lewia" = "#D14285",
                               "Leptosphaerulina" = "#652926",
                               "Other" = "Blue",
                               "Hannaella" ="pink",
                               "Gibberella" = "#673770",
                               "Fusarium" ="#673770",
                               "Didymella" = "#AD6F3B",
                              "Diaporthe" ="#CBD588",
                               "Davidiella" = "orangered2", 
                               "Cryptococcus" = "orange",
                               "Coniothyrium" = "#DA5724",
                               "Alternaria" = "#508578",
                               "Massilia" = "#CD9BCD",
                              "Neosetophoma" = "royalblue4",
                              "Dioszegia" = "deeppink",
                              "Bullera" = "darkorange3"
                               
                                
                        ))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme_classic()+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  ggtitle("Leaves")+
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  #theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance (Genera > 4%) \n") +
  xlab("")
plot(bar_ITS_leaves)



#combine barplots 
library(ggpubr)

# Figure 1-----------------------
ggarrange(bar_ITS_roots,bar_ITS_root,bar_ITS_stems,bar_ITS_leaves,
          labels = c("A", "B", "C", "D" ),
          widths = c(2.0, 2.0, 2.0, 2.0),
          align = "h", ncol = 4, nrow = 1)
ps.noncontam_leaves_obj1
ps.noncontam_roots_obj1
ps.noncontam_stems_obj1
ps.noncontam_leaves_obj1
library(tidyr)
# assesing community composition, retrieving abundances-------------------
# check abundances both by growth stage and management
# to get all of the percentages presented in the manuscript, just changed management and taxonomic levels presented below
ps.noncontam_soil_obj1_Phylum <- subset_samples(ps.noncontam_soil_obj1,Management%in%c("Organic"))
sample_data(ps.noncontam_soil_obj1_Phylum)
ps.noncontam_soil_obj1_Organic_Phylum= tax_glom(ps.noncontam_soil_obj1_Phylum, "Phylum")
otu_table(ps.noncontam_soil_obj1_Organic_Phylum)
tax_table <- tax_table(ps.noncontam_soil_obj1_Organic_Phylum)
tax_table
ps.noncontam_soil_obj1_Organic_Phylum= taxa_sums(ps.noncontam_soil_obj1_Organic_Phylum)/sum(taxa_sums(ps.noncontam_soil_obj1_Organic_Phylum))*100
ps.noncontam_soil_obj1_Organic_Phylum
ps.noncontam_soil_obj1_Organic_Phylum<- as.data.frame(ps.noncontam_soil_obj1_Organic_Phylum)
dim(ps.noncontam_soil_obj1_Organic_Phylum)
ps.noncontam_soil_obj1_Organic_Phylum<- setNames(ps.noncontam_soil_obj1_Organic_Phylum, c("abundance"))
# writing into excel to find to most abundant genera, will check with tax table above
write.csv(ps.noncontam_soil_obj1_Organic_Phylum, file = "Organic_Phylum.csv")

#remerge phyloseqs for beta diversity analysis and indicator species analysis

ps.noncontam_total_obj1= merge_phyloseq(ps.noncontam_leaves_obj1, ps.noncontam_soil_obj1, ps.noncontam_roots_obj1,ps.noncontam_stems_obj1)
sample_data(ps.noncontam_total_obj1)
ps.noncontam_above_obj1 = merge_phyloseq(ps.noncontam_leaves_obj1, ps.noncontam_stems_obj1)
ps.noncontam_below_obj1 = merge_phyloseq(ps.noncontam_roots_obj1, ps.noncontam_soil_obj1)

# indicator species analysis 
library("indicspecies")
library("ComplexHeatmap")
library("circlize")

#do above and below to make two heatmaps per fungal/prok dataset

#above 
#merging samples to combine replicate plots within a single Management at a singel growth stage for indicator species analysis
sample_data(ps.noncontam_above_obj1)
ps.noncontam_above_obj1_merged <- merge_samples(ps.noncontam_above_obj1, "Indicator_label")
sample_data(ps.noncontam_above_obj1_merged)


# make phyloseq components into dataframes
otu_above_obj1 <- as.data.frame(otu_table(ps.noncontam_above_obj1))
otu_above_obj1
tax_above_obj1 <- as.data.frame(as.matrix(tax_table(ps.noncontam_above_obj1)))
metadata_above_obj1 <- as.data.frame(as.matrix(sample_data(ps.noncontam_above_obj1)))
metadata_above_obj1
# perform indicator species analysis
isa_above_fungi <- multipatt(as.data.frame(t(otu_above_obj1)), metadata_above_obj1$Management_Indicator, control=how(nperm=9999))
summary(isa_above_fungi, indvalcomp=TRUE)
isa_above_fungi -> isa_above_fungi_Management_fdr
# perform indicator species analaysis filtering (fdr adjustment)
isa_above_fungi_Management_fdr$sign$p.value<-p.adjust(isa_above_fungi_Management_fdr$sign$p.value, "fdr")
isa_above_fungi_Management_fdr
summary(isa_above_fungi_Management_fdr)



sink(file="isa_above_fungi_Management_Origin.csv") 
summary(isa_above_fungi_Management_fdr)
sink()
isa_above_fungi_Management_fdr
# extracting ISA OTUs ----------------------------------------------------------------------------
# only keep otus with indicator p value below .05
result_isa_above_fungi_Management_fdr <- isa_above_fungi_Management_fdr$sign[which(isa_above_fungi_Management_fdr$sign$p.value <= 0.05), ]
result_isa_above_fungi_Management_fdr
dim(result_isa_above_fungi_Management_fdr)


result_isa_above_fungi_Management_fdr[result_isa_above_fungi_Management_fdr$s.No_Till==1 &
                                          result_isa_above_fungi_Management_fdr$s.Organic==0 &
                                          result_isa_above_fungi_Management_fdr$s.Conventional==0 ,] -> No_Till

result_isa_above_fungi_Management_fdr[result_isa_above_fungi_Management_fdr$s.No_Till==0 &
                                          result_isa_above_fungi_Management_fdr$s.Organic==1 &
                                          result_isa_above_fungi_Management_fdr$s.Conventional==0 ,] -> Organic

result_isa_above_fungi_Management_fdr[result_isa_above_fungi_Management_fdr$s.No_Till==0 &
                                          result_isa_above_fungi_Management_fdr$s.Organic==0 &
                                          result_isa_above_fungi_Management_fdr$s.Conventional==1,] -> Conventional

result_isa_above_fungi_Management_fdr[result_isa_above_fungi_Management_fdr$s.No_Till==1 &
                                          result_isa_above_fungi_Management_fdr$s.Organic==0 &
                                          result_isa_above_fungi_Management_fdr$s.Conventional==1 ,] -> No_Till_Conventional

result_isa_above_fungi_Management_fdr[result_isa_above_fungi_Management_fdr$s.No_Till==1&
                                          result_isa_above_fungi_Management_fdr$s.Organic==1 &
                                          result_isa_above_fungi_Management_fdr$s.Conventional==0 ,] -> No_Till_Organic

result_isa_above_fungi_Management_fdr[result_isa_above_fungi_Management_fdr$s.No_Till==0 &
                                          result_isa_above_fungi_Management_fdr$s.Organic==1 &
                                          result_isa_above_fungi_Management_fdr$s.Conventional==1 ,] -> Organic_Conventional

result_isa_above_fungi_Management_fdr
isa_above_Management_df <- rbind(No_Till,Organic,Conventional,No_Till_Conventional,No_Till_Organic,Organic_Conventional)
dim(isa_above_Management_df)
isa_above_Management_df

# phyloseq objets of ISA OTUs --------------------------------------------------------------------
ps.noncontam_above_obj1_merged -> ps_above_Management.isa
ps_above_Management.isa
ps_above_Management.isa = transform_sample_counts(ps_above_Management.isa, function(x) 100*x/sum(x)) # transform to relative abundances 
otu_table(ps_above_Management.isa) = otu_table(t(ps_above_Management.isa))
otu_table(ps_above_Management.isa)
otu_table(ps_above_Management.isa) <-otu_table(ps_above_Management.isa)[rownames(isa_above_Management_df),]
ps_above_Management.isa
sample_data(ps_above_Management.isa)
ps_above_Management.isa

install.packages("reltools")




isa_above_Management_otus <- as.data.frame(otu_table(ps_above_Management.isa))
dim(isa_above_Management_otus)

sample_data(ps_above_Management.isa)


# Creating a data.frame for plotting the heatmap 
identical(colnames(isa_above_Management_otus), rownames(sample_data(ps_above_Management.isa)))
isa_above_Management_otus
sample_data(ps_above_Management.isa)
colnames(isa_above_Management_otus) <- sample_data(ps_above_Management.isa)$Isa
colnames(isa_above_Management_df) <- c("Conventional", "No-Till", "Organic", "Index", "Stat", "p.value")
identical(rownames(isa_above_Management_df), rownames(isa_above_Management_otus))

isa_above_Management_obj <- cbind(isa_above_Management_otus, isa_above_Management_df)
isa_above_Management_obj$readNo <- rowSums(otu_above_obj1[rownames(isa_above_Management_df),])
isa_above_Management_obj$relAb <- (isa_above_Management_obj$readNo/sum(colSums(otu_above_obj1))) * 100
isa_above_Management_obj
isa_above_Management_obj$logAb <- log(isa_above_Management_obj$readNo)
isa_above_Management_obj$sqrtAb <- sqrt(isa_above_Management_obj$readNo)
isa_above_Management_obj <- isa_above_Management_obj[order(isa_above_Management_obj$relAb,decreasing = TRUE),]
isa_above_Management_obj
isa_above_Management_obj <- isa_above_Management_obj[1:30,]
dim(isa_above_Management_obj)
write.csv(sample_data(ps_above_Management.isa), file = "ps_above_Management.csv" )

# writing csv file and manually merging with taxonomy that was updated above in excel
# the taxonomy for the top 30 relatively abundant indicator species was checked using ncbi blast
# additionally, sample names were corrected in excel by consulting the ps_above_management.csv file
# and taxonomy was added to otu names
# after it was updated, the file was read back in
write.csv(isa_above_Management_obj, file = "isa_above_add_taxonomy.csv")
#isa_above_Management_obj <- read.csv("isa_above_add_taxonomy.csv", header=T, row.names = 1)
identical(rownames(tax_table(ps_above_Management.isa)), rownames(isa_above_Management_otus))

isa_above_Management_obj
sample_data(ps_above_Management.isa)

ht1_above_Management = Heatmap(as.matrix(sqrt(isa_above_Management_obj[,1:18]*10)), col = colorRamp2(c(0, 5), c("white","red")), 
                              cluster_rows = FALSE, cluster_columns = TRUE, name = "Abundance",
                              row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                              show_heatmap_legend = FALSE)
ht1_above_Management


ha_bar_above = HeatmapAnnotation("Rel. Abundance" = row_anno_barplot(isa_above_Management_obj$relAb, axis = FALSE, width = unit(.2, "%")), 
                                  which = "row", annotation_width = unit(1.75, "cm"), show_annotation_name = TRUE, annotation_name_gp = gpar(fontsize = 8), annotation_name_offset = unit(.75,"cm"),
                                  annotation_name_rot = c(0)) # annotation_width = unit(20, "cm"), annotation_height = unit(10, "cm")


ha_bar_above

ha_bar_above + ht1_above_Management



#below
#merging samples to combine replicate plots within a single Management at a singel growth stage for indicator species analysis
sample_data(ps.noncontam_below_obj1)
ps.noncontam_below_obj1_merged <- merge_samples(ps.noncontam_below_obj1, "Indicator_label")
sample_data(ps.noncontam_below_obj1_merged)


# make phyloseq components into dataframes
otu_below_obj1 <- as.data.frame(otu_table(ps.noncontam_below_obj1))
otu_below_obj1
tax_below_obj1 <- as.data.frame(as.matrix(tax_table(ps.noncontam_below_obj1)))
metadata_below_obj1 <- as.data.frame(as.matrix(sample_data(ps.noncontam_below_obj1)))
write.csv(metadata_below_obj1, file = "metadata_below_obj1.csv")
metadata_below_obj1 <- read.csv("metadata_below_obj1.csv")
write.csv(otu_below_obj1, file = "otu_below_obj1.csv")
# perform indicator species analysis
isa_below_fungi <- multipatt(as.data.frame(t(otu_below_obj1)), metadata_below_obj1$Management_Indicator, control=how(nperm=9999))
summary(isa_below_fungi, indvalcomp=TRUE)
# perform fdr transformation
isa_below_fungi -> isa_below_fungi_Management_fdr
isa_below_fungi_Management_fdr$sign$p.value<-p.adjust(isa_below_fungi_Management_fdr$sign$p.value, "fdr")
isa_below_fungi_Management_fdr
summary(isa_below_fungi_Management_fdr)



sink(file="isa_below_fungi_Management_Origin.csv") 
summary(isa_below_fungi_Management_fdr)
sink()
isa_below_fungi_Management_fdr
# extracting ISA OTUs ----------------------------------------------------------------------------
result_isa_below_fungi_Management_fdr <- isa_below_fungi_Management_fdr$sign[which(isa_below_fungi_Management_fdr$sign$p.value <= 0.05), ]
result_isa_below_fungi_Management_fdr
dim(result_isa_below_fungi_Management_fdr)


result_isa_below_fungi_Management_fdr[result_isa_below_fungi_Management_fdr$s.No_Till==1 &
                                        result_isa_below_fungi_Management_fdr$s.Organic==0 &
                                        result_isa_below_fungi_Management_fdr$s.Conventional==0 ,] -> No_Till

result_isa_below_fungi_Management_fdr[result_isa_below_fungi_Management_fdr$s.No_Till==0 &
                                        result_isa_below_fungi_Management_fdr$s.Organic==1 &
                                        result_isa_below_fungi_Management_fdr$s.Conventional==0 ,] -> Organic

result_isa_below_fungi_Management_fdr[result_isa_below_fungi_Management_fdr$s.No_Till==0 &
                                        result_isa_below_fungi_Management_fdr$s.Organic==0 &
                                        result_isa_below_fungi_Management_fdr$s.Conventional==1,] -> Conventional

result_isa_below_fungi_Management_fdr[result_isa_below_fungi_Management_fdr$s.No_Till==1 &
                                        result_isa_below_fungi_Management_fdr$s.Organic==0 &
                                        result_isa_below_fungi_Management_fdr$s.Conventional==1 ,] -> No_Till_Conventional

result_isa_below_fungi_Management_fdr[result_isa_below_fungi_Management_fdr$s.No_Till==1&
                                        result_isa_below_fungi_Management_fdr$s.Organic==1 &
                                        result_isa_below_fungi_Management_fdr$s.Conventional==0 ,] -> No_Till_Organic

result_isa_below_fungi_Management_fdr[result_isa_below_fungi_Management_fdr$s.No_Till==0 &
                                        result_isa_below_fungi_Management_fdr$s.Organic==1 &
                                        result_isa_below_fungi_Management_fdr$s.Conventional==1 ,] -> Organic_Conventional

result_isa_below_fungi_Management_fdr
isa_below_Management_df <- rbind(No_Till,Organic,Conventional,No_Till_Conventional,No_Till_Organic,Organic_Conventional)
dim(isa_below_Management_df)
isa_below_Management_df

# phyloseq objets of ISA OTUs --------------------------------------------------------------------
ps.noncontam_below_obj1_merged -> ps_below_Management.isa
ps_below_Management.isa
ps_below_Management.isa = transform_sample_counts(ps_below_Management.isa, function(x) 100*x/sum(x)) # transform to relativa abundances 
otu_table(ps_below_Management.isa) = otu_table(t(ps_below_Management.isa))
otu_table(ps_below_Management.isa)
otu_table(ps_below_Management.isa) <-otu_table(ps_below_Management.isa)[rownames(isa_below_Management_df),]
ps_below_Management.isa
sample_data(ps_below_Management.isa)
ps_below_Management.isa


### will determine if taxonomy needs to be improved later
isa_below_Management_otus <- as.data.frame(otu_table(ps_below_Management.isa))
dim(isa_below_Management_otus)

sample_data(ps_below_Management.isa)


# Creating a data.frame for plotting the heatmap 
identical(colnames(isa_below_Management_otus), rownames(sample_data(ps_below_Management.isa)))
isa_below_Management_otus
sample_data(ps_below_Management.isa)
colnames(isa_below_Management_otus) <- sample_data(ps_below_Management.isa)$Isa
colnames(isa_below_Management_df) <- c("Conventional", "No-Till", "Organic", "Index", "Stat", "p.value")
identical(rownames(isa_below_Management_df), rownames(isa_below_Management_otus))

isa_below_Management_obj <- cbind(isa_below_Management_otus, isa_below_Management_df)
isa_below_Management_obj$readNo <- rowSums(otu_below_obj1[rownames(isa_below_Management_df),])
isa_below_Management_obj$relAb <- (isa_below_Management_obj$readNo/sum(colSums(otu_below_obj1))) * 100
isa_below_Management_obj
isa_below_Management_obj$logAb <- log(isa_below_Management_obj$readNo)
isa_below_Management_obj$sqrtAb <- sqrt(isa_below_Management_obj$readNo)
isa_below_Management_obj <- isa_below_Management_obj[order(isa_below_Management_obj$relAb,decreasing = TRUE),]
isa_below_Management_obj
isa_below_Management_obj <- isa_below_Management_obj[1:30,]
dim(isa_below_Management_obj)
write.csv(sample_data(ps_below_Management.isa), file = "ps_below_Management.csv" )

# writing csv file and manually merging with taxonomy that was updated above in excel
# the taxonomy for the top 30 relatively abundant indicator species was checked using ncbi blast
# additionally, sample names were corrected in excel by consulting the ps_above_management.csv file
# and taxonomy was added to otu names
# after it was updated, the file was read back in
write.csv(isa_below_Management_obj, file = "isa_below_add_taxonomy.csv")
isa_below_Management_obj <- read.csv("isa_below_add_taxonomy.csv", header=T, row.names = 1)
identical(rownames(tax_table(ps_below_Management.isa)), rownames(isa_below_Management_otus))

isa_below_Management_obj
sample_data(ps_below_Management.isa)

ht1_below_Management = Heatmap(as.matrix(sqrt(isa_below_Management_obj[,1:18]*10)), col = colorRamp2(c(0, 5), c("white","red")), 
                               cluster_rows = FALSE, cluster_columns = TRUE, name = "Abundance",
                               row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 8),
                               show_heatmap_legend = FALSE)
ht1_below_Management


ha_bar_below = HeatmapAnnotation("Rel. Abundance" = row_anno_barplot(isa_below_Management_obj$relAb, axis = FALSE, width = unit(.2, "%")), 
                                 which = "row", annotation_width = unit(1.75, "cm"), show_annotation_name = TRUE, annotation_name_gp = gpar(fontsize = 8), annotation_name_offset = unit(.75,"cm"),
                                 annotation_name_rot = c(0)) # annotation_width = unit(20, "cm"), annotation_height = unit(10, "cm")


ha_bar_below

ha_bar_below + ht1_below_Management



# Heatmaps were arranged, along with prokaryotic indicator heatmaps outside of R


### normalizing with metagenome seq------------------------------------------------
source("http://bioconductor.org/biocLite.R")
biocLite("metagenomeSeq")
library("metagenomeSeq")

#total
# fitting into a Gaussian Model using metagenomeSeq-------------
# Normalizing before beta diversity analyses
write.csv(otu_table(ps.noncontam_total_obj1), file = "empty_sample_check.csv")
# overall

ps.noncontam_norm = phyloseq_to_metagenomeSeq(ps.noncontam_total_obj1)
p_biom<-cumNormStat(ps.noncontam_norm)
biom_quant<-cumNorm(ps.noncontam_norm, p=p_biom)
biom_quant
normFactors(biom_quant)
ps.noncontam_norm <-MRcounts(biom_quant, norm=T)
head(ps.noncontam_norm )
#create physeq object with normalized otu table
otu_table(ps.noncontam_total_obj1) <- otu_table(ps.noncontam_norm, taxa_are_rows = TRUE)

#above
# fitting into a Gaussian Model using metagenomeSeq-------------
write.csv(otu_table(ps.noncontam_above_obj1), file = "empty_sample_check.csv")
ps.noncontam_above_obj1_norm = phyloseq_to_metagenomeSeq(ps.noncontam_above_obj1)
p_biom_above<-cumNormStat(ps.noncontam_above_obj1_norm)
biom_quant_above<-cumNorm(ps.noncontam_above_obj1_norm, p=p_biom_above)
biom_quant_above
normFactors(biom_quant_above)
ps.noncontam_above_obj1_norm <-MRcounts(biom_quant_above, norm=T)
head(ps.noncontam_above_obj1_norm )
#create physeq object with normalized otu table
otu_table(ps.noncontam_above_obj1) <- otu_table(ps.noncontam_above_obj1_norm, taxa_are_rows = TRUE)

#below
# fitting into a Gaussian Model using metagenomeSeq-------------
write.csv(otu_table(ps.noncontam_below_obj1), file = "empty_sample_check.csv")
ps.noncontam_below_obj1_norm = phyloseq_to_metagenomeSeq(ps.noncontam_below_obj1)
p_biom_below<-cumNormStat(ps.noncontam_below_obj1_norm)
biom_quant_below<-cumNorm(ps.noncontam_below_obj1_norm, p=p_biom_below)
biom_quant_below
normFactors(biom_quant_below)
ps.noncontam_below_obj1_norm <-MRcounts(biom_quant_below, norm=T)
head(ps.noncontam_below_obj1_norm )
#create physeq object with normalized otu table
otu_table(ps.noncontam_below_obj1) <- otu_table(ps.noncontam_below_obj1_norm, taxa_are_rows = TRUE)








#leaves
ps.noncontam_leaves_norm = phyloseq_to_metagenomeSeq(ps.noncontam_leaves_obj1)
p_biom_leaves<-cumNormStat(ps.noncontam_leaves_norm)
biom_quant_leaves<-cumNorm(ps.noncontam_leaves_norm, p=p_biom_leaves)
biom_quant_leaves
normFactors(biom_quant_leaves)
ps.noncontam_leaves_norm <-MRcounts(biom_quant_leaves, norm=T)
head(ps.noncontam_leaves_norm )
#create physeq object with normalized otu table
otu_table(ps.noncontam_leaves_obj1) <- otu_table(ps.noncontam_leaves_norm, taxa_are_rows = TRUE)

#roots
ps.noncontam_roots_norm = phyloseq_to_metagenomeSeq(ps.noncontam_roots_obj1)
p_biom_roots<-cumNormStat(ps.noncontam_roots_norm)
biom_quant_roots<-cumNorm(ps.noncontam_roots_norm, p=p_biom_roots)
biom_quant_roots
normFactors(biom_quant_roots)
ps.noncontam_roots_norm <-MRcounts(biom_quant_roots, norm=T)
head(ps.noncontam_roots_norm )
#create physeq object with normalized otu table
otu_table(ps.noncontam_roots) <- otu_table(ps.noncontam_roots_norm, taxa_are_rows = TRUE)


#stems
ps.noncontam_stems_norm = phyloseq_to_metagenomeSeq(ps.noncontam_stems_obj1)
p_biom_stems<-cumNormStat(ps.noncontam_stems_norm)
biom_quant_stems<-cumNorm(ps.noncontam_stems_norm, p=p_biom_stems)
biom_quant_stems
normFactors(biom_quant_stems)
ps.noncontam_stems_norm <-MRcounts(biom_quant_stems, norm=T)
head(ps.noncontam_stems_norm )
#create physeq object with normalized otu table
otu_table(ps.noncontam_stems) <- otu_table(ps.noncontam_stems_norm, taxa_are_rows = TRUE)

#soil
ps.noncontam_soil_norm = phyloseq_to_metagenomeSeq(ps.noncontam_soil_obj1)
p_biom_soil<-cumNormStat(ps.noncontam_soil_norm)
biom_quant_soil<-cumNorm(ps.noncontam_soil_norm, p=p_biom_soil)
biom_quant_soil
normFactors(biom_quant_soil)
ps.noncontam_soil_norm <-MRcounts(biom_quant_soil, norm=T)
head(ps.noncontam_soil_norm )
#create physeq object with normalized otu table
otu_table(ps.noncontam_soil) <- otu_table(ps.noncontam_soil_norm, taxa_are_rows = TRUE)






# NMDS ------------------------------------
# Figure 3, fungal ordinations combined with prokaryotic ordinations outside of R
library("ggrepel")






# overall fungi
otu_table(ps.noncontam_total_obj1)
ord_ITS_obj_1 = ordinate(ps.noncontam_total_obj1, method ="NMDS", distance="bray", try=200)
ord_ITS_obj_1




# shape is stage
NMDS_ITS_Management = plot_ordination(ps.noncontam_total_obj1, ord_ITS_obj_1, color="origin", shape ="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  scale_colour_manual("origin",breaks = c("soil", "root","leaves", "stem"),
                      values = c("soil"="red", "root"="blue", "leaves" = "green", "stem" = "purple"))+
  #geom_text(aes(label=Description), size = 6)  +
  theme(legend.position="right")
plot(NMDS_ITS_Management)





#below, management

ord_ITS_obj_1_below = ordinate(ps.noncontam_below_obj1, method ="NMDS", distance="bray", try=200)
ord_ITS_obj_1_below


NMDS_ITS_below = plot_ordination(ps.noncontam_below_obj1, ord_ITS_obj_1_below, shape="Management", color ="origin") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  scale_colour_manual("origin",breaks = c("soil", "root"),
                      values = c("soil"="red", "root"="blue"))+
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_ITS_below)

#above, Management

ord_ITS_obj_1_above = ordinate(ps.noncontam_above_obj1, method ="NMDS", distance="bray", try=2000)
ord_ITS_obj_1_above


NMDS_ITS_above = plot_ordination(ps.noncontam_above_obj1, ord_ITS_obj_1_above, shape="Management", color ="origin") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  scale_colour_manual("origin",breaks = c("leaves", "stem"),
                      values = c( "leaves" = "green", "stem" = "purple"))+
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position = "none")
plot(NMDS_ITS_above)






library(ggpubr)

# to make supplemental figure 1, switch management as shape for Growth_Stage
ggarrange(NMDS_ITS_Management,NMDS_ITS_above, NMDS_ITS_below,labels = "AUTO", ncol =3, nrow=1,widths = c(2.5, 2.0,2.0))



###Split NMDS Supplemental Figure 2------------------------------------
#soil V2
otu_table(ps.noncontam_soil)
ps.noncontam_soil_V2 <- subset_samples(ps.noncontam_soil, Growth_Stage%in%c("V2"))
ord_soil_v2 = ordinate(ps.noncontam_soil_V2 , method ="NMDS", distance="bray", try=200)
ord_soil_v2 



# shape is management
NMDS_soil_V2= plot_ordination(ps.noncontam_soil_V2, ord_soil_v2 , color="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="left")
plot(NMDS_soil_V2)

# get R2 value for plot
otu_fungi_soil_V2<- as.data.frame(otu_table(ps.noncontam_soil_V2))
taxa_fungi_soil_V2<- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_V2)))
metadata_fungi_soil_V2 <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_V2)))
metadata_fungi_soil_V2
library("vegan")
library("RVAideMemoire")

adonis(t(otu_fungi_soil_V2) ~ Management, data=metadata_fungi_soil_V2, permutations=9999)
vegan::vegdist(t(otu_fungi_soil_V2), method="bray") -> dist_otu_fungi_soil_V2

permdisp_otu_fungi_soil_V2<- betadisper(dist_otu_fungi_soil_V2, metadata_fungi_soil_V2$Management)
anova(permdisp_otu_fungi_soil_V2, permutations = 9999)
#soil R2
ps.noncontam_soil_R2 <- subset_samples(ps.noncontam_soil, Growth_Stage%in%c("R2"))
ord_soil_R2 = ordinate(ps.noncontam_soil_R2 , method ="NMDS", distance="bray", try=200)
ord_soil_R2 



# shape is management
NMDS_soil_R2= plot_ordination(ps.noncontam_soil_R2, ord_soil_R2 , color="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_soil_R2)

# get R2 value for plot
otu_fungi_soil_R2<- as.data.frame(otu_table(ps.noncontam_soil_R2))
taxa_fungi_soil_R2<- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_R2)))
metadata_fungi_soil_R2 <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_R2)))
metadata_fungi_soil_R2


adonis(t(otu_fungi_soil_R2) ~ Management, data=metadata_fungi_soil_R2, permutations=9999)
vegan::vegdist(t(otu_fungi_soil_R2), method="bray") -> dist_otu_fungi_soil_R2

permdisp_otu_fungi_soil_R2<- betadisper(dist_otu_fungi_soil_R2, metadata_fungi_soil_R2$Management)
anova(permdisp_otu_fungi_soil_R2, permutations = 9999)

#soil R6
ps.noncontam_soil_R6 <- subset_samples(ps.noncontam_soil, Growth_Stage%in%c("R6"))
ord_soil_R6 = ordinate(ps.noncontam_soil_R6 , method ="NMDS", distance="bray", try=200)
ord_soil_R6 



# shape is management
NMDS_soil_R6= plot_ordination(ps.noncontam_soil_R6, ord_soil_R6 , color="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_soil_R6)

# get R2 value for plot
otu_fungi_soil_R6<- as.data.frame(otu_table(ps.noncontam_soil_R6))
taxa_fungi_soil_R6<- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_R6)))
metadata_fungi_soil_R6 <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_R6)))
metadata_fungi_soil_R6


adonis(t(otu_fungi_soil_R6) ~ Management, data=metadata_fungi_soil_R6, permutations=9999)
vegan::vegdist(t(otu_fungi_soil_R6), method="bray") -> dist_otu_fungi_soil_R6

permdisp_otu_fungi_soil_R6<- betadisper(dist_otu_fungi_soil_R6, metadata_fungi_soil_R6$Management)
anova(permdisp_otu_fungi_soil_R6, permutations = 9999)


library(ggpubr)
ggarrange(NMDS_soil_V2,NMDS_soil_R2,NMDS_soil_R6,
          widths = c(4.1, 2.7,2.7),
          align = "h", ncol = 3, nrow = 1)

# soil by management
#soil Conventional
ps.noncontam_soil_Conventional <- subset_samples(ps.noncontam_soil, Management%in%c("Conventional"))
ord_soil_Conventional = ordinate(ps.noncontam_soil_Conventional , method ="NMDS", distance="bray", try=200)
ord_soil_Conventional 



# shape is management
NMDS_soil_Conventional= plot_ordination(ps.noncontam_soil_Conventional, ord_soil_Conventional , color="Growth_Stage") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="left")
plot(NMDS_soil_Conventional)

# get R2 value for plot
otu_fungi_soil_Conventional<- as.data.frame(otu_table(ps.noncontam_soil_Conventional))
taxa_fungi_soil_Conventional<- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_Conventional)))
metadata_fungi_soil_Conventional <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_Conventional)))
metadata_fungi_soil_Conventional


adonis(t(otu_fungi_soil_Conventional) ~ Growth_Stage, data=metadata_fungi_soil_Conventional, permutations=9999)
vegan::vegdist(t(otu_fungi_soil_Conventional), method="bray") -> dist_otu_fungi_soil_Conventional

permdisp_otu_fungi_soil_Conventional<- betadisper(dist_otu_fungi_soil_Conventional, metadata_fungi_soil_Conventional$Growth_Stage)
anova(permdisp_otu_fungi_soil_Conventional, permutations = 9999)
#soil No-Till
ps.noncontam_soil_No_Till <- subset_samples(ps.noncontam_soil, Management%in%c("No-Till"))
ord_soil_No_Till = ordinate(ps.noncontam_soil_No_Till , method ="NMDS", distance="bray", try=200)
ord_soil_No_Till 



# shape is Growth_Stage
NMDS_soil_No_Till= plot_ordination(ps.noncontam_soil_No_Till, ord_soil_No_Till , color="Growth_Stage") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_soil_No_Till)

# get No-Till value for plot
otu_fungi_soil_No_Till<- as.data.frame(otu_table(ps.noncontam_soil_No_Till))
taxa_fungi_soil_No_Till<- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_No_Till)))
metadata_fungi_soil_No_Till <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_No_Till)))
metadata_fungi_soil_No_Till


adonis(t(otu_fungi_soil_No_Till) ~ Growth_Stage, data=metadata_fungi_soil_No_Till, permutations=9999)
vegan::vegdist(t(otu_fungi_soil_No_Till), method="bray") -> dist_otu_fungi_soil_No_Till

permdisp_otu_fungi_soil_No_Till<- betadisper(dist_otu_fungi_soil_No_Till, metadata_fungi_soil_No_Till$Growth_Stage)
anova(permdisp_otu_fungi_soil_No_Till, permutations = 9999)

#soil Organic
ps.noncontam_soil_Organic <- subset_samples(ps.noncontam_soil, Management%in%c("Organic"))
ord_soil_Organic = ordinate(ps.noncontam_soil_Organic , method ="NMDS", distance="bray", try=200)
ord_soil_Organic 



# shape is management
NMDS_soil_Organic= plot_ordination(ps.noncontam_soil_Organic, ord_soil_Organic , color="Growth_Stage") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_soil_Organic)

# get R2 value for plot
otu_fungi_soil_Organic<- as.data.frame(otu_table(ps.noncontam_soil_Organic))
taxa_fungi_soil_Organic<- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_Organic)))
metadata_fungi_soil_Organic <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_Organic)))
metadata_fungi_soil_Organic


adonis(t(otu_fungi_soil_Organic) ~ Growth_Stage, data=metadata_fungi_soil_Organic, permutations=9999)
vegan::vegdist(t(otu_fungi_soil_Organic), method="bray") -> dist_otu_fungi_soil_Organic

permdisp_otu_fungi_soil_Organic<- betadisper(dist_otu_fungi_soil_Organic, metadata_fungi_soil_Organic$Growth_Stage)
anova(permdisp_otu_fungi_soil_Organic, permutations = 9999)



ggarrange(NMDS_soil_Conventional,NMDS_soil_No_Till,NMDS_soil_Organic,
          widths = c(4.1, 2.7,2.7),
          align = "h", ncol = 3, nrow = 1)

#Roots
#roots V2
otu_table(ps.noncontam_roots)
ps.noncontam_roots_V2 <- subset_samples(ps.noncontam_roots, Growth_Stage%in%c("V2"))
ord_roots_v2 = ordinate(ps.noncontam_roots_V2 , method ="NMDS", distance="bray", try=200)
ord_roots_v2 



# shape is management
NMDS_roots_V2= plot_ordination(ps.noncontam_roots_V2, ord_roots_v2 , color="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="left")
plot(NMDS_roots_V2)

# get R2 value for plot
otu_fungi_roots_V2<- as.data.frame(otu_table(ps.noncontam_roots_V2))
taxa_fungi_roots_V2<- as.data.frame(as.matrix(tax_table(ps.noncontam_roots_V2)))
metadata_fungi_roots_V2 <- as.data.frame(as.matrix(sample_data(ps.noncontam_roots_V2)))
metadata_fungi_roots_V2


adonis(t(otu_fungi_roots_V2) ~ Management, data=metadata_fungi_roots_V2, permutations=9999)
vegan::vegdist(t(otu_fungi_roots_V2), method="bray") -> dist_otu_fungi_roots_V2

permdisp_otu_fungi_roots_V2<- betadisper(dist_otu_fungi_roots_V2, metadata_fungi_roots_V2$Management)
anova(permdisp_otu_fungi_roots_V2, permutations = 9999)
#roots R2
ps.noncontam_roots_R2 <- subset_samples(ps.noncontam_roots, Growth_Stage%in%c("R2"))
ord_roots_R2 = ordinate(ps.noncontam_roots_R2 , method ="NMDS", distance="bray", try=200)
ord_roots_R2 



# shape is management
NMDS_roots_R2= plot_ordination(ps.noncontam_roots_R2, ord_roots_R2 , color="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_roots_R2)

# get R2 value for plot
otu_fungi_roots_R2<- as.data.frame(otu_table(ps.noncontam_roots_R2))
taxa_fungi_roots_R2<- as.data.frame(as.matrix(tax_table(ps.noncontam_roots_R2)))
metadata_fungi_roots_R2 <- as.data.frame(as.matrix(sample_data(ps.noncontam_roots_R2)))
metadata_fungi_roots_R2


adonis(t(otu_fungi_roots_R2) ~ Management, data=metadata_fungi_roots_R2, permutations=9999)
vegan::vegdist(t(otu_fungi_roots_R2), method="bray") -> dist_otu_fungi_roots_R2

permdisp_otu_fungi_roots_R2<- betadisper(dist_otu_fungi_roots_R2, metadata_fungi_roots_R2$Management)
anova(permdisp_otu_fungi_roots_R2, permutations = 9999)

#roots R6
ps.noncontam_roots_R6 <- subset_samples(ps.noncontam_roots, Growth_Stage%in%c("R6"))
ord_roots_R6 = ordinate(ps.noncontam_roots_R6 , method ="NMDS", distance="bray", try=200)
ord_roots_R6 



# shape is management
NMDS_roots_R6= plot_ordination(ps.noncontam_roots_R6, ord_roots_R6 , color="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_roots_R6)

# get R2 value for plot
otu_fungi_roots_R6<- as.data.frame(otu_table(ps.noncontam_roots_R6))
taxa_fungi_roots_R6<- as.data.frame(as.matrix(tax_table(ps.noncontam_roots_R6)))
metadata_fungi_roots_R6 <- as.data.frame(as.matrix(sample_data(ps.noncontam_roots_R6)))
metadata_fungi_roots_R6


adonis(t(otu_fungi_roots_R6) ~ Management, data=metadata_fungi_roots_R6, permutations=9999)
vegan::vegdist(t(otu_fungi_roots_R6), method="bray") -> dist_otu_fungi_roots_R6

permdisp_otu_fungi_roots_R6<- betadisper(dist_otu_fungi_roots_R6, metadata_fungi_roots_R6$Management)
anova(permdisp_otu_fungi_roots_R6, permutations = 9999)


ggarrange(NMDS_roots_V2,NMDS_roots_R2,NMDS_roots_R6,
          widths = c(4.1, 2.7,2.7),
          align = "h", ncol = 3, nrow = 1)

# roots by management
#roots Conventional
ps.noncontam_roots_Conventional <- subset_samples(ps.noncontam_roots, Management%in%c("Conventional"))
ord_roots_Conventional = ordinate(ps.noncontam_roots_Conventional , method ="NMDS", distance="bray", try=200)
ord_roots_Conventional 



# shape is management
NMDS_roots_Conventional= plot_ordination(ps.noncontam_roots_Conventional, ord_roots_Conventional , color="Growth_Stage") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="left")
plot(NMDS_roots_Conventional)

# get R2 value for plot
otu_fungi_roots_Conventional<- as.data.frame(otu_table(ps.noncontam_roots_Conventional))
taxa_fungi_roots_Conventional<- as.data.frame(as.matrix(tax_table(ps.noncontam_roots_Conventional)))
metadata_fungi_roots_Conventional <- as.data.frame(as.matrix(sample_data(ps.noncontam_roots_Conventional)))
metadata_fungi_roots_Conventional


adonis(t(otu_fungi_roots_Conventional) ~ Growth_Stage, data=metadata_fungi_roots_Conventional, permutations=9999)
vegan::vegdist(t(otu_fungi_roots_Conventional), method="bray") -> dist_otu_fungi_roots_Conventional

permdisp_otu_fungi_roots_Conventional<- betadisper(dist_otu_fungi_roots_Conventional, metadata_fungi_roots_Conventional$Growth_Stage)
anova(permdisp_otu_fungi_roots_Conventional, permutations = 9999)
#roots No-Till
ps.noncontam_roots_No_Till <- subset_samples(ps.noncontam_roots, Management%in%c("No-Till"))
ord_roots_No_Till = ordinate(ps.noncontam_roots_No_Till , method ="NMDS", distance="bray", try=200)
ord_roots_No_Till 



# shape is Growth_Stage
NMDS_roots_No_Till= plot_ordination(ps.noncontam_roots_No_Till, ord_roots_No_Till , color="Growth_Stage") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_roots_No_Till)

# get No-Till value for plot
otu_fungi_roots_No_Till<- as.data.frame(otu_table(ps.noncontam_roots_No_Till))
taxa_fungi_roots_No_Till<- as.data.frame(as.matrix(tax_table(ps.noncontam_roots_No_Till)))
metadata_fungi_roots_No_Till <- as.data.frame(as.matrix(sample_data(ps.noncontam_roots_No_Till)))
metadata_fungi_roots_No_Till


adonis(t(otu_fungi_roots_No_Till) ~ Growth_Stage, data=metadata_fungi_roots_No_Till, permutations=9999)
vegan::vegdist(t(otu_fungi_roots_No_Till), method="bray") -> dist_otu_fungi_roots_No_Till

permdisp_otu_fungi_roots_No_Till<- betadisper(dist_otu_fungi_roots_No_Till, metadata_fungi_roots_No_Till$Growth_Stage)
anova(permdisp_otu_fungi_roots_No_Till, permutations = 9999)

#roots Organic
ps.noncontam_roots_Organic <- subset_samples(ps.noncontam_roots, Management%in%c("Organic"))
ord_roots_Organic = ordinate(ps.noncontam_roots_Organic , method ="NMDS", distance="bray", try=200)
ord_roots_Organic 



# shape is management
NMDS_roots_Organic= plot_ordination(ps.noncontam_roots_Organic, ord_roots_Organic , color="Growth_Stage") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_roots_Organic)

# get R2 value for plot
otu_fungi_roots_Organic<- as.data.frame(otu_table(ps.noncontam_roots_Organic))
taxa_fungi_roots_Organic<- as.data.frame(as.matrix(tax_table(ps.noncontam_roots_Organic)))
metadata_fungi_roots_Organic <- as.data.frame(as.matrix(sample_data(ps.noncontam_roots_Organic)))
metadata_fungi_roots_Organic


adonis(t(otu_fungi_roots_Organic) ~ Growth_Stage, data=metadata_fungi_roots_Organic, permutations=9999)
vegan::vegdist(t(otu_fungi_roots_Organic), method="bray") -> dist_otu_fungi_roots_Organic

permdisp_otu_fungi_roots_Organic<- betadisper(dist_otu_fungi_roots_Organic, metadata_fungi_roots_Organic$Growth_Stage)
anova(permdisp_otu_fungi_roots_Organic, permutations = 9999)



ggarrange(NMDS_roots_Conventional,NMDS_roots_No_Till,NMDS_roots_Organic,
          widths = c(4.1, 2.7,2.7),
          align = "h", ncol = 3, nrow = 1)

#Stems
#stems V2
otu_table(ps.noncontam_stems)
ps.noncontam_stems_V2 <- subset_samples(ps.noncontam_stems, Growth_Stage%in%c("V2"))
ord_stems_v2 = ordinate(ps.noncontam_stems_V2 , method ="NMDS", distance="bray", try=200)
ord_stems_v2 



# shape is management
NMDS_stems_V2= plot_ordination(ps.noncontam_stems_V2, ord_stems_v2 , color="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="left")
plot(NMDS_stems_V2)

# get R2 value for plot
otu_fungi_stems_V2<- as.data.frame(otu_table(ps.noncontam_stems_V2))
taxa_fungi_stems_V2<- as.data.frame(as.matrix(tax_table(ps.noncontam_stems_V2)))
metadata_fungi_stems_V2 <- as.data.frame(as.matrix(sample_data(ps.noncontam_stems_V2)))
metadata_fungi_stems_V2


adonis(t(otu_fungi_stems_V2) ~ Management, data=metadata_fungi_stems_V2, permutations=9999)
vegan::vegdist(t(otu_fungi_stems_V2), method="bray") -> dist_otu_fungi_stems_V2

permdisp_otu_fungi_stems_V2<- betadisper(dist_otu_fungi_stems_V2, metadata_fungi_stems_V2$Management)
anova(permdisp_otu_fungi_stems_V2, permutations = 9999)
#stems R2
ps.noncontam_stems_R2 <- subset_samples(ps.noncontam_stems, Growth_Stage%in%c("R2"))
ord_stems_R2 = ordinate(ps.noncontam_stems_R2 , method ="NMDS", distance="bray", try=200)
ord_stems_R2 



# shape is management
NMDS_stems_R2= plot_ordination(ps.noncontam_stems_R2, ord_stems_R2 , color="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_stems_R2)

# get R2 value for plot
otu_fungi_stems_R2<- as.data.frame(otu_table(ps.noncontam_stems_R2))
taxa_fungi_stems_R2<- as.data.frame(as.matrix(tax_table(ps.noncontam_stems_R2)))
metadata_fungi_stems_R2 <- as.data.frame(as.matrix(sample_data(ps.noncontam_stems_R2)))
metadata_fungi_stems_R2


adonis(t(otu_fungi_stems_R2) ~ Management, data=metadata_fungi_stems_R2, permutations=9999)
vegan::vegdist(t(otu_fungi_stems_R2), method="bray") -> dist_otu_fungi_stems_R2

permdisp_otu_fungi_stems_R2<- betadisper(dist_otu_fungi_stems_R2, metadata_fungi_stems_R2$Management)
anova(permdisp_otu_fungi_stems_R2, permutations = 9999)

#stems R6
ps.noncontam_stems_R6 <- subset_samples(ps.noncontam_stems, Growth_Stage%in%c("R6"))
ord_stems_R6 = ordinate(ps.noncontam_stems_R6 , method ="NMDS", distance="bray", try=200)
ord_stems_R6 



# shape is management
NMDS_stems_R6= plot_ordination(ps.noncontam_stems_R6, ord_stems_R6 , color="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_stems_R6)

# get R2 value for plot
otu_fungi_stems_R6<- as.data.frame(otu_table(ps.noncontam_stems_R6))
taxa_fungi_stems_R6<- as.data.frame(as.matrix(tax_table(ps.noncontam_stems_R6)))
metadata_fungi_stems_R6 <- as.data.frame(as.matrix(sample_data(ps.noncontam_stems_R6)))
metadata_fungi_stems_R6


adonis(t(otu_fungi_stems_R6) ~ Management, data=metadata_fungi_stems_R6, permutations=9999)
vegan::vegdist(t(otu_fungi_stems_R6), method="bray") -> dist_otu_fungi_stems_R6

permdisp_otu_fungi_stems_R6<- betadisper(dist_otu_fungi_stems_R6, metadata_fungi_stems_R6$Management)
anova(permdisp_otu_fungi_stems_R6, permutations = 9999)


ggarrange(NMDS_stems_V2,NMDS_stems_R2,NMDS_stems_R6,
          widths = c(4.1, 2.7,2.7),
          align = "h", ncol = 3, nrow = 1)

# stems by management
#stems Conventional
ps.noncontam_stems_Conventional <- subset_samples(ps.noncontam_stems, Management%in%c("Conventional"))
ord_stems_Conventional = ordinate(ps.noncontam_stems_Conventional , method ="NMDS", distance="bray", try=200)
ord_stems_Conventional 



# shape is management
NMDS_stems_Conventional= plot_ordination(ps.noncontam_stems_Conventional, ord_stems_Conventional , color="Growth_Stage") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="left")
plot(NMDS_stems_Conventional)

# get R2 value for plot
otu_fungi_stems_Conventional<- as.data.frame(otu_table(ps.noncontam_stems_Conventional))
taxa_fungi_stems_Conventional<- as.data.frame(as.matrix(tax_table(ps.noncontam_stems_Conventional)))
metadata_fungi_stems_Conventional <- as.data.frame(as.matrix(sample_data(ps.noncontam_stems_Conventional)))
metadata_fungi_stems_Conventional


adonis(t(otu_fungi_stems_Conventional) ~ Growth_Stage, data=metadata_fungi_stems_Conventional, permutations=9999)
vegan::vegdist(t(otu_fungi_stems_Conventional), method="bray") -> dist_otu_fungi_stems_Conventional

permdisp_otu_fungi_stems_Conventional<- betadisper(dist_otu_fungi_stems_Conventional, metadata_fungi_stems_Conventional$Growth_Stage)
anova(permdisp_otu_fungi_stems_Conventional, permutations = 9999)
#stems No-Till
ps.noncontam_stems_No_Till <- subset_samples(ps.noncontam_stems, Management%in%c("No-Till"))
ord_stems_No_Till = ordinate(ps.noncontam_stems_No_Till , method ="NMDS", distance="bray", try=200)
ord_stems_No_Till 



# shape is Growth_Stage
NMDS_stems_No_Till= plot_ordination(ps.noncontam_stems_No_Till, ord_stems_No_Till , color="Growth_Stage") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_stems_No_Till)

# get No-Till value for plot
otu_fungi_stems_No_Till<- as.data.frame(otu_table(ps.noncontam_stems_No_Till))
taxa_fungi_stems_No_Till<- as.data.frame(as.matrix(tax_table(ps.noncontam_stems_No_Till)))
metadata_fungi_stems_No_Till <- as.data.frame(as.matrix(sample_data(ps.noncontam_stems_No_Till)))
metadata_fungi_stems_No_Till


adonis(t(otu_fungi_stems_No_Till) ~ Growth_Stage, data=metadata_fungi_stems_No_Till, permutations=9999)
vegan::vegdist(t(otu_fungi_stems_No_Till), method="bray") -> dist_otu_fungi_stems_No_Till

permdisp_otu_fungi_stems_No_Till<- betadisper(dist_otu_fungi_stems_No_Till, metadata_fungi_stems_No_Till$Growth_Stage)
anova(permdisp_otu_fungi_stems_No_Till, permutations = 9999)

#stems Organic
ps.noncontam_stems_Organic <- subset_samples(ps.noncontam_stems, Management%in%c("Organic"))
ord_stems_Organic = ordinate(ps.noncontam_stems_Organic , method ="NMDS", distance="bray", try=200)
ord_stems_Organic 



# shape is management
NMDS_stems_Organic= plot_ordination(ps.noncontam_stems_Organic, ord_stems_Organic , color="Growth_Stage") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_stems_Organic)

# get R2 value for plot
otu_fungi_stems_Organic<- as.data.frame(otu_table(ps.noncontam_stems_Organic))
taxa_fungi_stems_Organic<- as.data.frame(as.matrix(tax_table(ps.noncontam_stems_Organic)))
metadata_fungi_stems_Organic <- as.data.frame(as.matrix(sample_data(ps.noncontam_stems_Organic)))
metadata_fungi_stems_Organic


adonis(t(otu_fungi_stems_Organic) ~ Growth_Stage, data=metadata_fungi_stems_Organic, permutations=9999)
vegan::vegdist(t(otu_fungi_stems_Organic), method="bray") -> dist_otu_fungi_stems_Organic

permdisp_otu_fungi_stems_Organic<- betadisper(dist_otu_fungi_stems_Organic, metadata_fungi_stems_Organic$Growth_Stage)
anova(permdisp_otu_fungi_stems_Organic, permutations = 9999)



ggarrange(NMDS_stems_Conventional,NMDS_stems_No_Till,NMDS_stems_Organic,
          widths = c(4.1, 2.7,2.7),
          align = "h", ncol = 3, nrow = 1)

#leaves
#leaves V2
otu_table(ps.noncontam_leaves_obj1)
ps.noncontam_leaves_obj1_V2 <- subset_samples(ps.noncontam_leaves_obj1, Growth_Stage%in%c("V2"))
ord_leaves_v2 = ordinate(ps.noncontam_leaves_obj1_V2 , method ="NMDS", distance="bray", try=200)
ord_leaves_v2 



# shape is management
NMDS_leaves_V2= plot_ordination(ps.noncontam_leaves_obj1_V2, ord_leaves_v2 , color="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="left")
plot(NMDS_leaves_V2)

# get R2 value for plot
otu_fungi_leaves_V2<- as.data.frame(otu_table(ps.noncontam_leaves_obj1_V2))
taxa_fungi_leaves_V2<- as.data.frame(as.matrix(tax_table(ps.noncontam_leaves_obj1_V2)))
metadata_fungi_leaves_V2 <- as.data.frame(as.matrix(sample_data(ps.noncontam_leaves_obj1_V2)))
metadata_fungi_leaves_V2


adonis(t(otu_fungi_leaves_V2) ~ Management, data=metadata_fungi_leaves_V2, permutations=9999)
vegan::vegdist(t(otu_fungi_leaves_V2), method="bray") -> dist_otu_fungi_leaves_V2

permdisp_otu_fungi_leaves_V2<- betadisper(dist_otu_fungi_leaves_V2, metadata_fungi_leaves_V2$Management)
anova(permdisp_otu_fungi_leaves_V2, permutations = 9999)
#leaves R2
ps.noncontam_leaves_obj1_R2 <- subset_samples(ps.noncontam_leaves_obj1, Growth_Stage%in%c("R2"))
otu_table(ps.noncontam_leaves_obj1_R2)

ord_leaves_R2 = ordinate(ps.noncontam_leaves_obj1_R2 , method ="NMDS", distance="bray", try=200)
ord_leaves_R2 



# shape is management
NMDS_leaves_R2= plot_ordination(ps.noncontam_leaves_obj1_R2, ord_leaves_R2 , color="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_leaves_R2)

# get R2 value for plot
otu_fungi_leaves_R2<- as.data.frame(otu_table(ps.noncontam_leaves_obj1_R2))
taxa_fungi_leaves_R2<- as.data.frame(as.matrix(tax_table(ps.noncontam_leaves_obj1_R2)))
metadata_fungi_leaves_R2 <- as.data.frame(as.matrix(sample_data(ps.noncontam_leaves_obj1_R2)))
metadata_fungi_leaves_R2


adonis(t(otu_fungi_leaves_R2) ~ Management, data=metadata_fungi_leaves_R2, permutations=9999)
vegan::vegdist(t(otu_fungi_leaves_R2), method="bray") -> dist_otu_fungi_leaves_R2

permdisp_otu_fungi_leaves_R2<- betadisper(dist_otu_fungi_leaves_R2, metadata_fungi_leaves_R2$Management)
anova(permdisp_otu_fungi_leaves_R2, permutations = 9999)

#leaves R6
ps.noncontam_leaves_obj1_R6 <- subset_samples(ps.noncontam_leaves_obj1, Growth_Stage%in%c("R6"))
ord_leaves_R6 = ordinate(ps.noncontam_leaves_obj1_R6 , method ="NMDS", distance="bray", try=200)
ord_leaves_R6 



# shape is management
NMDS_leaves_R6= plot_ordination(ps.noncontam_leaves_obj1_R6, ord_leaves_R6 , color="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_leaves_R6)

# get R2 value for plot
otu_fungi_leaves_R6<- as.data.frame(otu_table(ps.noncontam_leaves_obj1_R6))
taxa_fungi_leaves_R6<- as.data.frame(as.matrix(tax_table(ps.noncontam_leaves_obj1_R6)))
metadata_fungi_leaves_R6 <- as.data.frame(as.matrix(sample_data(ps.noncontam_leaves_obj1_R6)))
metadata_fungi_leaves_R6


adonis(t(otu_fungi_leaves_R6) ~ Management, data=metadata_fungi_leaves_R6, permutations=9999)
vegan::vegdist(t(otu_fungi_leaves_R6), method="bray") -> dist_otu_fungi_leaves_R6

permdisp_otu_fungi_leaves_R6<- betadisper(dist_otu_fungi_leaves_R6, metadata_fungi_leaves_R6$Management)
anova(permdisp_otu_fungi_leaves_R6, permutations = 9999)


ggarrange(NMDS_leaves_V2,NMDS_leaves_R2,NMDS_leaves_R6,
          widths = c(4.1, 2.7,2.7),
          align = "h", ncol = 3, nrow = 1)

# leaves by management
#leaves Conventional
ps.noncontam_leaves_obj1_Conventional <- subset_samples(ps.noncontam_leaves_obj1, Management%in%c("Conventional"))
ord_leaves_Conventional = ordinate(ps.noncontam_leaves_obj1_Conventional , method ="NMDS", distance="bray", try=200)
ord_leaves_Conventional 



# shape is management
NMDS_leaves_Conventional= plot_ordination(ps.noncontam_leaves_obj1_Conventional, ord_leaves_Conventional , color="Growth_Stage") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="left")
plot(NMDS_leaves_Conventional)

# get R2 value for plot
otu_fungi_leaves_Conventional<- as.data.frame(otu_table(ps.noncontam_leaves_obj1_Conventional))
taxa_fungi_leaves_Conventional<- as.data.frame(as.matrix(tax_table(ps.noncontam_leaves_obj1_Conventional)))
metadata_fungi_leaves_Conventional <- as.data.frame(as.matrix(sample_data(ps.noncontam_leaves_obj1_Conventional)))
metadata_fungi_leaves_Conventional


adonis(t(otu_fungi_leaves_Conventional) ~ Growth_Stage, data=metadata_fungi_leaves_Conventional, permutations=9999)
vegan::vegdist(t(otu_fungi_leaves_Conventional), method="bray") -> dist_otu_fungi_leaves_Conventional

permdisp_otu_fungi_leaves_Conventional<- betadisper(dist_otu_fungi_leaves_Conventional, metadata_fungi_leaves_Conventional$Growth_Stage)
anova(permdisp_otu_fungi_leaves_Conventional, permutations = 9999)
#leaves No-Till
ps.noncontam_leaves_obj1_No_Till <- subset_samples(ps.noncontam_leaves_obj1, Management%in%c("No-Till"))
ord_leaves_No_Till = ordinate(ps.noncontam_leaves_obj1_No_Till , method ="NMDS", distance="bray", try=200)
ord_leaves_No_Till 



# shape is Growth_Stage
NMDS_leaves_No_Till= plot_ordination(ps.noncontam_leaves_obj1_No_Till, ord_leaves_No_Till , color="Growth_Stage") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_leaves_No_Till)

# get No-Till value for plot
otu_fungi_leaves_No_Till<- as.data.frame(otu_table(ps.noncontam_leaves_obj1_No_Till))
taxa_fungi_leaves_No_Till<- as.data.frame(as.matrix(tax_table(ps.noncontam_leaves_obj1_No_Till)))
metadata_fungi_leaves_No_Till <- as.data.frame(as.matrix(sample_data(ps.noncontam_leaves_obj1_No_Till)))
metadata_fungi_leaves_No_Till


adonis(t(otu_fungi_leaves_No_Till) ~ Growth_Stage, data=metadata_fungi_leaves_No_Till, permutations=9999)
vegan::vegdist(t(otu_fungi_leaves_No_Till), method="bray") -> dist_otu_fungi_leaves_No_Till

permdisp_otu_fungi_leaves_No_Till<- betadisper(dist_otu_fungi_leaves_No_Till, metadata_fungi_leaves_No_Till$Growth_Stage)
anova(permdisp_otu_fungi_leaves_No_Till, permutations = 9999)

#leaves Organic
ps.noncontam_leaves_obj1_Organic <- subset_samples(ps.noncontam_leaves_obj1, Management%in%c("Organic"))
ord_leaves_Organic = ordinate(ps.noncontam_leaves_obj1_Organic , method ="NMDS", distance="bray", try=200)
ord_leaves_Organic 



# shape is management
NMDS_leaves_Organic= plot_ordination(ps.noncontam_leaves_obj1_Organic, ord_leaves_Organic , color="Growth_Stage") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_leaves_Organic)

# get R2 value for plot
otu_fungi_leaves_Organic<- as.data.frame(otu_table(ps.noncontam_leaves_obj1_Organic))
taxa_fungi_leaves_Organic<- as.data.frame(as.matrix(tax_table(ps.noncontam_leaves_obj1_Organic)))
metadata_fungi_leaves_Organic <- as.data.frame(as.matrix(sample_data(ps.noncontam_leaves_obj1_Organic)))
metadata_fungi_leaves_Organic


adonis(t(otu_fungi_leaves_Organic) ~ Growth_Stage, data=metadata_fungi_leaves_Organic, permutations=9999)
vegan::vegdist(t(otu_fungi_leaves_Organic), method="bray") -> dist_otu_fungi_leaves_Organic

permdisp_otu_fungi_leaves_Organic<- betadisper(dist_otu_fungi_leaves_Organic, metadata_fungi_leaves_Organic$Growth_Stage)
anova(permdisp_otu_fungi_leaves_Organic, permutations = 9999)



ggarrange(NMDS_leaves_Conventional,NMDS_leaves_No_Till,NMDS_leaves_Organic,
          widths = c(4.1, 2.7,2.7),
          align = "h", ncol = 3, nrow = 1)

### PERMANOVA (TABLE 3)-------

# creating vegan objects for use in permanova
# will analyze each sample origin individiually to parse out variation by growth stage and management 

#soil
otu_fungi_soil <- as.data.frame(otu_table(ps.noncontam_soil_obj1))
taxa_fungi_soil <- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_obj1)))
metadata_fungi_soil <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_obj1)))


#roots
otu_fungi_roots <- as.data.frame(otu_table(ps.noncontam_roots_obj1))
taxa_fungi_roots <- as.data.frame(as.matrix(tax_table(ps.noncontam_roots_obj1)))
metadata_fungi_roots <- as.data.frame(as.matrix(sample_data(ps.noncontam_roots_obj1)))

#stems
otu_fungi_stems <- as.data.frame(otu_table(ps.noncontam_stems_obj1))
taxa_fungi_stems <- as.data.frame(as.matrix(tax_table(ps.noncontam_stems_obj1)))
metadata_fungi_stems <- as.data.frame(as.matrix(sample_data(ps.noncontam_stems_obj1)))


#leaves
otu_fungi_leaves <- as.data.frame(otu_table(ps.noncontam_leaves_obj1))
taxa_fungi_leaves <- as.data.frame(as.matrix(tax_table(ps.noncontam_leaves_obj1)))
metadata_fungi_leaves <- as.data.frame(as.matrix(sample_data(ps.noncontam_leaves_obj1)))
#below permanova
options(scipen = 999) 
library("vegan")
library("RVAideMemoire")


model.matrix(~ Growth_Stage * Management, data=metadata_fungi_soil)
model.matrix(~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_fungi_soil)

adonis(t(otu_fungi_soil) ~ Growth_Stage * Management, data=metadata_fungi_soil, permutations=9999) # by = "margin"

adonis(t(otu_fungi_soil) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_fungi_soil, permutations=9999) 
adonis(t(otu_fungi_soil) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_fungi_soil, permutations=9999)


# soil betadisper
vegan::vegdist(t(otu_fungi_soil), method="bray") -> dist_otu_fungi_soil

permdisp_otu_fungi_soil_M <- betadisper(dist_otu_fungi_soil, metadata_fungi_soil$Management)
permdisp_otu_fungi_soil_GS<- betadisper(dist_otu_fungi_soil, metadata_fungi_soil$Growth_Stage)

anova(permdisp_otu_fungi_soil_M, permutations = 9999)
permutest(permdisp_otu_fungi_soil_M, permutations = 9999, pairwise = T)
plot(permdisp_otu_fungi_soil_M)
plot(TukeyHSD(permdisp_otu_fungi_soil_M), las=1)
boxplot(permdisp_otu_fungi_soil_M)

anova(permdisp_otu_fungi_soil_GS, permutations = 9999)
permutest(permdisp_otu_fungi_soil_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_fungi_soil_GS)
plot(TukeyHSD(permdisp_otu_fungi_soil_GS), las=1)
boxplot(permdisp_otu_fungi_soil_GS)

#roots permanova
options(scipen = 999) 
library("vegan")
library("RVAideMemoire")


model.matrix(~ Growth_Stage * Management, data=metadata_fungi_roots)
model.matrix(~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_fungi_roots)

adonis(t(otu_fungi_roots) ~ Growth_Stage * Management, data=metadata_fungi_roots, permutations=9999) # by = "margin"

adonis(t(otu_fungi_roots) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_fungi_roots, permutations=9999) 
adonis(t(otu_fungi_roots) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_fungi_roots, permutations=9999)


# roots betadisper
vegan::vegdist(t(otu_fungi_roots), method="bray") -> dist_otu_fungi_roots

permdisp_otu_fungi_roots_M <- betadisper(dist_otu_fungi_roots, metadata_fungi_roots$Management)
permdisp_otu_fungi_roots_GS<- betadisper(dist_otu_fungi_roots, metadata_fungi_roots$Growth_Stage)

anova(permdisp_otu_fungi_roots_M, permutations = 9999)
permutest(permdisp_otu_fungi_roots_M, permutations = 9999, pairwise = T)
plot(permdisp_otu_fungi_roots_M)
plot(TukeyHSD(permdisp_otu_fungi_roots_M), las=1)
boxplot(permdisp_otu_fungi_roots_M)

anova(permdisp_otu_fungi_roots_GS, permutations = 9999)
permutest(permdisp_otu_fungi_roots_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_fungi_roots_GS)
plot(TukeyHSD(permdisp_otu_fungi_roots_GS), las=1)
boxplot(permdisp_otu_fungi_roots_GS)


#stems permanova
options(scipen = 999) 
library("vegan")
library("RVAideMemoire")


model.matrix(~ Growth_Stage * Management, data=metadata_fungi_stems)
model.matrix(~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_fungi_stems)

adonis(t(otu_fungi_stems) ~ Growth_Stage * Management, data=metadata_fungi_stems, permutations=9999) # by = "margin"

adonis(t(otu_fungi_stems) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_fungi_stems, permutations=9999) 
adonis(t(otu_fungi_stems) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_fungi_stems, permutations=9999)


# stems betadisper
vegan::vegdist(t(otu_fungi_stems), method="bray") -> dist_otu_fungi_stems

permdisp_otu_fungi_stems_M <- betadisper(dist_otu_fungi_stems, metadata_fungi_stems$Management)
permdisp_otu_fungi_stems_GS<- betadisper(dist_otu_fungi_stems, metadata_fungi_stems$Growth_Stage)

anova(permdisp_otu_fungi_stems_M, permutations = 9999)
permutest(permdisp_otu_fungi_stems_M, permutations = 9999, pairwise = T)
plot(permdisp_otu_fungi_stems_M)
plot(TukeyHSD(permdisp_otu_fungi_stems_M), las=1)
boxplot(permdisp_otu_fungi_stems_M)

anova(permdisp_otu_fungi_stems_GS, permutations = 9999)
permutest(permdisp_otu_fungi_stems_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_fungi_stems_GS)
plot(TukeyHSD(permdisp_otu_fungi_stems_GS), las=1)
boxplot(permdisp_otu_fungi_stems_GS)


#leaves permanova
options(scipen = 999) 
library("vegan")
library("RVAideMemoire")


model.matrix(~ Growth_Stage * Management, data=metadata_fungi_leaves)
model.matrix(~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_fungi_leaves)

adonis(t(otu_fungi_leaves) ~ Growth_Stage * Management, data=metadata_fungi_leaves, permutations=9999) # by = "margin"

adonis(t(otu_fungi_leaves) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_fungi_leaves, permutations=9999) 
adonis(t(otu_fungi_leaves) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_fungi_leaves, permutations=9999)


# leaves betadisper
vegan::vegdist(t(otu_fungi_leaves), method="bray") -> dist_otu_fungi_leaves

permdisp_otu_fungi_leaves_M <- betadisper(dist_otu_fungi_leaves, metadata_fungi_leaves$Management)
permdisp_otu_fungi_leaves_GS<- betadisper(dist_otu_fungi_leaves, metadata_fungi_leaves$Growth_Stage)

anova(permdisp_otu_fungi_leaves_M, permutations = 9999)
permutest(permdisp_otu_fungi_leaves_M, permutations = 9999, pairwise = T)
plot(permdisp_otu_fungi_leaves_M)
plot(TukeyHSD(permdisp_otu_fungi_leaves_M), las=1)
boxplot(permdisp_otu_fungi_leaves_M)

anova(permdisp_otu_fungi_leaves_GS, permutations = 9999)
permutest(permdisp_otu_fungi_leaves_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_fungi_leaves_GS)
plot(TukeyHSD(permdisp_otu_fungi_leaves_GS), las=1)
boxplot(permdisp_otu_fungi_leaves_GS)


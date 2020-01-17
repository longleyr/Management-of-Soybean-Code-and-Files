#Title: Crop Management Impacts the Soy (Glycine max) Microbiome
#Authors: Reid Longley, Zachary A. Noel, Gian Maria Niccolo Benucci, Martin I. Chilvers, Frances Trail, Gregory Bonito 
#Journal: Submitted to Frontiers in Microbiology

# Network Analysis Script - Reid Longley

set.seed(9294)
# >>> loading packages --------------------------------------------------------------
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2")
library("Biostrings")
library("SpiecEasi"); packageVersion("SpiecEasi")
library("igraph"); packageVersion("igraph")

# fungi phyloseq object
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

physeq_object_Fungi <- phyloseq(ITS_otus_phy,
                                ITS_metadata_phy,
                                ITS_taxonomy_phy,
                                ITS_sequences)
physeq_object_Fungi
tax_table(physeq_object_Fungi)
sample_data(physeq_object_Fungi)

###formatting taxonomy------------------------
colnames(tax_table(physeq_object_Fungi)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
tax_table(physeq_object_Fungi)


tax_table(physeq_object_Fungi)[, "Kingdom"] <- gsub("d:", "", tax_table(physeq_object_Fungi)[, "Kingdom"])
tax_table(physeq_object_Fungi)[, "Phylum"] <- gsub("p:", "", tax_table(physeq_object_Fungi)[, "Phylum"])
tax_table(physeq_object_Fungi)[, "Class"] <- gsub("c:", "", tax_table(physeq_object_Fungi)[, "Class"])
tax_table(physeq_object_Fungi)[, "Order"] <- gsub("o:", "", tax_table(physeq_object_Fungi)[, "Order"])
tax_table(physeq_object_Fungi)[, "Family"] <- gsub("f:", "", tax_table(physeq_object_Fungi)[, "Family"])
tax_table(physeq_object_Fungi)[, "Genus"] <- gsub("g:", "", tax_table(physeq_object_Fungi)[, "Genus"])
tax_table(physeq_object_Fungi)[, "Species"] <- gsub("s:", "", tax_table(physeq_object_Fungi)[, "Species"])
tax_table(physeq_object_Fungi)


#remove chloropast,mitochondria,cyanobacteria
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Phylum!="Chloroplast")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Class!="Chloroplast")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Order!="Chloroplast")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Family!="Chloroplast")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Genus!="Chloroplast")
tax_table(physeq_object_Fungi)

physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Phylum!="Mitochondria")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Class!="Mitochondria")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Order!="Mitochondria")

physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Family!="Mitochondria")
physeq_object_Fungi <- subset_taxa(physeq_object_Fungi, Genus!="Mitochondria")
tax_table(physeq_object_Fungi)





###remove contaminants--------------------------------------------------


#removing samples from soil dataset which were noted to not dry properly, and had mold growing in envelope. were analyzed, but only produced a few otus
#non organic samples on this list are being switched with the matching control sample from obj 2 
otu_table(physeq_object_Fungi) <- subset(otu_table(physeq_object_Fungi),
                                         select = -c(T4R5AR2S,T4R6BR2S,T4R2CR2S,T4R5CR2S,T2R6BR2S))




# must split by sample source so that contaminants are removed from each individual run
physeq_object_Fungi_stems <- subset_samples(physeq_object_Fungi, origin%in%c("stem"))

physeq_object_Fungi_leaves <- subset_samples(physeq_object_Fungi,origin%in%c("leaves"))



physeq_object_Fungi_roots <- subset_samples(physeq_object_Fungi,origin%in%c("root"))
physeq_object_Fungi_soil <- subset_samples(physeq_object_Fungi,origin%in%c("soil"))
sample_data(physeq_object_Fungi_soil)


library(devtools)
library(processx)
devtools::install_github("benjjneb/decontam", force = TRUE)
library(decontam)

set.seed(9294)
#soil
# check library size distribution
write.csv(sample_data(physeq_object_Fungi_soil), file = "sample_check1_soil.csv")
df_soil <- as.data.frame(sample_data(physeq_object_Fungi_soil)) # Put sample_data into a ggplot-friendly data.frame
df_soil$LibrarySize_soil <- sample_sums(physeq_object_Fungi_soil)

df_soil <- df_soil[order(df_soil$LibrarySize_soil),]
df_soil$Index <- seq(nrow(df_soil))
write.csv(df_soil, file = "rank_sums_soil.csv")
ggplot(data=df_soil, aes(x=Index, y=LibrarySize_soil, color=Sample_or_Control)) + geom_point()

# filter by prevelance 
sample_data(physeq_object_Fungi_soil)$is.neg <- sample_data(physeq_object_Fungi_soil)$Sample_or_Control == "Control Sample"
contamdf.prev_soil <- isContaminant(physeq_object_Fungi_soil, method="prevalence", neg="is.neg")
table(contamdf.prev_soil$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_soil <- transform_sample_counts(physeq_object_Fungi_soil, function(abund) 1*(abund>0))
ps.pa.neg_soil <- prune_samples(sample_data(ps.pa_soil)$Sample_or_Control == "Control Sample", ps.pa_soil)
ps.pa.pos_soil <- prune_samples(sample_data(ps.pa_soil)$Sample_or_Control == "True Sample", ps.pa_soil)
# Make data.frame of prevalence in positive and negative samples
df.pa_soil <- data.frame(pa.pos_soil=taxa_sums(ps.pa.pos_soil), pa.neg_soil=taxa_sums(ps.pa.neg_soil),
                         contaminant=contamdf.prev_soil$contaminant)
ggplot(data=df.pa_soil, aes(x=pa.neg_soil, y=pa.pos_soil, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_fungi_soil <- prune_taxa(!contamdf.prev_soil$contaminant, physeq_object_Fungi_soil)
# with contaminants removed
otu_table(ps.noncontam_fungi_soil)


#roots

# check library size distribution
write.csv(sample_data(physeq_object_Fungi_roots), file = "sample_check1_roots.csv")
df_roots <- as.data.frame(sample_data(physeq_object_Fungi_roots)) # Put sample_data into a ggplot-friendly data.frame
df_roots$LibrarySize_roots <- sample_sums(physeq_object_Fungi_roots)
df_roots <- df_roots[order(df_roots$LibrarySize_roots),]
df_roots$Index <- seq(nrow(df_roots))
write.csv(df_roots, file = "rank_sums_roots.csv")
ggplot(data=df_roots, aes(x=Index, y=LibrarySize_roots, color=Sample_or_Control)) + geom_point()

# filter by prevelance 
sample_data(physeq_object_Fungi_roots)$is.neg <- sample_data(physeq_object_Fungi_roots)$Sample_or_Control == "Control Sample"
contamdf.prev_roots <- isContaminant(physeq_object_Fungi_roots, method="prevalence", neg="is.neg")
table(contamdf.prev_roots$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_roots <- transform_sample_counts(physeq_object_Fungi_roots, function(abund) 1*(abund>0))
ps.pa.neg_roots <- prune_samples(sample_data(ps.pa_roots)$Sample_or_Control == "Control Sample", ps.pa_roots)
ps.pa.pos_roots <- prune_samples(sample_data(ps.pa_roots)$Sample_or_Control == "True Sample", ps.pa_roots)
# Make data.frame of prevalence in positive and negative samples
df.pa_roots <- data.frame(pa.pos_roots=taxa_sums(ps.pa.pos_roots), pa.neg_roots=taxa_sums(ps.pa.neg_roots),
                          contaminant=contamdf.prev_roots$contaminant)
ggplot(data=df.pa_roots, aes(x=pa.neg_roots, y=pa.pos_roots, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_fungi_roots <- prune_taxa(!contamdf.prev_roots$contaminant, physeq_object_Fungi_roots)
# with contaminants removed
otu_table(ps.noncontam_fungi_roots)

#stems

# check library size distribution
df_stems <- as.data.frame(sample_data(physeq_object_Fungi_stems)) # Put sample_data into a ggplot-friendly data.frame
df_stems$LibrarySize_stems <- sample_sums(physeq_object_Fungi_stems)
df_stems <- df_stems[order(df_stems$LibrarySize_stems),]
df_stems$Index <- seq(nrow(df_stems))
write.csv(df_stems, file = "rank_sums_stems.csv")
ggplot(data=df_stems, aes(x=Index, y=LibrarySize_stems, color=Sample_or_Control)) + geom_point()

# filter by prevelance 
sample_data(physeq_object_Fungi_stems)$is.neg <- sample_data(physeq_object_Fungi_stems)$Sample_or_Control == "Control Sample"
contamdf.prev_stems <- isContaminant(physeq_object_Fungi_stems, method="prevalence", neg="is.neg")
table(contamdf.prev_stems$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_stems <- transform_sample_counts(physeq_object_Fungi_stems, function(abund) 1*(abund>0))
ps.pa.neg_stems <- prune_samples(sample_data(ps.pa_stems)$Sample_or_Control == "Control Sample", ps.pa_stems)
ps.pa.pos_stems <- prune_samples(sample_data(ps.pa_stems)$Sample_or_Control == "True Sample", ps.pa_stems)
# Make data.frame of prevalence in positive and negative samples
df.pa_stems <- data.frame(pa.pos_stems=taxa_sums(ps.pa.pos_stems), pa.neg_stems=taxa_sums(ps.pa.neg_stems),
                          contaminant=contamdf.prev_stems$contaminant)
ggplot(data=df.pa_stems, aes(x=pa.neg_stems, y=pa.pos_stems, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_fungi_stems <- prune_taxa(!contamdf.prev_stems$contaminant, physeq_object_Fungi_stems)
ps.noncontam_fungi_stems# with contaminants removed
otu_table(ps.noncontam_fungi_stems)

# leaves

# check library size distribution
df_leaves <- as.data.frame(sample_data(physeq_object_Fungi_leaves)) # Put sample_data into a ggplot-friendly data.frame
df_leaves$LibrarySize_leaves <- sample_sums(physeq_object_Fungi_leaves)
df_leaves <- df_leaves[order(df_leaves$LibrarySize_leaves),]
df_leaves$Index <- seq(nrow(df_leaves))
write.csv(df_leaves, file = "rank_sums_leaves.csv")
ggplot(data=df_leaves, aes(x=Index, y=LibrarySize_leaves, color=Sample_or_Control)) + geom_point()

# filter by prevelance 
sample_data(physeq_object_Fungi_leaves)$is.neg <- sample_data(physeq_object_Fungi_leaves)$Sample_or_Control == "Control Sample"
contamdf.prev_leaves <- isContaminant(physeq_object_Fungi_leaves, method="prevalence", neg="is.neg")
table(contamdf.prev_leaves$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_leaves <- transform_sample_counts(physeq_object_Fungi_leaves, function(abund) 1*(abund>0))
ps.pa.neg_leaves <- prune_samples(sample_data(ps.pa_leaves)$Sample_or_Control == "Control Sample", ps.pa_leaves)
ps.pa.pos_leaves <- prune_samples(sample_data(ps.pa_leaves)$Sample_or_Control == "True Sample", ps.pa_leaves)
# Make data.frame of prevalence in positive and negative samples
df.pa_leaves <- data.frame(pa.pos_leaves=taxa_sums(ps.pa.pos_leaves), pa.neg_leaves=taxa_sums(ps.pa.neg_leaves),
                           contaminant=contamdf.prev_leaves$contaminant)
ggplot(data=df.pa_leaves, aes(x=pa.neg_leaves, y=pa.pos_leaves, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_fungi_leaves <- prune_taxa(!contamdf.prev_leaves$contaminant, physeq_object_Fungi_leaves)
ps.noncontam_fungi_leaves# with contaminants removed
otu_table(ps.noncontam_fungi_leaves)

# remove negative controls from soil
ps.noncontam_fungi_soil <- subset_samples(ps.noncontam_fungi_soil, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_fungi_soil) <- otu_table(ps.noncontam_fungi_soil)[which(rowSums(otu_table(ps.noncontam_fungi_soil)) >= 1),]
ps.noncontam_fungi_soil
otu_table(ps.noncontam_fungi_soil)
sample_data(ps.noncontam_fungi_soil)
ps.noncontam_fungi_soil

# export soil otu table to check samples
# removing samples with less than 1000 reads, including samples that had less than 1000 read following sampling
write.csv(otu_table(ps.noncontam_fungi_soil), file = "filtering_low_soil.csv")


otu_table(ps.noncontam_fungi_soil) <- subset(otu_table(ps.noncontam_fungi_soil),
                                             select = -c(T4R1CR6S,T1R1BR6S,T4R2CR6S,T4R6AR2S,T1R1FBR3S,T1R2CCR3S,T1R1AR6S,T1R1FCR3S,T2R2FBR4S,T2R2FCR4S,T1R5CR2S,T1R6CCR4S,T2R5CR2S))


# export roots otu table to check samples

# remove negative controls
ps.noncontam_fungi_roots <- subset_samples(ps.noncontam_fungi_roots, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_fungi_roots) <- otu_table(ps.noncontam_fungi_roots)[which(rowSums(otu_table(ps.noncontam_fungi_roots)) >= 1),]
ps.noncontam_fungi_roots



write.csv(otu_table(ps.noncontam_fungi_roots), file = "filtering_low_roots.csv")


otu_table(ps.noncontam_fungi_roots) <- subset(otu_table(ps.noncontam_fungi_roots),
                                              select = -c(T1R5CBR3R,T2R1CBR6R,T4R5CR2R,T4R2BR2R,T4R1AR2R,T1R6FCR4R,T2R5AR2R,T1R1BR6R,T1R2FCR4R,T1R2AR2R,T1R2FAR4R,T1R6CCR4R,T1R2AR2R,T1R5FAR3R,T1R1FAR3R,T1R6CR2R,T1R2FAR3R,T4R6BR2R))
ps.noncontam_fungi_roots


# leaves
# remove negative controls
ps.noncontam_fungi_leaves <- subset_samples(ps.noncontam_fungi_leaves, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_fungi_leaves) <- otu_table(ps.noncontam_fungi_leaves)[which(rowSums(otu_table(ps.noncontam_fungi_leaves)) >= 1),]
ps.noncontam_fungi_leaves


# export leaves otu table to check samples
write.csv(otu_table(ps.noncontam_fungi_leaves), file = "filtering_low_leaves.csv")


otu_table(ps.noncontam_fungi_leaves) <- subset(otu_table(ps.noncontam_fungi_leaves),
                                               select = -c(T1R5AR2L,T2R1FCR6L))
sample_data(ps.noncontam_fungi_leaves)
#stems
ps.noncontam_fungi_stems <- subset_samples(ps.noncontam_fungi_stems, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_fungi_stems) <- otu_table(ps.noncontam_fungi_stems)[which(rowSums(otu_table(ps.noncontam_fungi_stems)) >= 1),]
ps.noncontam_fungi_stems



# export stems otu table to check samples
write.csv(otu_table(ps.noncontam_fungi_stems), file = "filtering_low_stems.csv")


otu_table(ps.noncontam_fungi_stems) <- subset(otu_table(ps.noncontam_fungi_stems),
                                              select = -c(T4R1CR2ST,T4R6CR2ST,T1R1AR4ST,T1R6FAR4ST,T2R5FCR3ST,T2R6AR2ST,T2R2CBR6ST,T1R6CR2ST,T1R6FBR6ST,T4R5AR2ST,T1R2AR2ST,T4R6CR6ST,T4R5CR2ST,T1R6FCR3ST,T2R5CR6ST,T2R6FBR3ST,T2R6FCR6ST,T4R5BR6ST,T1R2FBR4ST,T4R1BR6ST,T4R6CR2ST,T4R1AR6ST,T4R6AR6ST,T2R5CBR6ST))
sample_data(ps.noncontam_fungi_stems)

otu_table(ps.noncontam_fungi_roots) <- subset(otu_table(ps.noncontam_fungi_roots),
                                              select = -c(T4R2BV2R, T1R2AV2R))


otu_table(ps.noncontam_fungi_stems) <- subset(otu_table(ps.noncontam_fungi_stems),
                                              select = -c(T4R5CV2ST,T4R5AV2ST,T4R6CV2ST,T4R1CV2ST,T2R1BV2ST))
# extract obj 1
ps.noncontam_fungi_soil <- subset_samples(ps.noncontam_fungi_soil, experiment%in%c("obj_1"))
ps.noncontam_fungi_soil
ps.noncontam_fungi_roots <- subset_samples(ps.noncontam_fungi_roots, experiment%in%c("obj_1"))
ps.noncontam_fungi_roots
ps.noncontam_fungi_stems <- subset_samples(ps.noncontam_fungi_stems, experiment%in%c("obj_1"))
ps.noncontam_fungi_stems
ps.noncontam_fungi_leaves <- subset_samples(ps.noncontam_fungi_leaves, experiment%in%c("obj_1"))
ps.noncontam_fungi_leaves

# recombine fungi
ps.noncontam_fungi_total= merge_phyloseq(ps.noncontam_fungi_soil, ps.noncontam_fungi_roots, ps.noncontam_fungi_stems,ps.noncontam_fungi_leaves)
# add f in front of fungi ordersin order to be able to be able to merge prokaryotes and fungi and still 
# distinguish between the two


write.csv(otu_table(ps.noncontam_fungi_total), file = "otu_check.csv")
write.csv(tax_table(ps.noncontam_fungi_total), file = "tax_check.csv")
ps.noncontam_fungi_total



#prokaryotes
# create a phyloseq object ----------
prok_otus<- read.delim("otu_bac_total.txt",
                       row.names=1) 
head(prok_otus)


prok_otus_phy <-otu_table(prok_otus,
                          taxa_are_rows = TRUE)
prok_otus_phy


prok_metadata <-read.delim("total_prokaryote_map.txt",
                           row.names=1)
prok_metadata
prok_metadata_phy <-sample_data(prok_metadata)

prok_taxonomy<-read.delim("total_bac_tax.txt",
                          header=TRUE, 
                          row.names=1)
prok_taxonomy
prok_taxonomy_phy <- tax_table(as.matrix(prok_taxonomy))


prok_sequences <- readDNAStringSet("total_otus_bac.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
prok_sequences

physeq_prok_object <- phyloseq(prok_otus_phy, 
                               prok_metadata_phy, 
                               prok_taxonomy_phy,
                               prok_sequences)
physeq_prok_object
write.csv(tax_table(physeq_prok_object), file = "tax_check.csv")

#removing samples from soil dataset which were noted to not dry properly, and had mold growing in envelope. were analyzed, but only produced a few otus
#non organic samples on this list are being switched with the matching control sample from obj 2 
otu_table(physeq_prok_object) <- subset(otu_table(physeq_prok_object),
                                        select = -c(T4R5AR2S,T4R6BR2S,T4R2CR2S,T4R5CR2S,T2R6BR2S))

###formatting taxonomy------------------------
colnames(tax_table(physeq_prok_object)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
tax_table(physeq_prok_object)

tax_table(physeq_prok_object)[, "Kingdom"] <- gsub("k:", "", tax_table(physeq_prok_object)[, "Kingdom"])
tax_table(physeq_prok_object)[, "Phylum"] <- gsub("p:", "", tax_table(physeq_prok_object)[, "Phylum"])
tax_table(physeq_prok_object)[, "Class"] <- gsub("c:", "", tax_table(physeq_prok_object)[, "Class"])
tax_table(physeq_prok_object)[, "Order"] <- gsub("o:", "", tax_table(physeq_prok_object)[, "Order"])
tax_table(physeq_prok_object)[, "Family"] <- gsub("f:", "", tax_table(physeq_prok_object)[, "Family"])
tax_table(physeq_prok_object)[, "Genus"] <- gsub("g:", "", tax_table(physeq_prok_object)[, "Genus"])
tax_table(physeq_prok_object)[, "Species"] <- gsub("s:", "", tax_table(physeq_prok_object)[, "Species"])
tax_table(physeq_prok_object)


#remove chloropast,mitochondria,cyanobacteria
physeq_prok_object <- subset_taxa(physeq_prok_object, Phylum!="Chloroplast")
physeq_prok_object <- subset_taxa(physeq_prok_object, Class!="Chloroplast")
physeq_prok_object <- subset_taxa(physeq_prok_object, Order!="Chloroplast")
physeq_prok_object <- subset_taxa(physeq_prok_object, Family!="Chloroplast")
physeq_prok_object <- subset_taxa(physeq_prok_object, Genus!="Chloroplast")
tax_table(physeq_prok_object)

physeq_prok_object <- subset_taxa(physeq_prok_object, Phylum!="Mitochondria")
physeq_prok_object <- subset_taxa(physeq_prok_object, Class!="Mitochondria")
physeq_prok_object <- subset_taxa(physeq_prok_object, Order!="Mitochondria")
physeq_prok_object <- subset_taxa(physeq_prok_object, Family!="Mitochondria")
physeq_prok_object <- subset_taxa(physeq_prok_object, Genus!="Mitochondria")
tax_table(physeq_prok_object)

###remove contaminants--------------------------------------------------
# must split by sample source so that contaminants are removed from each individual run
physeq_prok_object_stems <- subset_samples(physeq_prok_object, origin%in%c("stem"))
physeq_prok_object_leaves <- subset_samples(physeq_prok_object,origin%in%c("leaves"))
physeq_prok_object_roots <- subset_samples(physeq_prok_object,origin%in%c("root"))
physeq_prok_object_soil <- subset_samples(physeq_prok_object,origin%in%c("soil"))
sample_data(physeq_prok_object_soil)
library(devtools)
library(processx)
devtools::install_github("benjjneb/decontam", force = TRUE)
library(decontam)
sample_data(physeq_prok_object_leaves)
write.csv(sample_data(physeq_prok_object_leaves),file ="sample_check.csv")

#leaves

# check library size distribution
df_leaves <- as.data.frame(sample_data(physeq_prok_object_leaves)) # Put sample_data into a ggplot-friendly data.frame
df_leaves$LibrarySize_leaves <- sample_sums(physeq_prok_object_leaves)
df_leaves <- df_leaves[order(df_leaves$LibrarySize_leaves),]
df_leaves$Index <- seq(nrow(df_leaves))
write.csv(df_leaves, file = "rank_sums_leaves.csv")
ggplot(data=df_leaves, aes(x=Index, y=LibrarySize_leaves, color=Sample_or_Control)) + geom_point()

# filter by prevelance 
sample_data(physeq_prok_object_leaves)$is.neg <- sample_data(physeq_prok_object_leaves)$Sample_or_Control == "Control Sample"
contamdf.prev_leaves <- isContaminant(physeq_prok_object_leaves, method="prevalence", neg="is.neg")
table(contamdf.prev_leaves$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_leaves <- transform_sample_counts(physeq_prok_object_leaves, function(abund) 1*(abund>0))
ps.pa.neg_leaves <- prune_samples(sample_data(ps.pa_leaves)$Sample_or_Control == "Control Sample", ps.pa_leaves)
ps.pa.pos_leaves <- prune_samples(sample_data(ps.pa_leaves)$Sample_or_Control == "True Sample", ps.pa_leaves)
# Make data.frame of prevalence in positive and negative samples
df.pa_leaves <- data.frame(pa.pos_leaves=taxa_sums(ps.pa.pos_leaves), pa.neg_leaves=taxa_sums(ps.pa.neg_leaves),
                           contaminant=contamdf.prev_leaves$contaminant)
ggplot(data=df.pa_leaves, aes(x=pa.neg_leaves, y=pa.pos_leaves, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_prok_leaves <- prune_taxa(!contamdf.prev_leaves$contaminant, physeq_prok_object_leaves)
ps.noncontam_prok_leaves # with contaminants removed
otu_table(ps.noncontam_prok_leaves)

#roots


# check library size distribution
df_roots <- as.data.frame(sample_data(physeq_prok_object_roots)) # Put sample_data into a ggplot-friendly data.frame
df_roots$LibrarySize_roots <- sample_sums(physeq_prok_object_roots)
df_roots <- df_roots[order(df_roots$LibrarySize_roots),]
df_roots$Index <- seq(nrow(df_roots))
write.csv(df_roots, file = "rank_sums_roots.csv")
ggplot(data=df_roots, aes(x=Index, y=LibrarySize_roots, color=Sample_or_Control)) + geom_point()

# filter by prevelance 
sample_data(physeq_prok_object_roots)$is.neg <- sample_data(physeq_prok_object_roots)$Sample_or_Control == "Control Sample"
contamdf.prev_roots <- isContaminant(physeq_prok_object_roots, method="prevalence", neg="is.neg")
table(contamdf.prev_roots$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_roots <- transform_sample_counts(physeq_prok_object_roots, function(abund) 1*(abund>0))
ps.pa.neg_roots <- prune_samples(sample_data(ps.pa_roots)$Sample_or_Control == "Control Sample", ps.pa_roots)
ps.pa.pos_roots <- prune_samples(sample_data(ps.pa_roots)$Sample_or_Control == "True Sample", ps.pa_roots)
# Make data.frame of prevalence in positive and negative samples
df.pa_roots <- data.frame(pa.pos_roots=taxa_sums(ps.pa.pos_roots), pa.neg_roots=taxa_sums(ps.pa.neg_roots),
                          contaminant=contamdf.prev_roots$contaminant)
ggplot(data=df.pa_roots, aes(x=pa.neg_roots, y=pa.pos_roots, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_prok_roots <- prune_taxa(!contamdf.prev_roots$contaminant, physeq_prok_object_roots)
# with contaminants removed
otu_table(ps.noncontam_prok_roots)

#stems

# check library size distribution
df_stems <- as.data.frame(sample_data(physeq_prok_object_stems)) # Put sample_data into a ggplot-friendly data.frame
df_stems$LibrarySize_stems <- sample_sums(physeq_prok_object_stems)
df_stems <- df_stems[order(df_stems$LibrarySize_stems),]
df_stems$Index <- seq(nrow(df_stems))
write.csv(df_stems, file = "rank_sums_stems.csv")
ggplot(data=df_stems, aes(x=Index, y=LibrarySize_stems, color=Sample_or_Control)) + geom_point()

# filter by prevelance 
sample_data(physeq_prok_object_stems)$is.neg <- sample_data(physeq_prok_object_stems)$Sample_or_Control == "Control Sample"
contamdf.prev_stems <- isContaminant(physeq_prok_object_stems, method="prevalence", neg="is.neg")
table(contamdf.prev_stems$contaminant)




# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_stems <- transform_sample_counts(physeq_prok_object_stems, function(abund) 1*(abund>0))
ps.pa.neg_stems <- prune_samples(sample_data(ps.pa_stems)$Sample_or_Control == "Control Sample", ps.pa_stems)
ps.pa.pos_stems <- prune_samples(sample_data(ps.pa_stems)$Sample_or_Control == "True Sample", ps.pa_stems)
# Make data.frame of prevalence in positive and negative samples
df.pa_stems <- data.frame(pa.pos_stems=taxa_sums(ps.pa.pos_stems), pa.neg_stems=taxa_sums(ps.pa.neg_stems),
                          contaminant=contamdf.prev_stems$contaminant)
ggplot(data=df.pa_stems, aes(x=pa.neg_stems, y=pa.pos_stems, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_prok_stems <- prune_taxa(!contamdf.prev_stems$contaminant, physeq_prok_object_stems)
ps.noncontam_prok_stems# with contaminants removed
otu_table(ps.noncontam_prok_stems)

# soil

# check library size distribution
df_soil <- as.data.frame(sample_data(physeq_prok_object_soil)) # Put sample_data into a ggplot-friendly data.frame
df_soil$LibrarySize_soil <- sample_sums(physeq_prok_object_soil)
df_soil <- df_soil[order(df_soil$LibrarySize_soil),]
df_soil$Index <- seq(nrow(df_soil))
write.csv(df_soil, file = "rank_sums_soil.csv")
ggplot(data=df_soil, aes(x=Index, y=LibrarySize_soil, color=Sample_or_Control)) + geom_point()

# filter by prevelance 
sample_data(physeq_prok_object_soil)$is.neg <- sample_data(physeq_prok_object_soil)$Sample_or_Control == "Control Sample"
contamdf.prev_soil <- isContaminant(physeq_prok_object_soil, method="prevalence", neg="is.neg")
table(contamdf.prev_soil$contaminant)



# Make phyloseq object of presence-absence in negative controls and true samples
ps.pa_soil <- transform_sample_counts(physeq_prok_object_soil, function(abund) 1*(abund>0))
ps.pa.neg_soil <- prune_samples(sample_data(ps.pa_soil)$Sample_or_Control == "Control Sample", ps.pa_soil)
ps.pa.pos_soil <- prune_samples(sample_data(ps.pa_soil)$Sample_or_Control == "True Sample", ps.pa_soil)
# Make data.frame of prevalence in positive and negative samples
df.pa_soil <- data.frame(pa.pos_soil=taxa_sums(ps.pa.pos_soil), pa.neg_soil=taxa_sums(ps.pa.neg_soil),
                         contaminant=contamdf.prev_soil$contaminant)
ggplot(data=df.pa_soil, aes(x=pa.neg_soil, y=pa.pos_soil, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)")

# remove contaminants
ps.noncontam_prok_soil <- prune_taxa(!contamdf.prev_soil$contaminant, physeq_prok_object_soil)
ps.noncontam_prok_soil# with contaminants removed
otu_table(ps.noncontam_prok_soil)
sample_data(ps.noncontam_prok_soil)
# remove negative controls from soil
ps.noncontam_prok_soil <- subset_samples(ps.noncontam_prok_soil, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_prok_soil) <- otu_table(ps.noncontam_prok_soil)[which(rowSums(otu_table(ps.noncontam_prok_soil)) >= 1),]
ps.noncontam_prok_soil
otu_table(ps.noncontam_prok_soil)
sample_data(ps.noncontam_prok_soil)
ps.noncontam_prok_soil


# removing samples with less than 1000 reads, including samples that had less than 1000 read following sampling
write.csv(otu_table(ps.noncontam_prok_soil), file = "filtering_low_soil.csv")


otu_table(ps.noncontam_prok_soil) <- subset(otu_table(ps.noncontam_prok_soil),
                                            select = -c(T1R6CCR3S, T1R6CBR3S, T1R2CBR3S, T2R1CCR3S, T1R5CR2S, T4R5CR2S, T4R5CR2S, T2R6BR2S,T1R5CR2S))


# export roots otu table to check samples
write.csv(otu_table(ps.noncontam_prok_roots), file = "filtering_low_roots.csv")


otu_table(ps.noncontam_prok_roots) <- subset(otu_table(ps.noncontam_prok_roots),
                                             select = -c(T1R5FAR3R, T1R6CR2R, NEG11R, NEG8R, NEG2R, T1R2AV2R, T1R5CAR6R, T1R6CCR4R, T2R1CAR6R, T1R5FCR6R, T1R6CBR6R))

# export leaves otu table to check samples
write.csv(otu_table(ps.noncontam_prok_leaves), file = "filtering_low_leaves.csv")


otu_table(ps.noncontam_prok_leaves) <- subset(otu_table(ps.noncontam_prok_leaves),
                                              select = -c(T2R6FBR3L,NEG17L,T4R1CV2L,T1R1BV2L,NEG9L,NEG16L,Neg14L,T4R1BR6L,T1R5AR2L,NEG10L,NEG12L,NEG5L,T4R5BV2L,NEG15L,Neg2L,NEG3L,Neg8L,NEG13L,NEG11L,T1R1BR2L, T4R2BV2L, T4R1AV2L))

# export stems otu table to check samples
write.csv(otu_table(ps.noncontam_prok_stems), file = "filtering_low_stems.csv")


otu_table(ps.noncontam_prok_stems) <- subset(otu_table(ps.noncontam_prok_stems),
                                             select = -c(NEG10ST,T4R5BR6ST,NEG13ST,T4R1CV2ST,NEG6ST,T4R6AR6ST,T1R1CBR3ST,NEG15ST,NEG14ST,NEG5ST,T1R2AR6ST,T1R2BR6ST,T4R1AR6ST,NEG2ST,NEG3ST,T4R6CR2ST,NEG9ST,NEG11ST,T2R6FCR6ST,T4R1BR6ST,T1R2FBR4ST,Mock2ST, T1R5AV2ST))
otu_table(ps.noncontam_prok_soil) <- subset(otu_table(ps.noncontam_prok_soil),
                                            select = -c(T1R5CR2S))

# extract obj 1
ps.noncontam_prok_soil <- subset_samples(ps.noncontam_prok_soil, experiment%in%c("obj_1"))
ps.noncontam_prok_soil
ps.noncontam_prok_roots <- subset_samples(ps.noncontam_prok_roots, experiment%in%c("obj_1"))
ps.noncontam_prok_roots
ps.noncontam_prok_stems <- subset_samples(ps.noncontam_prok_stems, experiment%in%c("obj_1"))
ps.noncontam_prok_stems
ps.noncontam_prok_leaves <- subset_samples(ps.noncontam_prok_leaves, experiment%in%c("obj_1"))
ps.noncontam_prok_leaves

# recombine prok
ps.noncontam_prok_total= merge_phyloseq(ps.noncontam_prok_soil, ps.noncontam_prok_roots, ps.noncontam_prok_stems,ps.noncontam_prok_leaves)

ps.noncontam_prok_total
# Paste B in front so that bacterial and fungal reads can be distinguished from each other
otu_names_16s<- rownames(otu_table(ps.noncontam_prok_total))
otu_names_16s <- paste("B", otu_names_16s, sep = "")
taxa_names(ps.noncontam_prok_total) <- paste(otu_names_16s)


# adding the B in front of the OTU name allows for distinction from fungi
phyloseq_total <- merge_phyloseq(ps.noncontam_prok_total,ps.noncontam_fungi_total)


# ****** EXTRACTING conventional-----------------------------------------
sample_data(ps.noncontam_fungi_total)
ps.noncontam_fungi_total_conventional <- subset_samples(ps.noncontam_fungi_total, Management%in%c("Conventional"))
otu_table(ps.noncontam_fungi_total_conventional) <- otu_table(ps.noncontam_fungi_total_conventional)[which(rowSums(otu_table(ps.noncontam_fungi_total_conventional)) >= 1),] 




ps.noncontam_prok_total_conventional <- subset_samples(ps.noncontam_prok_total, Management_Indicator%in%c("Conventional"))
otu_table(ps.noncontam_prok_total_conventional) <- otu_table(ps.noncontam_prok_total_conventional)[which(rowSums(otu_table(ps.noncontam_prok_total_conventional)) >= 1),] 

network_prok_below_conventional <- subset_samples(ps.noncontam_prok_total_conventional, origin%in%c("soil","root"))
network_prok_below_conventional
network_fungi_below_conventional <- subset_samples(ps.noncontam_fungi_total_conventional, origin%in%c("soil","root"))
network_fungi_below_conventional
### below conventional---------------------------------
# extracting the core taxa
#try defining as 80%
network_prok_below_conventional
network_prok_below_conventional <- filter_taxa(network_prok_below_conventional, function(x) sum(x >= 1) >= (57), TRUE)
sample_names(network_prok_below_conventional) <- paste("sample", 1:nsamples(network_prok_below_conventional), sep="")
sample_data(network_prok_below_conventional)
network_prok_below_conventional

network_fungi_below_conventional <- filter_taxa(network_fungi_below_conventional, function(x) sum(x >= 1) >= (48), TRUE)
sample_names(network_fungi_below_conventional) <- paste("sample", 1:nsamples(network_fungi_below_conventional), sep="")
network_fungi_below_conventional

# merging the dataset  combines datasets to make the network-----------------------
otu_names_ITS_conventional_below <- rownames(otu_table(network_fungi_below_conventional))
otu_names_ITS_conventional_below <- paste("F", otu_names_ITS_conventional_below, sep="")
taxa_names(network_fungi_below_conventional) <- paste(otu_names_ITS_conventional_below)
head(otu_table(network_fungi_below_conventional))
head(tax_table(network_fungi_below_conventional))
head(refseq(network_fungi_below_conventional))
sample_data(network_fungi_below_conventional)

otu_names_16s_conventional_below <- rownames(otu_table(network_prok_below_conventional))
otu_names_16s_conventional_below <- paste("B", otu_names_16s_conventional_below, sep = "")
taxa_names(network_prok_below_conventional) <- paste(otu_names_16s_conventional_below)
head(otu_table(network_prok_below_conventional))
head(tax_table(network_prok_below_conventional))
head(refseq(network_prok_below_conventional))
sample_data(network_prok_below_conventional)

# Making total merged otu table for plotting hub barplot

otu_names_ITS <- rownames(otu_table(ps.noncontam_fungi_total))
otu_names_ITS <- paste("F", otu_names_ITS_conventional_below, sep="")
otu_names_ITS
taxa_names(ps.noncontam_fungi_total) <- paste(otu_names_ITS)


otu_names_16s_conventional_below <- rownames(otu_table(network_prok_below_conventional))
otu_names_16s_conventional_below <- paste("B", otu_names_16s_conventional_below, sep = "")
taxa_names(network_prok_below_conventional) <- paste(otu_names_16s_conventional_below)
head(otu_table(network_prok_below_conventional))
head(tax_table(network_prok_below_conventional))
head(refseq(network_prok_below_conventional))
sample_data(network_prok_below_conventional)

conventional_merged_below <- merge_phyloseq(network_prok_below_conventional, network_fungi_below_conventional)
conventional_merged_below

write.csv(otu_table(conventional_merged_below), file = "conventional_net_below.csv")
write.csv(sample_data(conventional_merged_below),file="conventional_net_below_map.csv")

merged_total <- merge_phyloseq(ps.noncontam_prok_total,ps.noncontam_fungi_total)

# >>> CREATING A SPIECEASI NETWORK -------------
# https://matthewlincoln.net/2014/12/20/adjacency-matrix-plots-with-r-and-ggplot2.html

library("SpiecEasi"); packageVersion("SpiecEasi")
library("igraph"); packageVersion("igraph")
library("MASS")
library("huge"); packageVersion("huge")
library("qgraph"); packageVersion("qgraph")
library("MCL")
memory.limit()
memory.limit(size=50000)
# this code creates the spieceasi network
spiec_funpro_conventional_below <- spiec.easi(conventional_merged_below, 
                                        method="mb",
                                        lambda.min.ratio=1e-2, # the higher the more edges, higher density (5e-1 to 1e-1, 1e-2, 1e-3)
                                        nlambda=50, # number of slices between minimum amount od edged to maximum
                                        sel.criterion ="stars", # selects the lambda and modules the most stable connections
                                        pulsar.select = TRUE,
                                        pulsar.params=list(rep.num=99)) #iteration of starts, number of times sample the 80% of the data within each lambda

spiec_funpro_conventional_below
str(spiec_funpro_conventional_below)
spiec_funpro_conventional_below$lambda
spiec_funpro_conventional_below$select$stars$summary
getOptLambda(spiec_funpro_conventional_below)
getOptInd(spiec_funpro_conventional_below)
getStability(spiec_funpro_conventional_below)
spiec_funpro_conventional_below


# adding weights to the graph is you like 
plot_spiec_funpro_conventional_below <- adj2igraph(Matrix::drop0(getRefit(spiec_funpro_conventional_below)),
                                             vertex.attr = list(name=taxa_names(conventional_merged_below)),
                                             rmEmptyNodes = FALSE)
plot_spiec_funpro_conventional_below


# Get clusters - identifies clusters which will later be labelled on the networks
wt_conventional_below <- walktrap.community(plot_spiec_funpro_conventional_below)
ml_conventional_below <- multilevel.community(plot_spiec_funpro_conventional_below)
wt_conventional_below
commSummary_conventional_below <- data.frame(
  wt_conventional_below$names,
  wt_conventional_below$membership,
  wt_conventional_below$modularity
)
write.csv(commSummary_conventional_below, file ="conv_below_modules.csv")
commSummary_conventional_below
options("max.print")


# Hub detection - the "hubscore" is another way to identify network hubs, but was not used in this network



net.cs_conventional_below <- closeness(plot_spiec_funpro_conventional_below)

net.bn_conventional_below <- betweenness(plot_spiec_funpro_conventional_below) 
net.pr_conventional_below <- page_rank(plot_spiec_funpro_conventional_below)$vector 
net.hs_conventional_below <- hub_score(plot_spiec_funpro_conventional_below)$vector
net.dg_conventional_below <- degree(plot_spiec_funpro_conventional_below)

hubs_conventional_below <- data.frame(hub_score(plot_spiec_funpro_conventional_below)$vector); colnames(hubs_conventional_below) <- c("hubscore")
hubs_conventional_below$OTU <- rownames(hubs_conventional_below)
hubs_conventional_below$betweeness <- net.bn_conventional_below
hubs_conventional_below$pr <- net.pr_conventional_below
hubs_conventional_below$degree <- net.dg_conventional_below
hubs_conventional_below$closness <- net.cs_conventional_below
# check if normal, if not log transform
hist(hubs_conventional_below$closness) # fine
hist(hubs_conventional_below$betweeness) # transform
hist(hubs_conventional_below$degree) # fine
library(dplyr)
library(tidyr)
# getting the mean relative abundance of each OTU
conventional_merged.ra_below <- conventional_merged_below %>%   
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%    # Melt to long format
  arrange(OTU)

mean.rel.abund_conventional_below <- conventional_merged.ra_below %>% 
  group_by(OTU) %>% 
  nest() %>%
  mutate(mean.relabund_conventional_below = purrr::map(data,~mean(.$Abundance*100))) %>%
  mutate(SE.relabund = purrr::map(data,~sd(.$Abundance*100)/sqrt(length(.$Abundance*100)))) %>%
  unnest(mean.relabund_conventional_below, SE.relabund)

# hub analysis - see methods section for detail
conventional_tax_below <- data.frame(tax_table(conventional_merged_below))
nodes.stats_conventional_below <- hubs_conventional_below
nodes.stats_conventional_below$mean.abundance <- mean.rel.abund_conventional_below$mean.relabund_conventional_below[match(nodes.stats_conventional_below$OTU, mean.rel.abund_conventional_below$OTU)]
nodes.stats_conventional_below$Label <- conventional_tax_below$Label[match(nodes.stats_conventional_below$OTU, conventional_tax_below$OTU)]
nodes.stats_conventional_below$Genus <- conventional_merged.ra_below$Genus[match(nodes.stats_conventional_below$OTU, conventional_merged.ra_below$OTU)]
nodes.stats_conventional_below$Kingdom <- conventional_merged.ra_below$Kingdom[match(nodes.stats_conventional_below$OTU, conventional_merged.ra_below$OTU)]
nodes.stats_conventional_below$Phylum <- conventional_merged.ra_below$Phylum[match(nodes.stats_conventional_below$OTU, conventional_merged.ra_below$OTU)]
mean.close_conventional_below <- mean(nodes.stats_conventional_below$closness)
mean.close_conventional_below
sd.close_conventional_below <- sd(nodes.stats_conventional_below$closness)

# this creates the thresholds for calling taxa hubs
hubline.close_conventional_below <- (mean.close_conventional_below + 1.1*sd.close_conventional_below)
hubline.close_conventional_below
z.score.close_conventional_below = (hubline.close_conventional_below - mean.close_conventional_below)/sd.close_conventional_below
pnorm(z.score.close_conventional_below) # line is above 90 % - equal to p = 0.1
nodes.stats_conventional_below
mean.degree_conventional_below <- mean(nodes.stats_conventional_below$degree)
sd.degree_conventional_below <- sd(nodes.stats_conventional_below$degree)
hubline.degree_conventional_below <- (mean.degree_conventional_below + 1.3*sd.degree_conventional_below)

z.score.degree_conventional_below = (hubline.degree_conventional_below - mean.degree_conventional_below)/sd.degree_conventional_below
pnorm(z.score.degree_conventional_below) # line is above 90 % - equal to p = 0.1

mean.between_conventional_below <- mean(log10(nodes.stats_conventional_below$betweeness[nodes.stats_conventional_below$betweeness > 0]))
sd.between_conventional_below <- sd(log10(nodes.stats_conventional_below$betweeness[nodes.stats_conventional_below$betweeness > 0]))
hubline.between_conventional_below <- (mean.between_conventional_below + 1.3*sd.between_conventional_below)

z.score.between_conventional_below = (hubline.between_conventional_below - mean.between_conventional_below)/sd.between_conventional_below
pnorm(z.score.between_conventional_below) # line is above 90 % - equal to p = 0.1

library(ggrepel)
# plotting hub detections
close_conventional_below <- ggplot() + 
  geom_point(data = nodes.stats_conventional_below, aes(size = mean.abundance, x = closness, y = degree), alpha = 0.6) +
  scale_size_continuous(name = "Relative Abundance") +
  theme_bw() + 
  geom_text_repel(data = subset(nodes.stats_conventional_below, closness > hubline.close_conventional_below & degree > hubline.degree_conventional_below), aes(x = closness, y = degree, label = OTU)) +
  xlab("Closeness Centrality") + 
  ylab("Degree") + 
  geom_vline(xintercept = hubline.close_conventional_below, linetype = 2) + 
  geom_hline(yintercept = hubline.degree_conventional_below, linetype = 2)
close_conventional_below
between_conventional_below <- ggplot() + 
  geom_point(data = nodes.stats_conventional_below, aes(size = mean.abundance, x = betweeness, y = degree), alpha = 0.6) +
  scale_size_continuous(name = "Relative Abundance") +
  theme_bw() + 
  geom_text_repel(data = subset(nodes.stats_conventional_below, betweeness > 10^hubline.between_conventional_below & degree > hubline.degree_conventional_below), aes(x = betweeness, y = degree, label = OTU)) +
  xlab("Betweeness Centrality") + 
  ylab("Degree") +
  geom_vline(xintercept = 10^hubline.between_conventional_below, linetype = 2) + 
  geom_hline(yintercept = hubline.degree_conventional_below, linetype = 2)

between_conventional_below

# Considered to be a hub if above the thresholds for degree, betweeness centrality and closeness centrality


library(cowplot)
#combining hubplots
hub.plot <- plot_grid(close_conventional_below, between_conventional_below, labels = "AUTO", nrow = 2)
plot(hub.plot)
dev.off()
### add pos-neg to table
# this labels edges between nodes as either positive or negative correlations
betaMat_conventional_below=as.matrix(symBeta(getOptBeta(spiec_funpro_conventional_below)))
otu.ids_conventional_below=colnames(spiec_funpro_conventional_below$est$data)
edges_conventional_below=E(plot_spiec_funpro_conventional_below)
edge.stats_conventional_below <- NULL
for(e.index in 1:length(edges_conventional_below)){
  adj.nodes=ends(plot_spiec_funpro_conventional_below,edges_conventional_below[e.index])
  xindex=which(otu.ids_conventional_below==adj.nodes[1])
  yindex=which(otu.ids_conventional_below==adj.nodes[2])
  beta=betaMat_conventional_below[xindex,yindex]
  wt_i_conventional_below <- cbind.data.frame(adj.nodes, beta)
  colnames(wt_i_conventional_below) <- c("source", "target", "beta")
  edge.stats_conventional_below <- rbind.data.frame(wt_i_conventional_below, edge.stats_conventional_below)
}

edge.stats_conventional_below$direction <- ifelse(edge.stats_conventional_below$beta > 0, "Positive", "Negative")
edge.stats_conventional_below$shared.name <- paste(edge.stats_conventional_below$source, "(interacts with)", edge.stats_conventional_below$target)



library(RCy3)
# exporting 
#exporting files to create cytoscape project. Will also combine node stats file with modules files
createNetworkFromIgraph(plot_spiec_funpro_conventional_below,"conventional_net_below")

#Node stats - import this first.
write.csv(nodes.stats_conventional_below, "conventional_node_stats_below.csv")

# edge stats
write.csv(edge.stats_conventional_below, "conventional_edge_stats_below.csv")

# Visualization done in cytoscape as described in the methods section
# Other networks built in the same way, some used different core thresholds (see methods)

# making hub barplots figure 5--------------
# subset hubs only
tax_table(merged_total)
sample_data(phyloseq_total)

# OTU 99 - was not classified at the genus level so before startng making these barplots
# I went back to the original taxonomy tabel used to create the first phyloseq object
# and  labelled it Ascomycete_sp.
physeq.subset_hubs <- subset_taxa(phyloseq_total, rownames(tax_table(phyloseq_total)) %in% c("OTU_3563", "BOTU_57","BOTU_17","OTU_99","BOTU_58","BOTU_10038","BOTU_418","BOTU_194","BOTU_21","OTU_56","OTU_7","BOTU_22","BOTU_27","BOTU_87"))

tax_table(physeq.subset_hubs)





sample_data(physeq.subset_hubs)$Indicator_label <- factor(sample_data(physeq.subset_hubs)$Indicator_label,
                                                    levels = c("soil Conventional R2",	"soil Conventional R6",	"soil Conventional V2", "soil No-Till R2",	"soil No-Till R6",	"soil No-Till V2",	"soil Organic R2",	"soil Organic R6",	"soil Organic V2",	"root Conventional R2",	"root Conventional R6",	"root Conventional V2",	"root No-Till R2",	"root No-Till R6",	"root No-Till V2",	"root Organic R2",	"root Organic R6",	"root Organic V2",	"stem Conventional R2",	"stem Conventional R6",	"stem Conventional V2",	"stem No-Till R2",	"stem No-Till R6",	"stem No-Till V2",	"stem Organic R2",	"stem Organic R6",	"stem Organic V2",	"leaves Conventional R2",	"leaves Conventional R6",	"leaves Conventional V2",	"leaves No-Till R2",	"leaves No-Till R6",	"leaves No-Till V2",	"leaves Organic R2",	"leaves Organic R6",	"leaves Organic V2"))
write.csv(sample_data(hub_barplots), file ="hub_barplots.csv")

hub_barplots <- merge_samples(physeq.subset_hubs, "Indicator_label")
hub_barplots <- hub_barplots%>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by class
hub_barplots

library(dplyr)
library(RCy3)
library(data.table)
dat_hub_bp <- data.table(hub_barplots)
dat_hub_bp[(Abundance <= 0.00000001), Genus:= "Other"]
write.csv(dat_hub_bp, file = "hub_df.csv")
# create new phyloseq object

dat_hub_bp




dat_hub_bp
library(ggplot2)
dev.off()
# Plot
bar_hubs= ggplot(dat_hub_bp, aes(x = Indicator_label, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("soil Conventional R2",	"soil Conventional R6",	"soil Conventional V2", "soil No-Till R2",	"soil No-Till R6",	"soil No-Till V2",	"soil Organic R2",	"soil Organic R6",	"soil Organic V2",	"root Conventional R2",	"root Conventional R6",	"root Conventional V2",	"root No-Till R2",	"root No-Till R6",	"root No-Till V2",	"root Organic R2",	"root Organic R6",	"root Organic V2",	"stem Conventional R2",	"stem Conventional R6",	"stem Conventional V2",	"stem No-Till R2",	"stem No-Till R6",	"stem No-Till V2",	"stem Organic R2",	"stem Organic R6",	"stem Organic V2",	"leaves Conventional R2",	"leaves Conventional R6",	"leaves Conventional V2",	"leaves No-Till R2",	"leaves No-Till R6",	"leaves No-Till V2",	"leaves Organic R2",	"leaves Organic R6",	"leaves Organic V2"))+
  scale_fill_manual(values = c("Acidobacteria_RB41" ="#CBD588",
                               "Aquincola" = "#5F7FC7", 
                               "Ascomycete_sp." = "coral",
                               "Blastococcus" = "#DA5724",
                               "Deinococcus" = "#508578",
                               "Kineococcus" = "#CD9BCD",
                               "Marmoricola" = "#AD6F3B",
                               "Massilia" = "#673770", 
                               "Pelomonas" = "#D14285",
                               "Phoma" = "#652926",
                               "Podospora" = "#C84248",
                               "Sphingomonas" = "#8569D5",
                               "Streptomyces" = "steelblue",
                               "Other" = "yellow"
                               
                               
  ))+
  # Remove x axis title
  theme(axis.title.x = element_blank()) + 
  theme(legend.key.height = unit(0.15, "cm"), legend.key.width = unit(0.25, "cm")) +
  theme(legend.title = element_text(size = 8, face = "bold"), legend.text = element_text(size = 7)) +
  theme(strip.text.x = element_text(size = 8, face = "bold")) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  theme(plot.title = element_text(size = 14, hjust = 0.5)) +
  theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, hjust = 1)) +
  #theme(axis.ticks.x = element_blank()) +
  theme(axis.title = element_text(angle = 0, size = 8, face = "bold")) +
  theme(legend.position="right")+
  
  guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
  ylab("Relative Abundance of Hub Genera") 
plot(bar_hubs)


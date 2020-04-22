

# Network Creation Following Reviewer Comments


# >>> loading packages --------------------------------------------------------------
library("phyloseq"); packageVersion("phyloseq")
library("ggplot2")
library("Biostrings")
library("SpiecEasi"); packageVersion("SpiecEasi")
library("igraph"); packageVersion("igraph")

# fungi phyloseq object
ITS_otus<- read.delim("otu_table_Fungi.txt",
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


ITS_sequences <- readDNAStringSet("otus_total_Fungi.fasta", format="fasta", seek.first.rec=TRUE, use.names=TRUE)
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



#removing samples from roots dataset which were noted to not dry properly, and had mold growing in envelope. were analyzed, but only produced a few otus
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



# export roots otu table to check samples

# remove negative controls
ps.noncontam_fungi_roots <- subset_samples(ps.noncontam_fungi_roots, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_fungi_roots) <- otu_table(ps.noncontam_fungi_roots)[which(rowSums(otu_table(ps.noncontam_fungi_roots)) >= 1),]
ps.noncontam_fungi_roots






# leaves
# remove negative controls
ps.noncontam_fungi_leaves <- subset_samples(ps.noncontam_fungi_leaves, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_fungi_leaves) <- otu_table(ps.noncontam_fungi_leaves)[which(rowSums(otu_table(ps.noncontam_fungi_leaves)) >= 1),]
ps.noncontam_fungi_leaves



#stems
ps.noncontam_fungi_stems <- subset_samples(ps.noncontam_fungi_stems, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_fungi_stems) <- otu_table(ps.noncontam_fungi_stems)[which(rowSums(otu_table(ps.noncontam_fungi_stems)) >= 1),]
ps.noncontam_fungi_stems







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




# extract obj 1
ps.noncontam_prok_soil <- subset_samples(ps.noncontam_prok_soil, experiment%in%c("obj_1"))
ps.noncontam_prok_soil
ps.noncontam_fungi_soil
# must remove an extra sample from the prokaryotes so that samples exactly match the fungi
otu_table(ps.noncontam_prok_soil) <- subset(otu_table(ps.noncontam_prok_soil),
                                            select = -c(T4R1BR6S))



# recombine prok
ps.noncontam_prok_total= merge_phyloseq(ps.noncontam_prok_soil, ps.noncontam_prok_roots, ps.noncontam_prok_stems,ps.noncontam_prok_leaves)

ps.noncontam_prok_total

ps.noncontam_fungi_total
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
network_prok_below_conventional <- filter_taxa(network_prok_below_conventional, function(x) sum(x >= 1) >= (58), TRUE)
network_prok_below_conventional
sample_names(network_prok_below_conventional) <- paste("sample", 1:nsamples(network_prok_below_conventional), sep="")
sample_data(network_prok_below_conventional)
network_prok_below_conventional
network_fungi_below_conventional
network_fungi_below_conventional <- filter_taxa(network_fungi_below_conventional, function(x) sum(x >= 1) >= (58), TRUE)
network_fungi_below_conventional
sample_names(network_fungi_below_conventional) <- paste("sample", 1:nsamples(network_fungi_below_conventional), sep="")
network_fungi_below_conventional

# changing otu labels to distinguish between fungi and bacteria
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
otu_table(network_fungi_below_conventional)
conventional_merged_below <- merge_phyloseq(network_prok_below_conventional, network_fungi_below_conventional)
otu_table(conventional_merged_below)


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
network_prok_below_conventional
network_fungi_below_conventional
spiec_funpro_conventional_below <- SpiecEasi::spiec.easi(list(network_prok_below_conventional, network_fungi_below_conventional),
                                                         method="mb",
                                                         lambda.min.ratio=1e-1, # the higher the more edges, higher density (5e-1 to 1e-1, 1e-2, 1e-3)
                                                         nlambda=50, # number of slices between minimum amount od edged to maximum
                                                         sel.criterion ="stars", # selects the lambda and modules the most stable connections
                                                         pulsar.select = TRUE,
                                                         
                                                         pulsar.params=list(rep.num=99)) #iteration of starts, number of times sample the 80% of the data within each lambda

otu_table(network_prok_below_conventional)

spiec_funpro_conventional_below
str(spiec_funpro_conventional_below)
spiec_funpro_conventional_below$lambda
spiec_funpro_conventional_below$select$stars$summary
getOptLambda(spiec_funpro_conventional_below)
getOptInd(spiec_funpro_conventional_below)
getStability(spiec_funpro_conventional_below)
spiec_funpro_conventional_below


# Compare produced network to random network with same number of nodes
g_con_below <- sample_pa(n = 139,start.graph =con_below_igraph) # this makes the random graph (start graph is an igraph object with your actual graph)
 x <- degree_distribution(g_con_below) # this gets the degree distribution of the random graph 
g <- plot.igraph(g_con_below) #to see random plor 
y <- degree_distribution(plot_spiec_funpro_conventional_below) # get degree distribution of actual network
ks.test(x,y) # test for statistical differences between random and real network
x2 <- ks.test(x,y)
degamat_con_below <- NULL
n <- 100
for(i in 1:n){
  newmatrix <- sample_pa(n =139, start.graph =con_below_igraph)
  degmat_con_below <- degree_distribution(newmatrix)
  degamat_con_below<-rbind(degamat_con_below,degmat_con_below)
}

degamat_con_below

# test against 100 random networks
ks.test_con_below <- NULL
for (i in 1:nrow(degamat_con_below)){
  ks.test_new <- ks.test(degamat_con_below[i,],y)
  ks.test_con_below$D[i]       <- ks.test_new$statistic
  ks.test_con_below$p[i]       <- ks.test_new$p.value
}
ks.test_con_below
write.csv(ks.test_con_below, file ="conventional_below_p.csv")


# making graph
con_below_igraph <- plot.igraph(plot_spiec_funpro_conventional_below)

# adding weights to the graph is you like 
plot_spiec_funpro_conventional_below <- adj2igraph(Matrix::drop0(getRefit(spiec_funpro_conventional_below)),
                                                   vertex.attr = list(name=taxa_names(conventional_merged_below)),
                                                   rmEmptyNodes = FALSE)





# Get clusters
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
# Hub detection
net.cs_conventional_below <- closeness(plot_spiec_funpro_conventional_below)

net.bn_conventional_below <- betweenness(plot_spiec_funpro_conventional_below) 
net.pr_conventional_below <- page_rank(plot_spiec_funpro_conventional_below)$vector 
net.hs_conventional_below <- hub_score(plot_spiec_funpro_conventional_below)$vector
net.hs_conventional_below
net.dg_conventional_below <- degree(plot_spiec_funpro_conventional_below)
net.dg_conventional_below
hubs_conventional_below <- data.frame(hub_score(plot_spiec_funpro_conventional_below)$vector); colnames(hubs_conventional_below) <- c("hubscore")
hubs_conventional_below$OTU <- rownames(conventional_merged_below)
hubs_conventional_below$betweeness <- net.bn_conventional_below
hubs_conventional_below$pr <- net.pr_conventional_below
hubs_conventional_below$degree <- net.dg_conventional_below
hubs_conventional_below$closness <- net.cs_conventional_below
hubs_conventional_below
# check if normal, if not log transform
hist(hubs_conventional_below$closness)
hist(log10(hubs_conventional_below$hubscore))
hist(hubs_conventional_below$betweeness) # transform
hist(hubs_conventional_below$degree) 
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
  unnest(c(mean.relabund_conventional_below, SE.relabund))
conventional_tax_below <- data.frame(tax_table(conventional_merged_below))
nodes.stats_conventional_below <- hubs_conventional_below
nodes.stats_conventional_below$mean.abundance <- mean.rel.abund_conventional_below$mean.relabund_conventional_below[match(nodes.stats_conventional_below$OTU, mean.rel.abund_conventional_below$OTU)]
nodes.stats_conventional_below$Label <- conventional_tax_below$Label[match(nodes.stats_conventional_below$OTU, conventional_tax_below$OTU)]
nodes.stats_conventional_below$Genus <- conventional_merged.ra_below$Genus[match(nodes.stats_conventional_below$OTU, conventional_merged.ra_below$OTU)]
nodes.stats_conventional_below$Kingdom <- conventional_merged.ra_below$Kingdom[match(nodes.stats_conventional_below$OTU, conventional_merged.ra_below$OTU)]
nodes.stats_conventional_below$Phylum <- conventional_merged.ra_below$Phylum[match(nodes.stats_conventional_below$OTU, conventional_merged.ra_below$OTU)]
mean.close_conventional_below <- mean(nodes.stats_conventional_below$closness)
mean.close_conventional_below

nodes.stats_conventional_below
mean.degree_conventional_below <- mean(nodes.stats_conventional_below$degree)
mean.degree_conventional_below
sd.degree_conventional_below <- sd(nodes.stats_conventional_below$degree)
hubline.degree_conventional_below <- (mean.degree_conventional_below + 1.3*sd.degree_conventional_below)

z.score.degree_conventional_below = (hubline.degree_conventional_below - mean.degree_conventional_below)/sd.degree_conventional_below
pnorm(z.score.degree_conventional_below) # line is above 90 % - equal to p = 0.1

mean.between_conventional_below <- mean(log10(nodes.stats_conventional_below$betweeness[nodes.stats_conventional_below$betweeness > 0]))
sd.between_conventional_below <- sd(log10(nodes.stats_conventional_below$betweeness[nodes.stats_conventional_below$betweeness > 0]))
hubline.between_conventional_below <- (mean.between_conventional_below + 1.3*sd.between_conventional_below)
nodes.stats_conventional_below



write.csv(nodes.stats_conventional_below, file = "node_check_conventional_below.csv")
nodes.stats_conventional_below <- read.csv("node_check_conventional_below.csv")
library(ggrepel)

between_conventional_below <- ggplot() + 
  geom_point(data = nodes.stats_conventional_below, aes( x = betweeness, y = degree), alpha = 0.6) +
  scale_size_continuous(name = "Relative Abundance") +
  theme_bw() + 
  geom_text_repel(data = subset(nodes.stats_conventional_below, betweeness > 10^hubline.between_conventional_below & degree > hubline.degree_conventional_below), aes(x = betweeness, y = degree, label = OTU)) +
  xlab("Betweeness Centrality") + 
  ylab("Degree") +
  geom_vline(xintercept = 10^hubline.between_conventional_below, linetype = 2) + 
  geom_hline(yintercept = hubline.degree_conventional_below, linetype = 2)

# Hubs were identified by between_conventional_below graph and by being in top 10% of hubscores for fungi or bacteria (assesed outside of R)


### add pos-neg to table - for edges
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


#taxonomy to add to node stats table (outside of R)
write.csv(tax_table(conventional_merged_below), file = "conventional_tax.csv")

library(RCy3)
#exporting files to create cytoscape project. Will also combine node stats file with modules files
createNetworkFromIgraph(plot_spiec_funpro_conventional_below,"conventional_net_below")

#Node stats - import this first.
write.csv(nodes.stats_conventional_below, "conventional_node_stats_below.csv")

# edge stats
write.csv(edge.stats_conventional_below, "conventional_edge_stats_below.csv")

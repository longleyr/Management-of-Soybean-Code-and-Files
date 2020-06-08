#Title: Crop Management Impacts the Soy (Glycine max) Microbiome
#Authors: Reid Longley, Zachary A. Noel, Gian Maria Niccolo Benucci, Martin I. Chilvers, Frances Trail, Gregory Bonito 
#Journal: Submitted to Frontiers in Microbiology

#Prokaryote Script - Reid Longley


options(scipen = 999)
set.seed(9279)
# load packages ---------------------
install.packages("phyloseq")
install.packages("colorspace")
library(colorspace)
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

physeq_object <- phyloseq(prok_otus_phy, 
                          prok_metadata_phy, 
                          prok_taxonomy_phy,
                          prok_sequences)
physeq_object
tax_table(physeq_object)

#removing samples from soil dataset which were noted to not dry properly, and had mold growing in envelope. were analyzed, but only produced a few otus
#non organic samples on this list are being switched with the matching control sample from obj 2 
otu_table(physeq_object) <- subset(otu_table(physeq_object),
                                   select = -c(T4R5AR2S,T4R6BR2S,T4R2CR2S,T4R5CR2S,T2R6BR2S))
###formatting taxonomy------------------------
colnames(tax_table(physeq_object)) = c(k = "Kingdom", p = "Phylum", c = "Class", o = "Order", f = "Family", g = "Genus", s = "Species")
tax_table(physeq_object)

tax_table(physeq_object)[, "Kingdom"] <- gsub("k:", "", tax_table(physeq_object)[, "Kingdom"])
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
# must split by sample source so that contaminants are removed from each individual run
# using decontnam package to remove contaminants in negative controls by comparing
# distributions to those found in normal samples
physeq_object_stems <- subset_samples(physeq_object, origin%in%c("stem"))
physeq_object_leaves <- subset_samples(physeq_object,origin%in%c("leaves"))
physeq_object_roots <- subset_samples(physeq_object,origin%in%c("root"))
physeq_object_soil <- subset_samples(physeq_object,origin%in%c("soil"))
sample_data(physeq_object_soil)
library(devtools)
library(processx)
devtools::install_github("benjjneb/decontam", force = TRUE)
library(decontam)
sample_data(physeq_object_leaves)
write.csv(sample_data(physeq_object_leaves),file ="sample_check.csv")

#leaves

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
ps.noncontam_leaves # with contaminants removed
otu_table(ps.noncontam_leaves)

#roots


# check library size distribution
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

# soil

# check library size distribution
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
ps.noncontam_soil# with contaminants removed
otu_table(ps.noncontam_soil)
sample_data(ps.noncontam_soil)
# remove negative controls from soil
ps.noncontam_soil <- subset_samples(ps.noncontam_soil, experiment%in%c("obj_1","obj_2"))
otu_table(ps.noncontam_soil) <- otu_table(ps.noncontam_soil)[which(rowSums(otu_table(ps.noncontam_soil)) >= 1),]
ps.noncontam_soil
otu_table(ps.noncontam_soil)
sample_data(ps.noncontam_soil)
ps.noncontam_soil

# export soil otu table to check samples
# removing samples with less than 1000 reads, including samples that had less than 1000 read following sampling
write.csv(otu_table(ps.noncontam_soil), file = "filtering_low_soil.csv")


otu_table(ps.noncontam_soil) <- subset(otu_table(ps.noncontam_soil),
                                  select = -c(T1R6CCR3S,T1R6CBR3S,T1R2CBR3S,T2R1CCR3S,T1R5CR2S
))

sample_data(ps.noncontam_soil)
# export roots otu table to check samples
write.csv(otu_table(ps.noncontam_roots), file = "filtering_low_roots.csv")


otu_table(ps.noncontam_roots) <- subset(otu_table(ps.noncontam_roots),
                                       select = -c(T1R5FAR3R, T1R6CR2R, NEG11R, NEG8R, NEG2R, T1R2AV2R, T1R5CAR6R, T1R6CCR4R, T2R1CAR6R, T1R5FCR6R, T1R6CBR6R))

# export leaves otu table to check samples
write.csv(otu_table(ps.noncontam_leaves), file = "filtering_low_leaves.csv")


otu_table(ps.noncontam_leaves) <- subset(otu_table(ps.noncontam_leaves),
                                        select = -c(T2R6FBR3L,NEG17L,T4R1CV2L,T1R1BV2L,NEG9L,NEG16L,Neg14L,T4R1BR6L,T1R5AR2L,NEG10L,NEG12L,NEG5L,T4R5BV2L,NEG15L,Neg2L,NEG3L,Neg8L,NEG13L,NEG11L,T1R1BR2L, T4R2BV2L, T4R1AV2L))

# export stems otu table to check samples
write.csv(otu_table(ps.noncontam_stems), file = "filtering_low_stems.csv")


otu_table(ps.noncontam_stems) <- subset(otu_table(ps.noncontam_stems),
                                         select = -c(NEG10ST,T4R5BR6ST,NEG13ST,T4R1CV2ST,NEG6ST,T4R6AR6ST,T1R1CBR3ST,NEG15ST,NEG14ST,NEG5ST,T1R2AR6ST,T1R2BR6ST,T4R1AR6ST,NEG2ST,NEG3ST,T4R6CR2ST,NEG9ST,NEG11ST,T2R6FCR6ST,T4R1BR6ST,T1R2FBR4ST, T1R5AV2ST))


###for alpha diversity need data which includes singletons etc. will split the current datasets by objective in order to get read distributions for each objective-----------------------


# soil read distributions
ps.noncontam_soil_obj1 <- subset_samples(ps.noncontam_soil, experiment%in%c("obj_1"))
sample_data(ps.noncontam_soil_obj1)


sums_soil_obj1 <- data.frame(colSums(otu_table(ps.noncontam_soil_obj1)))
colnames(sums_soil_obj1) <- "Sample_totalSeqs_soil"
sums_soil_obj1$Sample <- row.names(sums_soil_obj1)



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
# alpha diversity plots
### alpha diversity
# making graph following reviewers reccomendations
library("gridExtra")
library("grid")
library("cowplot")

alpha_supp <- ps.noncontam_soil_obj1
otu_prok <- as.data.frame(otu_table(alpha_supp ))

otu_prok
meta_prok <- as.data.frame(sample_data(alpha_supp ))
alpha_prok <- meta_prok
alpha_prok
alpha_prok$readNO <- sample_sums(alpha_supp)
alpha_prok$Observed <- specnumber(otu_prok, MARGIN = 2)
alpha_prok$Shannon <- diversity(otu_prok, index="shannon", MARGIN = 2)
jevenness_prok <- diversityresult(t(otu_prok), method = "each site", index = "Jevenness")
alpha_prok$Jevenness <- jevenness_prok$Jevenness
#alpha_prok <- alpha_prok[order(alpha_prok$ReadNO), ]
alpha_prok



alpha_prok  
alpha_prok$alpha_label <- factor(alpha_prok$alpha_label,
                                  level=c("Conventional V2","Conventional R2","Conventional R6","No_Till V2","No_Till R2", "No_Till R6","Organic V2", "Organic R2", "Organic R6"))
p <- ggplot(alpha_prok, aes(x=alpha_label, y=readNO,color=Growth_Stage)) + 
  theme_classic()+
  scale_colour_manual("Growth_Stage",breaks = c("V2","R2","R6"),
                      values = c("V2"="orange", "R2"="blue", "R6" = "red")) +
  geom_point(size = 2, shape = 16) +
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold"))+
  geom_boxplot()
p


label_names <- c(Observed="Richness", Shannon="Shannon")
label_names


ps.noncontam_alpha_soil <- ps.noncontam_soil_obj1
sample_data(ps.noncontam_alpha_soil)$alpha_label <- factor(sample_data(ps.noncontam_alpha_soil)$alpha_label,
                      level=c("Conventional V2","Conventional R2","Conventional R6","No_Till V2","No_Till R2", "No_Till R6","Organic V2", "Organic R2", "Organic R6"))
estimate_richness(ps.noncontam_alpha_soil, split = TRUE, measures = NULL)
alpha_soil_prok = plot_richness(ps.noncontam_alpha_soil, x= "alpha_label", 
                                color="Growth_Stage", measures = c( "Observed")) +
  
  ylim(0,9000) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Observed OTUs") +
  scale_x_discrete("Sample", labels = c("Conventional V2" = "Conventional V2","No_Till V2" = "No-Till V2","Organic V2" = "Organic V2","Conventional R2" = "Conventional R2", "No_Till R2" = "No-Till R2","Organic R2" = "Organic R2", "Conventional R6" = "Conventional R6","No_Till R6" = "No-Till R6", "Organic R6" = "Organic R6" )) +
  scale_colour_manual("Growth_Stage",breaks = c("V2","R2","R6"),
  values = c("V2"="orange", "R2"="blue", "R6" = "red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 

  theme(axis.text.x = element_text(angle = 90,vjust =1.5)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  #theme(strip.text.x = element_text(size = 9)) +
  
  theme_set(theme_classic())+
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  
  theme(legend.title=element_blank())
plot(alpha_soil_prok)

alpha_soil_prok_shan = plot_richness(ps.noncontam_alpha_soil, x= "alpha_label", 
                                color="Growth_Stage", measures = c( "Shannon")) +
  
  ylim(1,8) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Shannon Diversity") +
  scale_x_discrete("Sample", labels = c("Conventional V2" = "Conventional V2","No_Till V2" = "No-Till V2","Organic V2" = "Organic V2","Conventional R2" = "Conventional R2", "No_Till R2" = "No-Till R2","Organic R2" = "Organic R2", "Conventional R6" = "Conventional R6","No_Till R6" = "No-Till R6", "Organic R6" = "Organic R6" )) +
  scale_colour_manual("Growth_Stage",breaks = c("V2","R2","R6"),
                      values = c("V2"="orange", "R2"="blue", "R6" = "red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90,vjust =1.5)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  #theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="bottom") +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  
  theme(legend.title=element_blank())
plot(alpha_soil_prok_shan)

# roots
ps.noncontam_alpha_roots <- ps.noncontam_roots_obj1
sample_data(ps.noncontam_alpha_roots)$alpha_label <- factor(sample_data(ps.noncontam_alpha_roots)$alpha_label,
                                                           level=c("Conventional V2","Conventional R2","Conventional R6","No_Till V2","No_Till R2", "No_Till R6","Organic V2", "Organic R2", "Organic R6"))
estimate_richness(ps.noncontam_alpha_roots, split = TRUE, measures = NULL)
alpha_roots_prok = plot_richness(ps.noncontam_alpha_roots, x= "alpha_label", 
                                color="Growth_Stage", measures = c( "Observed")) +
  
  ylim(0,9000)+
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Observed OTUs") +
  scale_x_discrete("Sample", labels = c("Conventional V2" = "Conventional V2","No_Till V2" = "No-Till V2","Organic V2" = "Organic V2","Conventional R2" = "Conventional R2", "No_Till R2" = "No-Till R2","Organic R2" = "Organic R2", "Conventional R6" = "Conventional R6","No_Till R6" = "No-Till R6", "Organic R6" = "Organic R6" )) +
  scale_colour_manual("Growth_Stage",breaks = c("V2","R2","R6"),
                      values = c("V2"="orange", "R2"="blue", "R6" = "red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90,vjust =1.5)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  #theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="none") +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  theme(legend.title=element_blank())
plot(alpha_roots_prok)


alpha_roots_prok_shan = plot_richness(ps.noncontam_alpha_roots, x= "alpha_label", 
                                 color="Growth_Stage", measures = c( "Shannon")) +
  
  ylim(1,8)+
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Shannon Diversity") +
  scale_x_discrete("Sample", labels = c("Conventional V2" = "Conventional V2","No_Till V2" = "No-Till V2","Organic V2" = "Organic V2","Conventional R2" = "Conventional R2", "No_Till R2" = "No-Till R2","Organic R2" = "Organic R2", "Conventional R6" = "Conventional R6","No_Till R6" = "No-Till R6", "Organic R6" = "Organic R6" )) +
  scale_colour_manual("Growth_Stage",breaks = c("V2","R2","R6"),
                      values = c("V2"="orange", "R2"="blue", "R6" = "red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90,vjust =1.5)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  #theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="none") +
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.title=element_blank())
plot(alpha_roots_prok_shan)
#stems
ps.noncontam_alpha_stems <- ps.noncontam_stems_obj1
sample_data(ps.noncontam_alpha_stems)$alpha_label <- factor(sample_data(ps.noncontam_alpha_stems)$alpha_label,
                                                           level=c("Conventional V2","Conventional R2","Conventional R6","No_Till V2","No_Till R2", "No_Till R6","Organic V2", "Organic R2", "Organic R6"))
estimate_richness(ps.noncontam_alpha_stems, split = TRUE, measures = NULL)
alpha_stems_prok = plot_richness(ps.noncontam_alpha_stems, x= "alpha_label", 
                                color="Growth_Stage", measures = c("Observed")) +
  
  ylim(0, 1400) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Observed OTUs") +
  scale_x_discrete("Sample", labels = c("Conventional V2" = "Conventional V2","No_Till V2" = "No-Till V2","Organic V2" = "Organic V2","Conventional R2" = "Conventional R2", "No_Till R2" = "No-Till R2","Organic R2" = "Organic R2", "Conventional R6" = "Conventional R6","No_Till R6" = "No-Till R6", "Organic R6" = "Organic R6" )) +
  scale_colour_manual("Growth_Stage",breaks = c("V2","R2","R6"),
                      values = c("V2"="orange", "R2"="blue", "R6" = "red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  #theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="bottom") +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  theme(legend.title=element_blank())
plot(alpha_stems_prok)

alpha_stems_prok_shan = plot_richness(ps.noncontam_alpha_stems, x= "alpha_label", 
                                 color="Growth_Stage", measures = c("Shannon")) +
  
  ylim(0,6) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Shannon Diversity") +
  scale_x_discrete("Sample", labels = c("Conventional V2" = "Conventional V2","No_Till V2" = "No-Till V2","Organic V2" = "Organic V2","Conventional R2" = "Conventional R2", "No_Till R2" = "No-Till R2","Organic R2" = "Organic R2", "Conventional R6" = "Conventional R6","No_Till R6" = "No-Till R6", "Organic R6" = "Organic R6" )) +
  scale_colour_manual("Growth_Stage",breaks = c("V2","R2","R6"),
                      values = c("V2"="orange", "R2"="blue", "R6" = "red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  #theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="bottom") +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  theme(legend.title=element_blank())
plot(alpha_stems_prok_shan)

# leaves
ps.noncontam_alpha_leaves <- ps.noncontam_leaves_obj1
sample_data(ps.noncontam_alpha_leaves)$alpha_label <- factor(sample_data(ps.noncontam_alpha_leaves)$alpha_label,
                                                           level=c("Conventional V2","Conventional R2","Conventional R6","No_Till V2","No_Till R2", "No_Till R6","Organic V2", "Organic R2", "Organic R6"))
estimate_richness(ps.noncontam_alpha_leaves, split = TRUE, measures = NULL)
alpha_leaves_prok = plot_richness(ps.noncontam_alpha_leaves, x= "alpha_label", 
                                color="Growth_Stage", measures = c( "Observed")) +
  ylim(0,1400)+
  
  expand_limits(x = 0, y = 0) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Observed OTUs") +
  scale_x_discrete("Sample", labels = c("Conventional V2" = "Conventional V2","No_Till V2" = "No-Till V2","Organic V2" = "Organic V2","Conventional R2" = "Conventional R2", "No_Till R2" = "No-Till R2","Organic R2" = "Organic R2", "Conventional R6" = "Conventional R6","No_Till R6" = "No-Till R6", "Organic R6" = "Organic R6" )) +
  scale_colour_manual("Growth_Stage",breaks = c("V2","R2","R6"),
                      values = c("V2"="orange", "R2"="blue", "R6" = "red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title = element_text(size = 10, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  #theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="bottom") +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  theme(legend.title=element_blank())
plot(alpha_leaves_prok)

sample_data(ps.noncontam_alpha_leaves)$alpha_label <- factor(sample_data(ps.noncontam_alpha_leaves)$alpha_label,
                                                             level=c("Conventional V2","Conventional R2","Conventional R6","No_Till V2","No_Till R2", "No_Till R6","Organic V2", "Organic R2", "Organic R6"))
estimate_richness(ps.noncontam_alpha_leaves, split = TRUE, measures = NULL)
alpha_leaves_prok_shan = plot_richness(ps.noncontam_alpha_leaves, x= "alpha_label", 
                                  color="Growth_Stage", measures = c( "Shannon")) +
  ylim(0,6)+
  
  expand_limits(x = 0, y = 0) +
  geom_boxplot(outlier.colour="black", outlier.fill = "black") +
  geom_point(size = 2, shape = 16) +
  labs(title="", x="", y = "Shannon Diversity") +
  scale_x_discrete("Sample", labels = c("Conventional V2" = "Conventional V2","No_Till V2" = "No-Till V2","Organic V2" = "Organic V2","Conventional R2" = "Conventional R2", "No_Till R2" = "No-Till R2","Organic R2" = "Organic R2", "Conventional R6" = "Conventional R6","No_Till R6" = "No-Till R6", "Organic R6" = "Organic R6" )) +
  scale_colour_manual("Growth_Stage",breaks = c("V2","R2","R6"),
                      values = c("V2"="orange", "R2"="blue", "R6" = "red")) +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0)) + 
  theme_set(theme_classic())+
  theme(axis.text.x = element_text(angle = 90)) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  #theme(strip.text.x = element_text(size = 9)) +
  theme(legend.position="bottom") +
  theme(legend.position="bottom") +
  theme(axis.text.x = element_text(angle = 90,vjust =1.5,size = 11, face = "bold")) +
  theme(axis.text.y = element_text(size = 11, face = "bold")) +
  theme(axis.title = element_text(size = 12, face = "bold")) + 
  theme(legend.key = element_blank(), legend.title = element_text(size = 12)) +
  theme(strip.text.x = element_text(size = 12, face = "bold")) +
  theme(legend.position="none") +
  theme(legend.title=element_blank())
plot(alpha_leaves_prok_shan)

# combine plots
library(ggpubr)
ggarrange(alpha_soil_prok,alpha_soil_prok_shan,alpha_roots_prok,alpha_roots_prok_shan,alpha_stems_prok,alpha_stems_prok_shan,alpha_leaves_prok,alpha_leaves_prok_shan,
          #labels = c("A", "B", "C", "D" ),
          widths = c(2.0, 2.0, 2.0, 2.0,2.0,2.0,2.0,2.0),
          align = "h", ncol = 4, nrow = 2)

###alpha diversity_table---------------------------------------------------
# make data frames of metadata and otu tables for analyzing alpha diversity 
otu_soil_prok <- as.data.frame(otu_table(ps.noncontam_soil_obj1))
meta_soil_prok <- as.data.frame(sample_data(ps.noncontam_soil_obj1))
otu_roots_prok <- as.data.frame(otu_table(ps.noncontam_roots_obj1))
meta_roots_prok <- as.data.frame(sample_data(ps.noncontam_roots_obj1))
otu_stems_prok <- as.data.frame(otu_table(ps.noncontam_stems_obj1))
meta_stems_prok <- as.data.frame(sample_data(ps.noncontam_stems_obj1))
otu_leaves_prok <- as.data.frame(otu_table(ps.noncontam_leaves_obj1))
meta_leaves_prok <- as.data.frame(sample_data(ps.noncontam_leaves_obj1))
meta_soil_prok
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

#soil
alpha_div_soil_prok <- meta_soil_prok
alpha_div_soil_prok
alpha_div_soil_prok$readNO <- sample_sums(ps.noncontam_soil_obj1)
alpha_div_soil_prok$Observed <- specnumber(otu_soil_prok, MARGIN = 2)
alpha_div_soil_prok$Shannon <- diversity(otu_soil_prok, index="shannon", MARGIN = 2)
jevenness_prok <- diversityresult(t(otu_soil_prok), method = "each site", index = "Jevenness")
alpha_div_soil_prok$Jevenness <- jevenness_prok$Jevenness
#alpha_div_soil_prok <- alpha_div_soil_prok[order(alpha_div_soil_prok$ReadNO), ]
alpha_div_soil_prok



# get descriptive stats
install.packages("psych")
library("psych")
describeBy(alpha_div_soil_prok, alpha_div_soil_prok$Indicator_label)
warnings()


#roots
alpha_div_roots_prok <- meta_roots_prok
alpha_div_roots_prok
alpha_div_roots_prok$readNO <- sample_sums(ps.noncontam_roots_obj1)
alpha_div_roots_prok$Observed <- specnumber(otu_roots_prok, MARGIN = 2)
alpha_div_roots_prok$Shannon <- diversity(otu_roots_prok, index="shannon", MARGIN = 2)
jevenness_prok <- diversityresult(t(otu_roots_prok), method = "each site", index = "Jevenness")
alpha_div_roots_prok$Jevenness <- jevenness_prok$Jevenness
#alpha_div_roots_prok <- alpha_div_roots_prok[order(alpha_div_roots_prok$ReadNO), ]
alpha_div_roots_prok

# get descriptive stats
describeBy(alpha_div_roots_prok, alpha_div_roots_prok$Indicator_label)

#stems
alpha_div_stems_prok <- meta_stems_prok
alpha_div_stems_prok
alpha_div_stems_prok$readNO <- sample_sums(ps.noncontam_stems_obj1)
alpha_div_stems_prok$Observed <- specnumber(otu_stems_prok, MARGIN = 2)
alpha_div_stems_prok$Shannon <- diversity(otu_stems_prok, index="shannon", MARGIN = 2)
jevenness_prok <- diversityresult(t(otu_stems_prok), method = "each site", index = "Jevenness")
alpha_div_stems_prok$Jevenness <- jevenness_prok$Jevenness
#alpha_div_stems_prok <- alpha_div_stems_prok[order(alpha_div_stems_prok$ReadNO), ]
alpha_div_stems_prok

# get descriptive stats
describeBy(alpha_div_stems_prok, alpha_div_stems_prok$Indicator_label)

#leaves
alpha_div_leaves_prok <- meta_leaves_prok
alpha_div_leaves_prok
alpha_div_leaves_prok$readNO <- sample_sums(ps.noncontam_leaves_obj1)
alpha_div_leaves_prok$Observed <- specnumber(otu_leaves_prok, MARGIN = 2)
alpha_div_leaves_prok$Shannon <- diversity(otu_leaves_prok, index="shannon", MARGIN = 2)
jevenness_prok <- diversityresult(t(otu_leaves_prok), method = "each site", index = "Jevenness")
alpha_div_leaves_prok$Jevenness <- jevenness_prok$Jevenness
#alpha_div_leaves_prok <- alpha_div_leaves_prok[order(alpha_div_leaves_prok$ReadNO), ]
alpha_div_leaves_prok

# get descriptive stats
describeBy(alpha_div_leaves_prok, alpha_div_leaves_prok$Indicator_label)


#alpha_div_soil_prok_df <- read.csv(alpha_soil_)

# fix dataframes so they are recognized
write.csv(alpha_div_soil_prok, file = "alpha_soil_check.csv")
alpha_div_soil_prok_df <- read.csv("alpha_soil_check.csv")
write.csv(alpha_div_roots_prok, file = "alpha_roots_prok_check.csv")
alpha_div_roots_prok_df <- read.csv("alpha_roots_prok_check.csv")
write.csv(alpha_div_stems_prok, file = "alpha_stems_prok_check.csv")
alpha_div_stems_prok_df <- read.csv("alpha_stems_prok_check.csv")
write.csv(alpha_div_leaves_prok, file = "alpha_leaves_prok_check.csv")
alpha_div_leaves_prok_df <- read.csv("alpha_leaves_prok_check.csv")
#soil
# get significant differences
install.packages("agricolae")
library("agricolae")
library(dplyr)

alpha_div_soil_prok_df_v2 <- alpha_div_soil_prok_df[alpha_div_soil_prok_df$Growth_Stage == "V2",]
alpha_div_soil_prok_df_v2

#soil
# get significant differences
install.packages("agricolae")
library("agricolae")
library(dplyr)

alpha_div_soil_prok_df_v2 <- alpha_div_soil_prok_df[alpha_div_soil_prok_df$Growth_Stage == "V2",]
alpha_div_soil_prok_df_v2

#soil v2 rich
aov_prok_soil_v2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_soil_prok_df_v2)
summary(aov_prok_soil_v2_rich)
HSD.test(aov_prok_soil_v2_rich, "Indicator_label") -> tukeyHSD_prok_soil_v2_rich
tukeyHSD_prok_soil_v2_rich

#isa_below_fungi_Management_fdr$sign$p.value<-p.adjust(isa_below_fungi_Management_fdr$sign$p.value, "fdr")

kruskal.test(Observed ~ Indicator_label, data = alpha_div_soil_prok_df_v2)
pairwise.wilcox.test(alpha_div_soil_prok_df_v2$Observed, alpha_div_soil_prok_df_v2$Indicator_label,
                     p.adjust.method = "fdr")
#soil R2 rich
alpha_div_soil_prok_df_R2 <- alpha_div_soil_prok_df[alpha_div_soil_prok_df$Growth_Stage == "R2",]
alpha_div_soil_prok_df_R2
aov_prok_soil_R2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_soil_prok_df_R2)
summary(aov_prok_soil_R2_rich)
HSD.test(aov_prok_soil_R2_rich, "Indicator_label") -> tukeyHSD_prok_soil_R2_rich
tukeyHSD_prok_soil_R2_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_soil_prok_df_R2)
pairwise.wilcox.test(alpha_div_soil_prok_df_R2$Observed, alpha_div_soil_prok_df_R2$Indicator_label,
                     p.adjust.method = "fdr")
# soil R6 rich
alpha_div_soil_prok_df_R6 <- alpha_div_soil_prok_df[alpha_div_soil_prok_df$Growth_Stage == "R6",]
alpha_div_soil_prok_df_R6
aov_prok_soil_R6_rich <- aov(Observed ~ Indicator_label, data=alpha_div_soil_prok_df_R6)
summary(aov_prok_soil_R6_rich)
HSD.test(aov_prok_soil_R6_rich, "Indicator_label") -> tukeyHSD_prok_soil_R6_rich
tukeyHSD_prok_soil_R6_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_soil_prok_df_R6)
pairwise.wilcox.test(alpha_div_soil_prok_df_R6$Observed, alpha_div_soil_prok_df_R6$Indicator_label,
                     p.adjust.method = "fdr")

#soil v2 shan
aov_prok_soil_v2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_soil_prok_df_v2)
summary(aov_prok_soil_v2_shan)
HSD.test(aov_prok_soil_v2_shan, "Indicator_label") -> tukeyHSD_prok_soil_v2_shan
tukeyHSD_prok_soil_v2_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_soil_prok_df_v2)
pairwise.wilcox.test(alpha_div_soil_prok_df_v2$Shannon, alpha_div_soil_prok_df_v2$Indicator_label,
                     p.adjust.method = "fdr")
#soil R2 shan
alpha_div_soil_prok_df_R2 <- alpha_div_soil_prok_df[alpha_div_soil_prok_df$Growth_Stage == "R2",]
alpha_div_soil_prok_df_R2
aov_prok_soil_R2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_soil_prok_df_R2)
summary(aov_prok_soil_R2_shan)
HSD.test(aov_prok_soil_R2_shan, "Indicator_label") -> tukeyHSD_prok_soil_R2_shan
tukeyHSD_prok_soil_R2_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_soil_prok_df_R2)
pairwise.wilcox.test(alpha_div_soil_prok_df_R2$Shannon, alpha_div_soil_prok_df_R2$Indicator_label,
                     p.adjust.method = "fdr")
# soil R6 shan
alpha_div_soil_prok_df_R6 <- alpha_div_soil_prok_df[alpha_div_soil_prok_df$Growth_Stage == "R6",]
alpha_div_soil_prok_df_R6
aov_prok_soil_R6_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_soil_prok_df_R6)
summary(aov_prok_soil_R6_shan)
HSD.test(aov_prok_soil_R6_shan, "Indicator_label") -> tukeyHSD_prok_soil_R6_shan
tukeyHSD_prok_soil_R6_shan


kruskal.test(Shannon ~ Indicator_label, data = alpha_div_soil_prok_df_R6)
pairwise.wilcox.test(alpha_div_soil_prok_df_R6$Shannon, alpha_div_soil_prok_df_R6$Indicator_label,
                     p.adjust.method = "fdr")
#soil v2 even
aov_prok_soil_v2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_soil_prok_df_v2)
summary(aov_prok_soil_v2_even)
HSD.test(aov_prok_soil_v2_even, "Indicator_label") -> tukeyHSD_prok_soil_v2_even
tukeyHSD_prok_soil_v2_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_soil_prok_df_v2)
pairwise.wilcox.test(alpha_div_soil_prok_df_v2$Jevenness, alpha_div_soil_prok_df_v2$Indicator_label,
                     p.adjust.method = "fdr")
#soil R2 even
alpha_div_soil_prok_df_R2 <- alpha_div_soil_prok_df[alpha_div_soil_prok_df$Growth_Stage == "R2",]
alpha_div_soil_prok_df_R2
aov_prok_soil_R2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_soil_prok_df_R2)
summary(aov_prok_soil_R2_even)
HSD.test(aov_prok_soil_R2_even, "Indicator_label") -> tukeyHSD_prok_soil_R2_even
tukeyHSD_prok_soil_R2_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_soil_prok_df_R2)
pairwise.wilcox.test(alpha_div_soil_prok_df_R2$Jevenness, alpha_div_soil_prok_df_R2$Indicator_label,
                     p.adjust.method = "fdr")

# soil R6 even
alpha_div_soil_prok_df_R6 <- alpha_div_soil_prok_df[alpha_div_soil_prok_df$Growth_Stage == "R6",]
alpha_div_soil_prok_df_R6
aov_prok_soil_R6_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_soil_prok_df_R6)
summary(aov_prok_soil_R6_even)
HSD.test(aov_prok_soil_R6_even, "Indicator_label") -> tukeyHSD_prok_soil_R6_even
tukeyHSD_prok_soil_R6_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_soil_prok_df_R6)
pairwise.wilcox.test(alpha_div_soil_prok_df_R6$Jevenness, alpha_div_soil_prok_df_R6$Indicator_label,
                     p.adjust.method = "fdr")
#must also split into management systems 
alpha_div_soil_prok_df_Conventional <- alpha_div_soil_prok_df[alpha_div_soil_prok_df$Management == "Conventional",]
alpha_div_soil_prok_df_Conventional
alpha_div_soil_prok_df_No_Till <- alpha_div_soil_prok_df[alpha_div_soil_prok_df$Management == "No-Till",]
alpha_div_soil_prok_df_No_Till
alpha_div_soil_prok_df_Organic <- alpha_div_soil_prok_df[alpha_div_soil_prok_df$Management == "Organic",]
alpha_div_soil_prok_df_Organic

#soil Conventional rich
aov_prok_soil_Conventional_rich <- aov(Observed ~ Indicator_label, data=alpha_div_soil_prok_df_Conventional)
summary(aov_prok_soil_Conventional_rich)
HSD.test(aov_prok_soil_Conventional_rich, "Indicator_label") -> tukeyHSD_prok_soil_Conventional_rich
tukeyHSD_prok_soil_Conventional_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_soil_prok_df_Conventional)
pairwise.wilcox.test(alpha_div_soil_prok_df_Conventional$Observed, alpha_div_soil_prok_df_Conventional$Indicator_label,
                     p.adjust.method = "fdr")
#soil No_Till rich
aov_prok_soil_No_Till_rich <- aov(Observed ~ Indicator_label, data=alpha_div_soil_prok_df_No_Till)
summary(aov_prok_soil_No_Till_rich)
HSD.test(aov_prok_soil_No_Till_rich, "Indicator_label") -> tukeyHSD_prok_soil_No_Till_rich
tukeyHSD_prok_soil_No_Till_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_soil_prok_df_No_Till)
pairwise.wilcox.test(alpha_div_soil_prok_df_No_Till$Observed, alpha_div_soil_prok_df_No_Till$Indicator_label,
                     p.adjust.method = "fdr")

# soil Organic rich
alpha_div_soil_prok_df_Organic
aov_prok_soil_Organic_rich <- aov(Observed ~ Indicator_label, data=alpha_div_soil_prok_df_Organic)
summary(aov_prok_soil_Organic_rich)
HSD.test(aov_prok_soil_Organic_rich, "Indicator_label") -> tukeyHSD_prok_soil_Organic_rich
tukeyHSD_prok_soil_Organic_rich


kruskal.test(Observed ~ Indicator_label, data = alpha_div_soil_prok_df_Organic)
pairwise.wilcox.test(alpha_div_soil_prok_df_Organic$Observed, alpha_div_soil_prok_df_Organic$Indicator_label,
                     p.adjust.method = "fdr")

#soil Conventional shan
aov_prok_soil_Conventional_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_soil_prok_df_Conventional)
summary(aov_prok_soil_Conventional_shan)
HSD.test(aov_prok_soil_Conventional_shan, "Indicator_label") -> tukeyHSD_prok_soil_Conventional_shan
tukeyHSD_prok_soil_Conventional_shan


kruskal.test(Shannon ~ Indicator_label, data = alpha_div_soil_prok_df_Conventional)
pairwise.wilcox.test(alpha_div_soil_prok_df_Conventional$Shannon, alpha_div_soil_prok_df_Conventional$Indicator_label,
                     p.adjust.method = "fdr")
#soil No_Till shan
aov_prok_soil_No_Till_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_soil_prok_df_No_Till)
summary(aov_prok_soil_No_Till_shan)
HSD.test(aov_prok_soil_No_Till_shan, "Indicator_label") -> tukeyHSD_prok_soil_No_Till_shan
tukeyHSD_prok_soil_No_Till_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_soil_prok_df_No_Till)
pairwise.wilcox.test(alpha_div_soil_prok_df_No_Till$Shannon, alpha_div_soil_prok_df_No_Till$Indicator_label,
                     p.adjust.method = "fdr")
# soil Organic shan
alpha_div_soil_prok_df_Organic
aov_prok_soil_Organic_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_soil_prok_df_Organic)
summary(aov_prok_soil_Organic_shan)
HSD.test(aov_prok_soil_Organic_shan, "Indicator_label") -> tukeyHSD_prok_soil_Organic_shan
tukeyHSD_prok_soil_Organic_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_soil_prok_df_Organic)
pairwise.wilcox.test(alpha_div_soil_prok_df_Organic$Shannon, alpha_div_soil_prok_df_Organic$Indicator_label,
                     p.adjust.method = "fdr")

#soil Conventional even
aov_prok_soil_Conventional_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_soil_prok_df_Conventional)
summary(aov_prok_soil_Conventional_even)
HSD.test(aov_prok_soil_Conventional_even, "Indicator_label") -> tukeyHSD_prok_soil_Conventional_even
tukeyHSD_prok_soil_Conventional_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_soil_prok_df_Conventional)
pairwise.wilcox.test(alpha_div_soil_prok_df_Conventional$Jevenness, alpha_div_soil_prok_df_Conventional$Indicator_label,
                     p.adjust.method = "fdr")
#soil No_Till even
aov_prok_soil_No_Till_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_soil_prok_df_No_Till)
summary(aov_prok_soil_No_Till_even)
HSD.test(aov_prok_soil_No_Till_even, "Indicator_label") -> tukeyHSD_prok_soil_No_Till_even
tukeyHSD_prok_soil_No_Till_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_soil_prok_df_No_Till)
pairwise.wilcox.test(alpha_div_soil_prok_df_No_Till$Jevenness, alpha_div_soil_prok_df_No_Till$Indicator_label,
                     p.adjust.method = "fdr")
# soil Organic even
alpha_div_soil_prok_df_Organic
aov_prok_soil_Organic_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_soil_prok_df_Organic)
summary(aov_prok_soil_Organic_even)
HSD.test(aov_prok_soil_Organic_even, "Indicator_label") -> tukeyHSD_prok_soil_Organic_even
tukeyHSD_prok_soil_Organic_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_soil_prok_df_Organic)
pairwise.wilcox.test(alpha_div_soil_prok_df_Organic$Jevenness, alpha_div_soil_prok_df_Organic$Indicator_label,
                     p.adjust.method = "fdr")






alpha_div_roots_prok_df_v2 <- alpha_div_roots_prok_df[alpha_div_roots_prok_df$Growth_Stage == "V2",]
alpha_div_roots_prok_df_v2

#roots v2 rich
aov_prok_roots_v2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_roots_prok_df_v2)
summary(aov_prok_roots_v2_rich)
HSD.test(aov_prok_roots_v2_rich, "Indicator_label") -> tukeyHSD_prok_roots_v2_rich
tukeyHSD_prok_roots_v2_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_roots_prok_df_v2)
pairwise.wilcox.test(alpha_div_roots_prok_df_v2$Observed, alpha_div_roots_prok_df_v2$Indicator_label,
                     p.adjust.method = "fdr")
#roots R2 rich
alpha_div_roots_prok_df_R2 <- alpha_div_roots_prok_df[alpha_div_roots_prok_df$Growth_Stage == "R2",]
alpha_div_roots_prok_df_R2
aov_prok_roots_R2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_roots_prok_df_R2)
summary(aov_prok_roots_R2_rich)
HSD.test(aov_prok_roots_R2_rich, "Indicator_label") -> tukeyHSD_prok_roots_R2_rich
tukeyHSD_prok_roots_R2_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_roots_prok_df_R2)
pairwise.wilcox.test(alpha_div_roots_prok_df_R2$Observed, alpha_div_roots_prok_df_R2$Indicator_label,
                     p.adjust.method = "fdr")
# roots R6 rich
alpha_div_roots_prok_df_R6 <- alpha_div_roots_prok_df[alpha_div_roots_prok_df$Growth_Stage == "R6",]
alpha_div_roots_prok_df_R6
aov_prok_roots_R6_rich <- aov(Observed ~ Indicator_label, data=alpha_div_roots_prok_df_R6)
summary(aov_prok_roots_R6_rich)
HSD.test(aov_prok_roots_R6_rich, "Indicator_label") -> tukeyHSD_prok_roots_R6_rich
tukeyHSD_prok_roots_R6_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_roots_prok_df_R6)
pairwise.wilcox.test(alpha_div_roots_prok_df_R6$Observed, alpha_div_roots_prok_df_R6$Indicator_label,
                     p.adjust.method = "fdr")

#roots v2 shan
aov_prok_roots_v2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_roots_prok_df_v2)
summary(aov_prok_roots_v2_shan)
HSD.test(aov_prok_roots_v2_shan, "Indicator_label") -> tukeyHSD_prok_roots_v2_shan
tukeyHSD_prok_roots_v2_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_roots_prok_df_v2)
pairwise.wilcox.test(alpha_div_roots_prok_df_v2$Shannon, alpha_div_roots_prok_df_v2$Indicator_label,
                     p.adjust.method = "fdr")
#roots R2 shan
alpha_div_roots_prok_df_R2 <- alpha_div_roots_prok_df[alpha_div_roots_prok_df$Growth_Stage == "R2",]
alpha_div_roots_prok_df_R2
aov_prok_roots_R2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_roots_prok_df_R2)
summary(aov_prok_roots_R2_shan)
HSD.test(aov_prok_roots_R2_shan, "Indicator_label") -> tukeyHSD_prok_roots_R2_shan
tukeyHSD_prok_roots_R2_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_roots_prok_df_R2)
pairwise.wilcox.test(alpha_div_roots_prok_df_R2$Shannon, alpha_div_roots_prok_df_R2$Indicator_label,
                     p.adjust.method = "fdr")
# roots R6 shan
alpha_div_roots_prok_df_R6 <- alpha_div_roots_prok_df[alpha_div_roots_prok_df$Growth_Stage == "R6",]
alpha_div_roots_prok_df_R6
aov_prok_roots_R6_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_roots_prok_df_R6)
summary(aov_prok_roots_R6_shan)
HSD.test(aov_prok_roots_R6_shan, "Indicator_label") -> tukeyHSD_prok_roots_R6_shan
tukeyHSD_prok_roots_R6_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_roots_prok_df_R6)
pairwise.wilcox.test(alpha_div_roots_prok_df_R6$Shannon, alpha_div_roots_prok_df_R6$Indicator_label,
                     p.adjust.method = "fdr")

#roots v2 even
aov_prok_roots_v2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_roots_prok_df_v2)
summary(aov_prok_roots_v2_even)
HSD.test(aov_prok_roots_v2_even, "Indicator_label") -> tukeyHSD_prok_roots_v2_even
tukeyHSD_prok_roots_v2_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_roots_prok_df_v2)
pairwise.wilcox.test(alpha_div_roots_prok_df_v2$Jevenness, alpha_div_roots_prok_df_v2$Indicator_label,
                     p.adjust.method = "fdr")
#roots R2 even
alpha_div_roots_prok_df_R2 <- alpha_div_roots_prok_df[alpha_div_roots_prok_df$Growth_Stage == "R2",]
alpha_div_roots_prok_df_R2
aov_prok_roots_R2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_roots_prok_df_R2)
summary(aov_prok_roots_R2_even)
HSD.test(aov_prok_roots_R2_even, "Indicator_label") -> tukeyHSD_prok_roots_R2_even
tukeyHSD_prok_roots_R2_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_roots_prok_df_R2)
pairwise.wilcox.test(alpha_div_roots_prok_df_R2$Jevenness, alpha_div_roots_prok_df_R2$Indicator_label,
                     p.adjust.method = "fdr")
# roots R6 even
alpha_div_roots_prok_df_R6 <- alpha_div_roots_prok_df[alpha_div_roots_prok_df$Growth_Stage == "R6",]
alpha_div_roots_prok_df_R6
aov_prok_roots_R6_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_roots_prok_df_R6)
summary(aov_prok_roots_R6_even)
HSD.test(aov_prok_roots_R6_even, "Indicator_label") -> tukeyHSD_prok_roots_R6_even
tukeyHSD_prok_roots_R6_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_roots_prok_df_R6)
pairwise.wilcox.test(alpha_div_roots_prok_df_R6$Jevenness, alpha_div_roots_prok_df_R6$Indicator_label,
                     p.adjust.method = "fdr")
#must also split into management systems 
alpha_div_roots_prok_df_Conventional <- alpha_div_roots_prok_df[alpha_div_roots_prok_df$Management == "Conventional",]
alpha_div_roots_prok_df_Conventional
alpha_div_roots_prok_df_No_Till <- alpha_div_roots_prok_df[alpha_div_roots_prok_df$Management == "No-Till",]
alpha_div_roots_prok_df_No_Till
alpha_div_roots_prok_df_Organic <- alpha_div_roots_prok_df[alpha_div_roots_prok_df$Management == "Organic",]
alpha_div_roots_prok_df_Organic

#roots Conventional rich
aov_prok_roots_Conventional_rich <- aov(Observed ~ Indicator_label, data=alpha_div_roots_prok_df_Conventional)
summary(aov_prok_roots_Conventional_rich)
HSD.test(aov_prok_roots_Conventional_rich, "Indicator_label") -> tukeyHSD_prok_roots_Conventional_rich
tukeyHSD_prok_roots_Conventional_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_roots_prok_df_Conventional)
pairwise.wilcox.test(alpha_div_roots_prok_df_Conventional$Observed, alpha_div_roots_prok_df_Conventional$Indicator_label,
                     p.adjust.method = "fdr")
#roots No_Till rich
aov_prok_roots_No_Till_rich <- aov(Observed ~ Indicator_label, data=alpha_div_roots_prok_df_No_Till)
summary(aov_prok_roots_No_Till_rich)
HSD.test(aov_prok_roots_No_Till_rich, "Indicator_label") -> tukeyHSD_prok_roots_No_Till_rich
tukeyHSD_prok_roots_No_Till_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_roots_prok_df_No_Till)
pairwise.wilcox.test(alpha_div_roots_prok_df_No_Till$Observed, alpha_div_roots_prok_df_No_Till$Indicator_label,
                     p.adjust.method = "fdr")
# roots Organic rich
alpha_div_roots_prok_df_Organic
aov_prok_roots_Organic_rich <- aov(Observed ~ Indicator_label, data=alpha_div_roots_prok_df_Organic)
summary(aov_prok_roots_Organic_rich)
HSD.test(aov_prok_roots_Organic_rich, "Indicator_label") -> tukeyHSD_prok_roots_Organic_rich
tukeyHSD_prok_roots_Organic_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_roots_prok_df_Organic)
pairwise.wilcox.test(alpha_div_roots_prok_df_Organic$Observed, alpha_div_roots_prok_df_Organic$Indicator_label,
                     p.adjust.method = "fdr")


#roots Conventional shan
aov_prok_roots_Conventional_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_roots_prok_df_Conventional)
summary(aov_prok_roots_Conventional_shan)
HSD.test(aov_prok_roots_Conventional_shan, "Indicator_label") -> tukeyHSD_prok_roots_Conventional_shan
tukeyHSD_prok_roots_Conventional_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_roots_prok_df_Conventional)
pairwise.wilcox.test(alpha_div_roots_prok_df_Conventional$Shannon, alpha_div_roots_prok_df_Conventional$Indicator_label,
                     p.adjust.method = "fdr")
#roots No_Till shan
aov_prok_roots_No_Till_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_roots_prok_df_No_Till)
summary(aov_prok_roots_No_Till_shan)
HSD.test(aov_prok_roots_No_Till_shan, "Indicator_label") -> tukeyHSD_prok_roots_No_Till_shan
tukeyHSD_prok_roots_No_Till_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_roots_prok_df_No_Till)
pairwise.wilcox.test(alpha_div_roots_prok_df_No_Till$Shannon, alpha_div_roots_prok_df_No_Till$Indicator_label,
                     p.adjust.method = "fdr")

# roots Organic shan
alpha_div_roots_prok_df_Organic
aov_prok_roots_Organic_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_roots_prok_df_Organic)
summary(aov_prok_roots_Organic_shan)
HSD.test(aov_prok_roots_Organic_shan, "Indicator_label") -> tukeyHSD_prok_roots_Organic_shan
tukeyHSD_prok_roots_Organic_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_roots_prok_df_Organic)
pairwise.wilcox.test(alpha_div_roots_prok_df_Organic$Shannon, alpha_div_roots_prok_df_Organic$Indicator_label,
                     p.adjust.method = "fdr")

#roots Conventional even
aov_prok_roots_Conventional_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_roots_prok_df_Conventional)
summary(aov_prok_roots_Conventional_even)
HSD.test(aov_prok_roots_Conventional_even, "Indicator_label") -> tukeyHSD_prok_roots_Conventional_even
tukeyHSD_prok_roots_Conventional_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_roots_prok_df_Conventional)
pairwise.wilcox.test(alpha_div_roots_prok_df_Conventional$Jevenness, alpha_div_roots_prok_df_Conventional$Indicator_label,
                     p.adjust.method = "fdr")
#roots No_Till even
aov_prok_roots_No_Till_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_roots_prok_df_No_Till)
summary(aov_prok_roots_No_Till_even)
HSD.test(aov_prok_roots_No_Till_even, "Indicator_label") -> tukeyHSD_prok_roots_No_Till_even
tukeyHSD_prok_roots_No_Till_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_roots_prok_df_No_Till)
pairwise.wilcox.test(alpha_div_roots_prok_df_No_Till$Jevenness, alpha_div_roots_prok_df_No_Till$Indicator_label,
                     p.adjust.method = "fdr")
# roots Organic even
alpha_div_roots_prok_df_Organic
aov_prok_roots_Organic_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_roots_prok_df_Organic)
summary(aov_prok_roots_Organic_even)
HSD.test(aov_prok_roots_Organic_even, "Indicator_label") -> tukeyHSD_prok_roots_Organic_even
tukeyHSD_prok_roots_Organic_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_roots_prok_df_Organic)
pairwise.wilcox.test(alpha_div_roots_prok_df_Organic$Jevenness, alpha_div_roots_prok_df_Organic$Indicator_label,
                     p.adjust.method = "fdr")








alpha_div_stems_prok_df_v2 <- alpha_div_stems_prok_df[alpha_div_stems_prok_df$Growth_Stage == "V2",]
alpha_div_stems_prok_df_v2

#stems v2 rich
aov_prok_stems_v2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_stems_prok_df_v2)
summary(aov_prok_stems_v2_rich)
HSD.test(aov_prok_stems_v2_rich, "Indicator_label") -> tukeyHSD_prok_stems_v2_rich
tukeyHSD_prok_stems_v2_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_stems_prok_df_v2)
pairwise.wilcox.test(alpha_div_stems_prok_df_v2$Observed, alpha_div_stems_prok_df_v2$Indicator_label,
                     p.adjust.method = "fdr")
#stems R2 rich
alpha_div_stems_prok_df_R2 <- alpha_div_stems_prok_df[alpha_div_stems_prok_df$Growth_Stage == "R2",]
alpha_div_stems_prok_df_R2
aov_prok_stems_R2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_stems_prok_df_R2)
summary(aov_prok_stems_R2_rich)
HSD.test(aov_prok_stems_R2_rich, "Indicator_label") -> tukeyHSD_prok_stems_R2_rich
tukeyHSD_prok_stems_R2_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_stems_prok_df_R2)
pairwise.wilcox.test(alpha_div_stems_prok_df_R2$Observed, alpha_div_stems_prok_df_R2$Indicator_label,
                     p.adjust.method = "fdr")
# stems R6 rich
alpha_div_stems_prok_df_R6 <- alpha_div_stems_prok_df[alpha_div_stems_prok_df$Growth_Stage == "R6",]
alpha_div_stems_prok_df_R6
aov_prok_stems_R6_rich <- aov(Observed ~ Indicator_label, data=alpha_div_stems_prok_df_R6)
summary(aov_prok_stems_R6_rich)
HSD.test(aov_prok_stems_R6_rich, "Indicator_label") -> tukeyHSD_prok_stems_R6_rich
tukeyHSD_prok_stems_R6_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_stems_prok_df_R6)
pairwise.wilcox.test(alpha_div_stems_prok_df_R6$Observed, alpha_div_stems_prok_df_R6$Indicator_label,
                     p.adjust.method = "fdr")


#stems v2 shan
aov_prok_stems_v2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_stems_prok_df_v2)
summary(aov_prok_stems_v2_shan)
HSD.test(aov_prok_stems_v2_shan, "Indicator_label") -> tukeyHSD_prok_stems_v2_shan
tukeyHSD_prok_stems_v2_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_stems_prok_df_v2)
pairwise.wilcox.test(alpha_div_stems_prok_df_v2$Shannon, alpha_div_stems_prok_df_v2$Indicator_label,
                     p.adjust.method = "fdr")
#stems R2 shan
alpha_div_stems_prok_df_R2 <- alpha_div_stems_prok_df[alpha_div_stems_prok_df$Growth_Stage == "R2",]
alpha_div_stems_prok_df_R2
aov_prok_stems_R2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_stems_prok_df_R2)
summary(aov_prok_stems_R2_shan)
HSD.test(aov_prok_stems_R2_shan, "Indicator_label") -> tukeyHSD_prok_stems_R2_shan
tukeyHSD_prok_stems_R2_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_stems_prok_df_R2)
pairwise.wilcox.test(alpha_div_stems_prok_df_R2$Shannon, alpha_div_stems_prok_df_R2$Indicator_label,
                     p.adjust.method = "fdr")
# stems R6 shan
alpha_div_stems_prok_df_R6 <- alpha_div_stems_prok_df[alpha_div_stems_prok_df$Growth_Stage == "R6",]
alpha_div_stems_prok_df_R6
aov_prok_stems_R6_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_stems_prok_df_R6)
summary(aov_prok_stems_R6_shan)
HSD.test(aov_prok_stems_R6_shan, "Indicator_label") -> tukeyHSD_prok_stems_R6_shan
tukeyHSD_prok_stems_R6_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_stems_prok_df_R6)
pairwise.wilcox.test(alpha_div_stems_prok_df_R6$Shannon, alpha_div_stems_prok_df_R6$Indicator_label,
                     p.adjust.method = "fdr")

#stems v2 even
aov_prok_stems_v2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_stems_prok_df_v2)
summary(aov_prok_stems_v2_even)
HSD.test(aov_prok_stems_v2_even, "Indicator_label") -> tukeyHSD_prok_stems_v2_even
tukeyHSD_prok_stems_v2_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_stems_prok_df_v2)
pairwise.wilcox.test(alpha_div_stems_prok_df_v2$Jevenness, alpha_div_stems_prok_df_v2$Indicator_label,
                     p.adjust.method = "fdr")
#stems R2 even
alpha_div_stems_prok_df_R2 <- alpha_div_stems_prok_df[alpha_div_stems_prok_df$Growth_Stage == "R2",]
alpha_div_stems_prok_df_R2
aov_prok_stems_R2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_stems_prok_df_R2)
summary(aov_prok_stems_R2_even)
HSD.test(aov_prok_stems_R2_even, "Indicator_label") -> tukeyHSD_prok_stems_R2_even
tukeyHSD_prok_stems_R2_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_stems_prok_df_R2)
pairwise.wilcox.test(alpha_div_stems_prok_df_R2$Jevenness, alpha_div_stems_prok_df_R2$Indicator_label,
                     p.adjust.method = "fdr")
# stems R6 even
alpha_div_stems_prok_df_R6 <- alpha_div_stems_prok_df[alpha_div_stems_prok_df$Growth_Stage == "R6",]
alpha_div_stems_prok_df_R6
aov_prok_stems_R6_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_stems_prok_df_R6)
summary(aov_prok_stems_R6_even)
HSD.test(aov_prok_stems_R6_even, "Indicator_label") -> tukeyHSD_prok_stems_R6_even
tukeyHSD_prok_stems_R6_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_stems_prok_df_R6)
pairwise.wilcox.test(alpha_div_stems_prok_df_R6$Jevenness, alpha_div_stems_prok_df_R6$Indicator_label,
                     p.adjust.method = "fdr")
#must also split into management systems 
alpha_div_stems_prok_df_Conventional <- alpha_div_stems_prok_df[alpha_div_stems_prok_df$Management == "Conventional",]
alpha_div_stems_prok_df_Conventional
alpha_div_stems_prok_df_No_Till <- alpha_div_stems_prok_df[alpha_div_stems_prok_df$Management == "No-Till",]
alpha_div_stems_prok_df_No_Till
alpha_div_stems_prok_df_Organic <- alpha_div_stems_prok_df[alpha_div_stems_prok_df$Management == "Organic",]
alpha_div_stems_prok_df_Organic

#stems Conventional rich
aov_prok_stems_Conventional_rich <- aov(Observed ~ Indicator_label, data=alpha_div_stems_prok_df_Conventional)
summary(aov_prok_stems_Conventional_rich)
HSD.test(aov_prok_stems_Conventional_rich, "Indicator_label") -> tukeyHSD_prok_stems_Conventional_rich
tukeyHSD_prok_stems_Conventional_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_stems_prok_df_Conventional)
pairwise.wilcox.test(alpha_div_stems_prok_df_Conventional$Observed, alpha_div_stems_prok_df_Conventional$Indicator_label,
                     p.adjust.method = "fdr")
#stems No_Till rich
aov_prok_stems_No_Till_rich <- aov(Observed ~ Indicator_label, data=alpha_div_stems_prok_df_No_Till)
summary(aov_prok_stems_No_Till_rich)
HSD.test(aov_prok_stems_No_Till_rich, "Indicator_label") -> tukeyHSD_prok_stems_No_Till_rich
tukeyHSD_prok_stems_No_Till_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_stems_prok_df_No_Till)
pairwise.wilcox.test(alpha_div_stems_prok_df_No_Till$Observed, alpha_div_stems_prok_df_No_Till$Indicator_label,
                     p.adjust.method = "fdr")
# stems Organic rich
alpha_div_stems_prok_df_Organic
aov_prok_stems_Organic_rich <- aov(Observed ~ Indicator_label, data=alpha_div_stems_prok_df_Organic)
summary(aov_prok_stems_Organic_rich)
HSD.test(aov_prok_stems_Organic_rich, "Indicator_label") -> tukeyHSD_prok_stems_Organic_rich
tukeyHSD_prok_stems_Organic_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_stems_prok_df_Organic)
pairwise.wilcox.test(alpha_div_stems_prok_df_Organic$Observed, alpha_div_stems_prok_df_Organic$Indicator_label,
                     p.adjust.method = "fdr")


#stems Conventional shan
aov_prok_stems_Conventional_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_stems_prok_df_Conventional)
summary(aov_prok_stems_Conventional_shan)
HSD.test(aov_prok_stems_Conventional_shan, "Indicator_label") -> tukeyHSD_prok_stems_Conventional_shan
tukeyHSD_prok_stems_Conventional_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_stems_prok_df_Conventional)
pairwise.wilcox.test(alpha_div_stems_prok_df_Conventional$Shannon, alpha_div_stems_prok_df_Conventional$Indicator_label,
                     p.adjust.method = "fdr")
#stems No_Till shan
aov_prok_stems_No_Till_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_stems_prok_df_No_Till)
summary(aov_prok_stems_No_Till_shan)
HSD.test(aov_prok_stems_No_Till_shan, "Indicator_label") -> tukeyHSD_prok_stems_No_Till_shan
tukeyHSD_prok_stems_No_Till_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_stems_prok_df_No_Till)
pairwise.wilcox.test(alpha_div_stems_prok_df_No_Till$Shannon, alpha_div_stems_prok_df_No_Till$Indicator_label,
                     p.adjust.method = "fdr")
# stems Organic shan
alpha_div_stems_prok_df_Organic
aov_prok_stems_Organic_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_stems_prok_df_Organic)
summary(aov_prok_stems_Organic_shan)
HSD.test(aov_prok_stems_Organic_shan, "Indicator_label") -> tukeyHSD_prok_stems_Organic_shan
tukeyHSD_prok_stems_Organic_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_stems_prok_df_Organic)
pairwise.wilcox.test(alpha_div_stems_prok_df_Organic$Shannon, alpha_div_stems_prok_df_Organic$Indicator_label,
                     p.adjust.method = "fdr")

#stems Conventional even
aov_prok_stems_Conventional_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_stems_prok_df_Conventional)
summary(aov_prok_stems_Conventional_even)
HSD.test(aov_prok_stems_Conventional_even, "Indicator_label") -> tukeyHSD_prok_stems_Conventional_even
tukeyHSD_prok_stems_Conventional_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_stems_prok_df_Conventional)
pairwise.wilcox.test(alpha_div_stems_prok_df_Conventional$Jevenness, alpha_div_stems_prok_df_Conventional$Indicator_label,
                     p.adjust.method = "fdr")
#stems No_Till even
aov_prok_stems_No_Till_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_stems_prok_df_No_Till)
summary(aov_prok_stems_No_Till_even)
HSD.test(aov_prok_stems_No_Till_even, "Indicator_label") -> tukeyHSD_prok_stems_No_Till_even
tukeyHSD_prok_stems_No_Till_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_stems_prok_df_No_Till)
pairwise.wilcox.test(alpha_div_stems_prok_df_No_Till$Jevenness, alpha_div_stems_prok_df_No_Till$Indicator_label,
                     p.adjust.method = "fdr")
# stems Organic even
alpha_div_stems_prok_df_Organic
aov_prok_stems_Organic_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_stems_prok_df_Organic)
summary(aov_prok_stems_Organic_even)
HSD.test(aov_prok_stems_Organic_even, "Indicator_label") -> tukeyHSD_prok_stems_Organic_even
tukeyHSD_prok_stems_Organic_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_stems_prok_df_Organic)
pairwise.wilcox.test(alpha_div_stems_prok_df_Organic$Jevenness, alpha_div_stems_prok_df_Organic$Indicator_label,
                     p.adjust.method = "fdr")










alpha_div_leaves_prok_df_v2 <- alpha_div_leaves_prok_df[alpha_div_leaves_prok_df$Growth_Stage == "V2",]
alpha_div_leaves_prok_df_v2
#leaves v2 rich
aov_prok_leaves_v2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_leaves_prok_df_v2)
summary(aov_prok_leaves_v2_rich)
HSD.test(aov_prok_leaves_v2_rich, "Indicator_label") -> tukeyHSD_prok_leaves_v2_rich
tukeyHSD_prok_leaves_v2_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_leaves_prok_df_v2)
pairwise.wilcox.test(alpha_div_leaves_prok_df_v2$Observed, alpha_div_leaves_prok_df_v2$Indicator_label,
                     p.adjust.method = "fdr")
#leaves R2 rich
alpha_div_leaves_prok_df_R2 <- alpha_div_leaves_prok_df[alpha_div_leaves_prok_df$Growth_Stage == "R2",]
alpha_div_leaves_prok_df_R2
aov_prok_leaves_R2_rich <- aov(Observed ~ Indicator_label, data=alpha_div_leaves_prok_df_R2)
summary(aov_prok_leaves_R2_rich)
HSD.test(aov_prok_leaves_R2_rich, "Indicator_label") -> tukeyHSD_prok_leaves_R2_rich
tukeyHSD_prok_leaves_R2_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_leaves_prok_df_R2)
pairwise.wilcox.test(alpha_div_leaves_prok_df_R2$Observed, alpha_div_leaves_prok_df_R2$Indicator_label,
                     p.adjust.method = "fdr")
# leaves R6 rich
alpha_div_leaves_prok_df_R6 <- alpha_div_leaves_prok_df[alpha_div_leaves_prok_df$Growth_Stage == "R6",]
alpha_div_leaves_prok_df_R6
aov_prok_leaves_R6_rich <- aov(Observed ~ Indicator_label, data=alpha_div_leaves_prok_df_R6)
summary(aov_prok_leaves_R6_rich)
HSD.test(aov_prok_leaves_R6_rich, "Indicator_label") -> tukeyHSD_prok_leaves_R6_rich
tukeyHSD_prok_leaves_R6_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_leaves_prok_df_R6)
pairwise.wilcox.test(alpha_div_leaves_prok_df_R6$Observed, alpha_div_leaves_prok_df_R6$Indicator_label,
                     p.adjust.method = "fdr")


#leaves v2 shan
aov_prok_leaves_v2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_leaves_prok_df_v2)
summary(aov_prok_leaves_v2_shan)
HSD.test(aov_prok_leaves_v2_shan, "Indicator_label") -> tukeyHSD_prok_leaves_v2_shan
tukeyHSD_prok_leaves_v2_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_leaves_prok_df_v2)
pairwise.wilcox.test(alpha_div_leaves_prok_df_v2$Shannon, alpha_div_leaves_prok_df_v2$Indicator_label,
                     p.adjust.method = "fdr")
#leaves R2 shan
alpha_div_leaves_prok_df_R2 <- alpha_div_leaves_prok_df[alpha_div_leaves_prok_df$Growth_Stage == "R2",]
alpha_div_leaves_prok_df_R2
aov_prok_leaves_R2_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_leaves_prok_df_R2)
summary(aov_prok_leaves_R2_shan)
HSD.test(aov_prok_leaves_R2_shan, "Indicator_label") -> tukeyHSD_prok_leaves_R2_shan
tukeyHSD_prok_leaves_R2_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_leaves_prok_df_R2)
pairwise.wilcox.test(alpha_div_leaves_prok_df_R2$Shannon, alpha_div_leaves_prok_df_R2$Indicator_label,
                     p.adjust.method = "fdr")
# leaves R6 shan
alpha_div_leaves_prok_df_R6 <- alpha_div_leaves_prok_df[alpha_div_leaves_prok_df$Growth_Stage == "R6",]
alpha_div_leaves_prok_df_R6
aov_prok_leaves_R6_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_leaves_prok_df_R6)
summary(aov_prok_leaves_R6_shan)
HSD.test(aov_prok_leaves_R6_shan, "Indicator_label") -> tukeyHSD_prok_leaves_R6_shan
tukeyHSD_prok_leaves_R6_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_leaves_prok_df_R6)
pairwise.wilcox.test(alpha_div_leaves_prok_df_R6$Shannon, alpha_div_leaves_prok_df_R6$Indicator_label,
                     p.adjust.method = "fdr")

#leaves v2 even
aov_prok_leaves_v2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_leaves_prok_df_v2)
summary(aov_prok_leaves_v2_even)
HSD.test(aov_prok_leaves_v2_even, "Indicator_label") -> tukeyHSD_prok_leaves_v2_even
tukeyHSD_prok_leaves_v2_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_leaves_prok_df_v2)
pairwise.wilcox.test(alpha_div_leaves_prok_df_v2$Jevenness, alpha_div_leaves_prok_df_v2$Indicator_label,
                     p.adjust.method = "fdr")
#leaves R2 even
alpha_div_leaves_prok_df_R2 <- alpha_div_leaves_prok_df[alpha_div_leaves_prok_df$Growth_Stage == "R2",]
alpha_div_leaves_prok_df_R2
aov_prok_leaves_R2_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_leaves_prok_df_R2)
summary(aov_prok_leaves_R2_even)
HSD.test(aov_prok_leaves_R2_even, "Indicator_label") -> tukeyHSD_prok_leaves_R2_even
tukeyHSD_prok_leaves_R2_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_leaves_prok_df_R2)
pairwise.wilcox.test(alpha_div_leaves_prok_df_R2$Jevenness, alpha_div_leaves_prok_df_R2$Indicator_label,
                     p.adjust.method = "fdr")

# leaves R6 even
alpha_div_leaves_prok_df_R6 <- alpha_div_leaves_prok_df[alpha_div_leaves_prok_df$Growth_Stage == "R6",]
alpha_div_leaves_prok_df_R6
aov_prok_leaves_R6_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_leaves_prok_df_R6)
summary(aov_prok_leaves_R6_even)
HSD.test(aov_prok_leaves_R6_even, "Indicator_label") -> tukeyHSD_prok_leaves_R6_even
tukeyHSD_prok_leaves_R6_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_leaves_prok_df_R6)
pairwise.wilcox.test(alpha_div_leaves_prok_df_R6$Jevenness, alpha_div_leaves_prok_df_R6$Indicator_label,
                     p.adjust.method = "fdr")
#must also split into management syleaves 
alpha_div_leaves_prok_df_Conventional <- alpha_div_leaves_prok_df[alpha_div_leaves_prok_df$Management == "Conventional",]
alpha_div_leaves_prok_df_Conventional
alpha_div_leaves_prok_df_No_Till <- alpha_div_leaves_prok_df[alpha_div_leaves_prok_df$Management == "No-Till",]
alpha_div_leaves_prok_df_No_Till
alpha_div_leaves_prok_df_Organic <- alpha_div_leaves_prok_df[alpha_div_leaves_prok_df$Management == "Organic",]
alpha_div_leaves_prok_df_Organic

#leaves Conventional rich
aov_prok_leaves_Conventional_rich <- aov(Observed ~ Indicator_label, data=alpha_div_leaves_prok_df_Conventional)
summary(aov_prok_leaves_Conventional_rich)
HSD.test(aov_prok_leaves_Conventional_rich, "Indicator_label") -> tukeyHSD_prok_leaves_Conventional_rich
tukeyHSD_prok_leaves_Conventional_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_leaves_prok_df_Conventional)
pairwise.wilcox.test(alpha_div_leaves_prok_df_Conventional$Observed, alpha_div_leaves_prok_df_Conventional$Indicator_label,
                     p.adjust.method = "fdr")
#leaves No_Till rich
aov_prok_leaves_No_Till_rich <- aov(Observed ~ Indicator_label, data=alpha_div_leaves_prok_df_No_Till)
summary(aov_prok_leaves_No_Till_rich)
HSD.test(aov_prok_leaves_No_Till_rich, "Indicator_label") -> tukeyHSD_prok_leaves_No_Till_rich
tukeyHSD_prok_leaves_No_Till_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_leaves_prok_df_No_Till)
pairwise.wilcox.test(alpha_div_leaves_prok_df_No_Till$Observed, alpha_div_leaves_prok_df_No_Till$Indicator_label,
                     p.adjust.method = "fdr")
# leaves Organic rich
alpha_div_leaves_prok_df_Organic
aov_prok_leaves_Organic_rich <- aov(Observed ~ Indicator_label, data=alpha_div_leaves_prok_df_Organic)
summary(aov_prok_leaves_Organic_rich)
HSD.test(aov_prok_leaves_Organic_rich, "Indicator_label") -> tukeyHSD_prok_leaves_Organic_rich
tukeyHSD_prok_leaves_Organic_rich

kruskal.test(Observed ~ Indicator_label, data = alpha_div_leaves_prok_df_Organic)
pairwise.wilcox.test(alpha_div_leaves_prok_df_Organic$Observed, alpha_div_leaves_prok_df_Organic$Indicator_label,
                     p.adjust.method = "fdr")


#leaves Conventional shan
aov_prok_leaves_Conventional_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_leaves_prok_df_Conventional)
summary(aov_prok_leaves_Conventional_shan)
HSD.test(aov_prok_leaves_Conventional_shan, "Indicator_label") -> tukeyHSD_prok_leaves_Conventional_shan
tukeyHSD_prok_leaves_Conventional_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_leaves_prok_df_Conventional)
pairwise.wilcox.test(alpha_div_leaves_prok_df_Conventional$Shannon, alpha_div_leaves_prok_df_Conventional$Indicator_label,
                     p.adjust.method = "fdr")
#leaves No_Till shan
aov_prok_leaves_No_Till_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_leaves_prok_df_No_Till)
summary(aov_prok_leaves_No_Till_shan)
HSD.test(aov_prok_leaves_No_Till_shan, "Indicator_label") -> tukeyHSD_prok_leaves_No_Till_shan
tukeyHSD_prok_leaves_No_Till_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_leaves_prok_df_No_Till)
pairwise.wilcox.test(alpha_div_leaves_prok_df_No_Till$Shannon, alpha_div_leaves_prok_df_No_Till$Indicator_label,
                     p.adjust.method = "fdr")
# leaves Organic shan
alpha_div_leaves_prok_df_Organic
aov_prok_leaves_Organic_shan <- aov(Shannon ~ Indicator_label, data=alpha_div_leaves_prok_df_Organic)
summary(aov_prok_leaves_Organic_shan)
HSD.test(aov_prok_leaves_Organic_shan, "Indicator_label") -> tukeyHSD_prok_leaves_Organic_shan
tukeyHSD_prok_leaves_Organic_shan

kruskal.test(Shannon ~ Indicator_label, data = alpha_div_leaves_prok_df_Organic)
pairwise.wilcox.test(alpha_div_leaves_prok_df_Organic$Shannon, alpha_div_leaves_prok_df_Organic$Indicator_label,
                     p.adjust.method = "fdr")

#leaves Conventional even
aov_prok_leaves_Conventional_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_leaves_prok_df_Conventional)
summary(aov_prok_leaves_Conventional_even)
HSD.test(aov_prok_leaves_Conventional_even, "Indicator_label") -> tukeyHSD_prok_leaves_Conventional_even
tukeyHSD_prok_leaves_Conventional_even

kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_leaves_prok_df_Conventional)
pairwise.wilcox.test(alpha_div_leaves_prok_df_Conventional$Jevenness, alpha_div_leaves_prok_df_Conventional$Indicator_label,
                     p.adjust.method = "fdr")
#leaves No_Till even
aov_prok_leaves_No_Till_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_leaves_prok_df_No_Till)
summary(aov_prok_leaves_No_Till_even)
HSD.test(aov_prok_leaves_No_Till_even, "Indicator_label") -> tukeyHSD_prok_leaves_No_Till_even
tukeyHSD_prok_leaves_No_Till_even


kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_leaves_prok_df_No_Till)
pairwise.wilcox.test(alpha_div_leaves_prok_df_No_Till$Jevenness, alpha_div_leaves_prok_df_No_Till$Indicator_label,
                     p.adjust.method = "fdr")

# leaves Organic even
alpha_div_leaves_prok_df_Organic
aov_prok_leaves_Organic_even <- aov(Jevenness ~ Indicator_label, data=alpha_div_leaves_prok_df_Organic)
summary(aov_prok_leaves_Organic_even)
HSD.test(aov_prok_leaves_Organic_even, "Indicator_label") -> tukeyHSD_prok_leaves_Organic_even
tukeyHSD_prok_leaves_Organic_even


kruskal.test(Jevenness ~ Indicator_label, data = alpha_div_leaves_prok_df_Organic)
pairwise.wilcox.test(alpha_div_leaves_prok_df_Organic$Jevenness, alpha_div_leaves_prok_df_Organic$Indicator_label,
                     p.adjust.method = "fdr")





# filtering otus---------------------------------------------------------
# will filter now before creating barplots

# any sample with less than 5 reads for a particular otu will be placed to 0
ps.noncontam_soil_obj1
otu_table(ps.noncontam_soil_obj1)[otu_table(ps.noncontam_soil_obj1) <= 4] <- 0 ### tag switching
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
sample_data(ps.noncontam_soil_obj1)


# check to make sure additioonal filtering didn't make any samples go below 1000 reads
write.csv(otu_table(ps.noncontam_stems_obj1), file = "stems_check.csv")
write.csv(otu_table(ps.noncontam_roots_obj1), file = "roots_check.csv")
write.csv(otu_table(ps.noncontam_soil_obj1), file = "soil_check.csv")
write.csv(otu_table(ps.noncontam_leaves_obj1), file = "leaves_check.csv")


# no additional samples below 1000




###barplots by prokaryotic class for soil, genus for others---------------------------
library(data.table)
library(dplyr)


#soil
soil_obj1_barplots <- merge_samples(ps.noncontam_soil_obj1, "bar_label")
sample_data(soil_obj1_barplots)$bar_label <- factor(sample_data(soil_obj1_barplots)$bar_label,
                                                  levels = c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))
write.csv(sample_data(soil_obj1_barplots), file ="sample_data_soil_bp.csv")

soil_obj1_barplots <- soil_obj1_barplots %>%
  tax_glom(taxrank = "Class") %>%                     # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Class)           # Sort data frame alphabetically by class
soil_obj1_barplots


dat_soil_bp <- data.table(soil_obj1_barplots)
dat_soil_bp[(Abundance <= 0.02), Class:= "Other"]

# Plot
bar_16s_soil= ggplot(dat_soil_bp, aes(x = Sample, y = Abundance, fill = Class)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))+
  scale_fill_manual(values = c("Thermoleophilia" ="#CBD588",
                               "Spartobacteria" = "#5F7FC7", 
                               "KD4-96" = "orange",
                               "Gemmatimonadetes" = "#DA5724",
                               "Gammaproteobacteria" = "#508578",
                               "Cyanobacteria" = "#CD9BCD",
                               "Betaproteobacteria" = "#AD6F3B",
                               "Alphaproteobacteria" = "#673770", 
                               "Actinobacteria" = "#D14285",
                               "Acidobacteria" = "#652926",
                               "Bradyrhizobium" = "#C84248",
                               "Aureimonas" = "#8569D5",
                               "Bacteroidia" = "steelblue",
                               "Other" = "blue",
                               "Unclassified" = "yellow",
                               "Wolbachia" = "#5E738F",
                               "Nocardioides" = "#D1A33D",
                               "Modestobacter" = "#8A7C64",
                               "Marmoricola" = "#599861",
                               "Dyadobacter" = "coral1",
                               "Arthrobacter" = "greenyellow",
                               "Streptomyces" = "wheat",
                               "Rhizobium" = "magenta1",
                               "Phenylobacterium" = "seagreen2",
                               "Novosphingobium" = "darksalmon",
                               "Niastella" = "firebrick2",
                               "Haemophilus" = "lightgoldenrod2",
                               "Fusobacterium" = "darkmagenta",
                               "Burkholderia" = "lightgrey",
                               "Amycolatopsis" = "olivedrab4",
                               "Deltaproteobacteria" = "yellow",
                               "Sphingobacteriia" = "darkviolet",
                               "TK10" = "pink",
                               "Chloroflexia" = "gray40",
                               "Holophagae" = "aquamarine2",
                               "Planctomycetacia" = "yellowgreen",
                               "Soil_Crenarchaeotic" = "navyblue"
                               
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
  ylab("Relative Abundance (Classes > 2%) \n") +
  ggtitle("Soil")
plot(bar_16s_soil)


#roots
roots_obj1_barplots <- merge_samples(ps.noncontam_roots_obj1, "bar_label")
sample_data(roots_obj1_barplots)$bar_label <- factor(sample_data(roots_obj1_barplots)$bar_label,
                                                  levels = c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))
write.csv(sample_data(roots_obj1_barplots), file ="sample_data_roots_bp.csv")

roots_obj1_barplots <- roots_obj1_barplots %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by class
roots_obj1_barplots


dat_roots_bp <- data.table(roots_obj1_barplots)
dat_roots_bp[(Abundance <= 0.02), Genus := "Other"]

# Plot
bar_16s_roots= ggplot(dat_roots_bp, aes(x = Sample, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))+
  scale_fill_manual(values = c("Methylobacterium" ="#CBD588",
                               "Sphingomonas" = "#5F7FC7", 
                               "Pseudomonas" = "orange",
                               "Streptococcus" = "green",
                               "Paracocccus" = "#508578",
                               "Massilia" = "#CD9BCD",
                               "Hymenobacter" = "#AD6F3B",
                               "Exiguobacterium" = "#673770", 
                               "Enterobacter" = "#D14285",
                               "Buchnera" = "#652926",
                               "Bradyrhizobium" = "#C84248",
                               "Aureimonas" = "#8569D5",
                               "Bacteroidia" = "steelblue",
                               "Other" = "blue",
                               "Unclassified" = "yellow",
                               "Wolbachia" = "#5E738F",
                               "Nocardioides" = "#D1A33D",
                               "Modestobacter" = "#8A7C64",
                               "Marmoricola" = "#599861",
                               "Dyadobacter" = "coral1",
                               "Arthrobacter" = "greenyellow",
                               "Streptomyces" = "wheat",
                               "Rhizobium" = "magenta1",
                               "Phenylobacterium" = "seagreen2",
                               "Novosphingobium" = "darksalmon",
                               "Niastella" = "firebrick2",
                               "Haemophilus" = "lightgoldenrod2",
                               "Fusobacterium" = "darkmagenta",
                               "Burkholderia" = "lightgrey",
                               "Amycolatopsis" = "olivedrab4"
                               
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
  ylab("Relative Abundance (Genera > 2%) \n") +
  ggtitle("Roots")
plot(bar_16s_roots)

  




#stems
stems_obj1_barplots <- merge_samples(ps.noncontam_stems_obj1, "bar_label")
sample_data(stems_obj1_barplots)$bar_label <- factor(sample_data(stems_obj1_barplots)$bar_label,
                                                   levels = c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))
write.csv(sample_data(stems_obj1_barplots), file ="sample_data_stems_bp.csv")

stems_obj1_barplots <- stems_obj1_barplots %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by class
stems_obj1_barplots


dat_stems_bp <- data.table(stems_obj1_barplots)
dat_stems_bp[(Abundance <= 0.02), Genus := "Other"]

# Plot
bar_16s_stems= ggplot(dat_stems_bp, aes(x = Sample, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))+
  scale_fill_manual(values = c("Methylobacterium" ="#CBD588",
                               "Sphingomonas" = "#5F7FC7", 
                               "Pseudomonas" = "orange",
                               "Streptococcus" = "#DA5724",
                               "Paracocccus" = "#508578",
                               "Massilia" = "#CD9BCD",
                               "Hymenobacter" = "#AD6F3B",
                               "Exiguobacterium" = "#673770", 
                               "Enterobacter" = "#D14285",
                               "Buchnera" = "#652926",
                               "Bradyrhizobium" = "#C84248",
                               "Aureimonas" = "#8569D5",
                               "Bacteroidia" = "steelblue",
                               "Other" = "blue",
                               "Unclassified" = "yellow",
                               "Wolbachia" = "#5E738F",
                               "Nocardioides" = "#D1A33D",
                               "Modestobacter" = "#8A7C64",
                               "Marmoricola" = "#599861",
                               "Dyadobacter" = "coral1",
                               "Arthrobacter" = "greenyellow",
                               "Deinococcus" = "lightskyblue",
                               "Roseomonas" = "cornsilk1",
                               "Kineococcus" = "gray40",
                               "Spirosoma" = "aquamarine2"
                               
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
  ylab("Relative Abundance (Genera > 2%) \n") +
  ggtitle("Stems")
plot(bar_16s_stems)


#leaves
leaves_obj1_barplots <- merge_samples(ps.noncontam_leaves_obj1, "bar_label")
sample_data(leaves_obj1_barplots)$bar_label <- factor(sample_data(leaves_obj1_barplots)$bar_label,
                                                  levels = c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))
write.csv(sample_data(leaves_obj1_barplots), file ="sample_data_leaves_bp.csv")

leaves_obj1_barplots <- leaves_obj1_barplots %>%
  tax_glom(taxrank = "Genus") %>%                     # agglomerate at class level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance
  psmelt() %>%                                         # Melt to long format                       # Filter out low abundance taxa
  arrange(Genus)           # Sort data frame alphabetically by class
leaves_obj1_barplots


dat_leaves_bp <- data.table(leaves_obj1_barplots)
dat_leaves_bp[(Abundance <= 0.02), Genus := "Other"]

# Plot
bar_16s_leaves= ggplot(dat_leaves_bp, aes(x = Sample, y = Abundance, fill = Genus)) + 
  #facet_wrap(~Stage, strip.position = "bottom") +
  theme(axis.text.x = element_text(angle = 90))+
  geom_bar(stat = "identity") +
  scale_x_discrete(limits= c("Conventional V2", "Conventional R2", "Conventional R6", "No-Till V2", "No-Till R2", "No-Till R6","Organic V2", "Organic R2", "Organic R6"))+
  scale_fill_manual(values = c("Methylobacterium" ="#CBD588",
                               "Sphingomonas" = "#5F7FC7", 
                               "Pseudomonas" = "orange",
                               "Streptococcus" = "#DA5724",
                               "Paracocccus" = "#508578",
                               "Massilia" = "#CD9BCD",
                               "Hymenobacter" = "#AD6F3B",
                               "Exiguobacterium" = "#673770", 
                               "Enterobacter" = "#D14285",
                               "Buchnera" = "#652926",
                               "Bradyrhizobium" = "#C84248",
                               "Aureimonas" = "#8569D5",
                               "Bacteroidia" = "steelblue",
                               "Other" = "blue",
                               "Unclassified" = "yellow"
                               
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
  ylab("Relative Abundance (Genera > 2%) \n")+
  ggtitle("Leaves")
plot(bar_16s_leaves)



#combine barplots - To produce Figure 2
library(ggpubr)
ggarrange(bar_16s_soil,bar_16s_roots,bar_16s_stems,bar_16s_leaves,
          labels = c("A", "B" ,"C","D"),
          widths = c(2.0, 2.0,2.0,2.0),
          align = "h", ncol = 4, nrow = 1)

sample_data(ps.noncontam_soil_obj1)
sample_data(ps.noncontam_roots_obj1)
sample_data(ps.noncontam_stems_obj1)
sample_data(ps.noncontam_leaves_obj1)

# assesing community composition, retrieving abundances-------------------
# This code was used to get taxonomic percentages shown in the paper
# just change taxonomic level and management or replace management with Growth_Stage.
# check abundances both by growth stage and management
ps.noncontam_soil_obj1_Phylum <- subset_samples(ps.noncontam_soil_obj1,Management%in%c("Conventional", " Conventional"))
sample_data(ps.noncontam_soil_obj1_Phylum)
ps.noncontam_soil_obj1_traditional_Phylum= tax_glom(ps.noncontam_soil_obj1_Phylum, "Phylum")
otu_table(ps.noncontam_soil_obj1_traditional_Phylum)
tax_table <- tax_table(ps.noncontam_soil_obj1_traditional_Phylum)
tax_table
ps.noncontam_soil_obj1_traditional_Phylum= taxa_sums(ps.noncontam_soil_obj1_traditional_Phylum)/sum(taxa_sums(ps.noncontam_soil_obj1_traditional_Phylum))*100
ps.noncontam_soil_obj1_traditional_Phylum
ps.noncontam_soil_obj1_traditional_Phylum<- as.data.frame(ps.noncontam_soil_obj1_traditional_Phylum)
dim(ps.noncontam_soil_obj1_traditional_Phylum)
ps.noncontam_soil_obj1_traditional_Phylum<- setNames(ps.noncontam_soil_obj1_traditional_Phylum, c("abundance"))
# writing into excel to find to most abundant genera, will check with tax table above
write.csv(ps.noncontam_soil_obj1_traditional_Phylum, file = "traditional_Phylum.csv")


#remerge phyloseqs for beta diversity analysis and indicator species analysis
write.csv(tax_table(ps.noncontam_soil_obj1), "tax_check_soil.csv")
ps.noncontam_total_obj1= merge_phyloseq(ps.noncontam_soil_obj1, ps.noncontam_soil_obj1, ps.noncontam_roots_obj1, ps.noncontam_leaves_obj1,ps.noncontam_stems_obj1)
ps.noncontam_total_obj1
sample_data(ps.noncontam_total_obj1)
ps.noncontam_above = merge_phyloseq(ps.noncontam_stems_obj1,ps.noncontam_leaves_obj1)
ps.noncontam_below = merge_phyloseq(ps.noncontam_soil_obj1,ps.noncontam_roots_obj1)



# indicator species analysis 
library("indicspecies")
library("ComplexHeatmap")
library("circlize")

#do above and below to make two heatmaps per fungal/prok dataset

#above 
#merging samples to combine replicate plots within a single Management at a singel growth stage for indicator species analysis
sample_data(ps.noncontam_above)
ps.noncontam_above_merged <- merge_samples(ps.noncontam_above, "Indicator_label")
sample_data(ps.noncontam_above_merged)


# make phyloseq components into dataframes
otu_above_obj1 <- as.data.frame(otu_table(ps.noncontam_above))
otu_above_obj1
tax_above_obj1 <- as.data.frame(as.matrix(tax_table(ps.noncontam_above)))
metadata_above_obj1 <- as.data.frame(as.matrix(sample_data(ps.noncontam_above)))
metadata_above_obj1
# perform indicator species analysis
isa_above_fungi <- multipatt(as.data.frame(t(otu_above_obj1)), metadata_above_obj1$Management_Indicator, control=how(nperm=9999))
summary(isa_above_fungi, indvalcomp=TRUE)
isa_above_fungi -> isa_above_fungi_Management_fdr
isa_above_fungi_Management_fdr$sign$p.value<-p.adjust(isa_above_fungi_Management_fdr$sign$p.value, "fdr")
isa_above_fungi_Management_fdr
summary(isa_above_fungi_Management_fdr)



sink(file="isa_above_fungi_Management_Origin.csv") 
summary(isa_above_fungi_Management_fdr)
sink()
isa_above_fungi_Management_fdr
# extracting ISA OTUs ----------------------------------------------------------------------------
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
ps.noncontam_above_merged -> ps_above_Management.isa
ps_above_Management.isa
ps_above_Management.isa = transform_sample_counts(ps_above_Management.isa, function(x) 100*x/sum(x)) # transform to relativa abundances 
otu_table(ps_above_Management.isa) = otu_table(t(ps_above_Management.isa))
otu_table(ps_above_Management.isa)
otu_table(ps_above_Management.isa) <-otu_table(ps_above_Management.isa)[rownames(isa_above_Management_df),]
ps_above_Management.isa
sample_data(ps_above_Management.isa)
ps_above_Management.isa

install.packages("reltools")







### will determine if taxonomy needs to be improved later
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
# taxonomy was corrected for otus that weren't classified to the genus level
# samples were also renamed using the ps_above_management file
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
sample_data(ps.noncontam_below)
ps.noncontam_below_merged <- merge_samples(ps.noncontam_below, "Indicator_label")
sample_data(ps.noncontam_below_merged)


# make phyloseq components into dataframes
otu_below_obj1 <- as.data.frame(otu_table(ps.noncontam_below))
otu_below_obj1
tax_below_obj1 <- as.data.frame(as.matrix(tax_table(ps.noncontam_below)))
metadata_below_obj1 <- as.data.frame(as.matrix(sample_data(ps.noncontam_below)))
write.csv(metadata_below_obj1, file = "metadata_below_obj1.csv")
metadata_below_obj1 <- read.csv("metadata_below_obj1.csv")
write.csv(otu_below_obj1, file = "otu_below_obj1.csv")
isa_below_fungi <- multipatt(as.data.frame(t(otu_below_obj1)), metadata_below_obj1$Management_Indicator, control=how(nperm=9999))
summary(isa_below_fungi, indvalcomp=TRUE)
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
ps.noncontam_below_merged -> ps_below_Management.isa
ps_below_Management.isa
ps_below_Management.isa = transform_sample_counts(ps_below_Management.isa, function(x) 100*x/sum(x)) # transform to relativa abundances 
otu_table(ps_below_Management.isa) = otu_table(t(ps_below_Management.isa))
otu_table(ps_below_Management.isa)
otu_table(ps_below_Management.isa) <-otu_table(ps_below_Management.isa)[rownames(isa_below_Management_df),]
ps_below_Management.isa
sample_data(ps_below_Management.isa)
ps_below_Management.isa

install.packages("reltools")







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

# writing csv file and manually merging with taxonomy that was updated below in excel
# taxonomy was corrected for otus that weren't classified to the genus level
# samples were also renamed using the ps_above_management file
#write.csv(isa_below_Management_obj, file = "isa_below_add_taxonomy.csv")
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



# These heatmaps as well as the fungal heatmaps were combined to produce figure 4




### normalizing with metagenome seq------------------------------------------------
#NMDS below were combined with fungal NMDS outside of R to produce figure 3
source("http://bioconductor.org/biocLite.R")
biocLite("metagenomeSeq")
library("metagenomeSeq")

#total
# fitting into a Gaussian Model usinFg metagenomeSeq-------------
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
write.csv(otu_table(ps.noncontam_above), file = "empty_sample_check.csv")
ps.noncontam_above_norm = phyloseq_to_metagenomeSeq(ps.noncontam_above)
p_biom_above<-cumNormStat(ps.noncontam_above_norm)
biom_quant_above<-cumNorm(ps.noncontam_above_norm, p=p_biom_above)
biom_quant_above
normFactors(biom_quant_above)
ps.noncontam_above_norm <-MRcounts(biom_quant_above, norm=T)
head(ps.noncontam_above_norm )
#create physeq object with normalized otu table
otu_table(ps.noncontam_above) <- otu_table(ps.noncontam_above_norm, taxa_are_rows = TRUE)

#below
# fitting into a Gaussian Model using metagenomeSeq-------------
write.csv(otu_table(ps.noncontam_below), file = "empty_sample_check.csv")
ps.noncontam_below_norm = phyloseq_to_metagenomeSeq(ps.noncontam_below)
p_biom_below<-cumNormStat(ps.noncontam_below_norm)
biom_quant_below<-cumNorm(ps.noncontam_below_norm, p=p_biom_below)
biom_quant_below
normFactors(biom_quant_below)
ps.noncontam_below_norm <-MRcounts(biom_quant_below, norm=T)
head(ps.noncontam_below_norm )
#create physeq object with normalized otu table
otu_table(ps.noncontam_below) <- otu_table(ps.noncontam_below_norm, taxa_are_rows = TRUE)








#soil
otu_table(ps.noncontam_soil_obj1)

ps.noncontam_soil_norm = phyloseq_to_metagenomeSeq(ps.noncontam_soil_obj1)
p_biom_soil<-cumNormStat(ps.noncontam_soil_norm)
biom_quant_soil<-cumNorm(ps.noncontam_soil_norm, p=p_biom_soil)
biom_quant_soil
normFactors(biom_quant_soil)
ps.noncontam_soil_norm <-MRcounts(biom_quant_soil, norm=T)
head(ps.noncontam_soil_norm )
#create physeq object with normalized otu table
otu_table(ps.noncontam_soil_obj1) <- otu_table(ps.noncontam_soil_norm, taxa_are_rows = TRUE)
otu_table(ps.noncontam_soil_obj1)
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

#leaves
ps.noncontam_leaves_norm = phyloseq_to_metagenomeSeq(ps.noncontam_leaves_obj1)
p_biom_leaves<-cumNormStat(ps.noncontam_leaves_norm)
biom_quant_leaves<-cumNorm(ps.noncontam_leaves_norm, p=p_biom_leaves)
biom_quant_leaves
normFactors(biom_quant_leaves)
ps.noncontam_leaves_norm <-MRcounts(biom_quant_leaves, norm=T)
head(ps.noncontam_leaves_norm )
#create physeq object with normalized otu table
otu_table(ps.noncontam_leaves) <- otu_table(ps.noncontam_leaves_norm, taxa_are_rows = TRUE)






# NMDS - Figure 3 - Prokaryote Plots combined with Fungal NMDS plots outside of R to 
# create final code
library("ggrepel")



# overall bacteria
otu_table(ps.noncontam_total_obj1)
?ordinate()

ord_16s_obj_1 = ordinate(ps.noncontam_total_obj1, method ="NMDS", distance="bray", try=200)
ord_16s_obj_1
otu_table(ps.noncontam_total_obj1)






# shape is stage
NMDS_16s_management = plot_ordination(ps.noncontam_total_obj1, ord_16s_obj_1, color="origin", shape ="Management") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  scale_colour_manual("origin",breaks = c("soil", "root","leaves", "stem"),
                      values = c("soil"="red", "root"="blue", "leaves" = "green", "stem" = "purple"))+
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="right")
plot(NMDS_16s_management)




#below, management
otu_table(ps.noncontam_below)
ord_16s_obj_1_below = ordinate(ps.noncontam_below, method ="NMDS", distance="bray", try=200)
ord_16s_obj_1_below

otu_table(ps.noncontam_below)
NMDS_16S_below = plot_ordination(ps.noncontam_below,ord_16s_obj_1_below,  shape="Management", color ="origin") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  geom_point(size=2.5, alpha=0.9) +
  scale_colour_manual("origin",breaks = c("soil", "root"),
                      values = c("soil"="red", "root"="blue"))+# ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_16S_below)

#above, Growth Stage
otu_table(ps.noncontam_above)
ord_16s_obj_1_above = ordinate(ps.noncontam_above, method ="NMDS", distance="bray", try=1000)
ord_16s_obj_1_above


NMDS_16s_above = plot_ordination(ps.noncontam_above, ord_16s_obj_1_above, shape="Management", color ="origin") + 
  #labs(title="ITS Outdoor NMDS") +
  theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5)) +
  scale_colour_manual("origin",breaks = c("leaves", "stem"),
                      values = c("leaves" = "green", "stem" = "purple"))+
  geom_point(size=2.5, alpha=0.9) + # ,aes(shape=Age))
  #scale_shape_manual(values=1:7) +
  #scale_shape_manual(values=c(0,1,2,3,4,5,6)) +
  theme_classic() +
  #geom_text(aes(label=Description), size = 6) +
  theme(legend.position="none")
plot(NMDS_16s_above)

#pub figure 3 - combined with fungal NMDS outside of R
# to make supplemental fig 1 just switch out Growth_Stage for Management
library(ggpubr)
ggarrange(NMDS_16s_management,NMDS_16s_above, NMDS_16S_below,labels = c("D", "E" ,"F"), ncol =3, nrow=1,widths = c(2.5, 2.0,2.0))






### PERMANOVA (TABLE 3)-------

# creating vegan objects for use in permanova
# will analyze each sample origin individiually to parse out variation by growth stage and management 

#soil
otu_prok_soil <- as.data.frame(otu_table(ps.noncontam_soil_obj1))
taxa_prok_soil <- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_obj1)))
metadata_prok_soil <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_obj1)))
metadata_prok_soil

#roots
otu_prok_roots <- as.data.frame(otu_table(ps.noncontam_roots_obj1))
taxa_prok_roots <- as.data.frame(as.matrix(tax_table(ps.noncontam_roots_obj1)))
metadata_prok_roots <- as.data.frame(as.matrix(sample_data(ps.noncontam_roots_obj1)))

#stems
otu_prok_stems <- as.data.frame(otu_table(ps.noncontam_stems_obj1))
taxa_prok_stems <- as.data.frame(as.matrix(tax_table(ps.noncontam_stems_obj1)))
metadata_prok_stems <- as.data.frame(as.matrix(sample_data(ps.noncontam_stems_obj1)))


#leaves
otu_prok_leaves <- as.data.frame(otu_table(ps.noncontam_leaves_obj1))
taxa_prok_leaves <- as.data.frame(as.matrix(tax_table(ps.noncontam_leaves_obj1)))
metadata_prok_leaves <- as.data.frame(as.matrix(sample_data(ps.noncontam_leaves_obj1)))
#soil permanova
options(scipen = 999) 
library("vegan")
library("RVAideMemoire")

metadata_prok_soil
model.matrix(~ Growth_Stage* Management, data=metadata_prok_soil)
model.matrix(~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_prok_soil)

adonis(t(otu_prok_soil) ~ Growth_Stage * Management, data=metadata_prok_soil, permutations=9999) # by = "margin"

adonis(t(otu_prok_soil) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_prok_soil, permutations=9999) 
adonis(t(otu_prok_soil) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_prok_soil, permutations=9999)


# soil betadisper
vegan::vegdist(t(otu_prok_soil), method="bray") -> dist_otu_prok_soil

permdisp_otu_prok_soil_M <- betadisper(dist_otu_prok_soil, metadata_prok_soil$Management)
permdisp_otu_prok_soil_GS<- betadisper(dist_otu_prok_soil, metadata_prok_soil$Growth_Stage)

anova(permdisp_otu_prok_soil_M, permutations = 9999)
permutest(permdisp_otu_prok_soil_M, permutations = 9999, pairwise = T)
plot(permdisp_otu_prok_soil_M)
plot(TukeyHSD(permdisp_otu_prok_soil_M), las=1)
boxplot(permdisp_otu_prok_soil_M)

anova(permdisp_otu_prok_soil_GS, permutations = 9999)
permutest(permdisp_otu_prok_soil_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_prok_soil_GS)
plot(TukeyHSD(permdisp_otu_prok_soil_GS), las=1)
boxplot(permdisp_otu_prok_soil_GS)

#roots permanova
options(scipen = 999) 
library("vegan")
library("RVAideMemoire")


model.matrix(~ Growth_Stage * Management, data=metadata_prok_roots)
model.matrix(~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_prok_roots)

adonis(t(otu_prok_roots) ~ Growth_Stage * Management, data=metadata_prok_roots, permutations=9999) # by = "margin"

adonis(t(otu_prok_roots) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_prok_roots, permutations=9999) 
adonis(t(otu_prok_roots) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_prok_roots, permutations=9999)


# roots betadisper
vegan::vegdist(t(otu_prok_roots), method="bray") -> dist_otu_prok_roots

permdisp_otu_prok_roots_M <- betadisper(dist_otu_prok_roots, metadata_prok_roots$Management)
permdisp_otu_prok_roots_GS<- betadisper(dist_otu_prok_roots, metadata_prok_roots$Growth_Stage)

anova(permdisp_otu_prok_roots_M, permutations = 9999)
permutest(permdisp_otu_prok_roots_M, permutations = 9999, pairwise = T)
plot(permdisp_otu_prok_roots_M)
plot(TukeyHSD(permdisp_otu_prok_roots_M), las=1)
boxplot(permdisp_otu_prok_roots_M)

anova(permdisp_otu_prok_roots_GS, permutations = 9999)
permutest(permdisp_otu_prok_roots_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_prok_roots_GS)
plot(TukeyHSD(permdisp_otu_prok_roots_GS), las=1)
boxplot(permdisp_otu_prok_roots_GS)


#stems permanova
options(scipen = 999) 
library("vegan")
library("RVAideMemoire")


model.matrix(~ Growth_Stage * Management, data=metadata_prok_stems)
model.matrix(~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_prok_stems)

adonis(t(otu_prok_stems) ~ Growth_Stage * Management, data=metadata_prok_stems, permutations=9999) # by = "margin"

adonis(t(otu_prok_stems) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_prok_stems, permutations=9999) 
adonis(t(otu_prok_stems) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_prok_stems, permutations=9999)


# stems betadisper
vegan::vegdist(t(otu_prok_stems), method="bray") -> dist_otu_prok_stems

permdisp_otu_prok_stems_M <- betadisper(dist_otu_prok_stems, metadata_prok_stems$Management)
permdisp_otu_prok_stems_GS<- betadisper(dist_otu_prok_stems, metadata_prok_stems$Growth_Stage)

anova(permdisp_otu_prok_stems_M, permutations = 9999)
permutest(permdisp_otu_prok_stems_M, permutations = 9999, pairwise = T)
plot(permdisp_otu_prok_stems_M)
plot(TukeyHSD(permdisp_otu_prok_stems_M), las=1)
boxplot(permdisp_otu_prok_stems_M)

anova(permdisp_otu_prok_stems_GS, permutations = 9999)
permutest(permdisp_otu_prok_stems_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_prok_stems_GS)
plot(TukeyHSD(permdisp_otu_prok_stems_GS), las=1)
boxplot(permdisp_otu_prok_stems_GS)


#leaves permanova
options(scipen = 999) 
library("vegan")
library("RVAideMemoire")


model.matrix(~ Growth_Stage * Management, data=metadata_prok_leaves)
model.matrix(~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_prok_leaves)

adonis(t(otu_prok_leaves) ~ Growth_Stage * Management, data=metadata_prok_leaves, permutations=9999) # by = "margin"

adonis(t(otu_prok_leaves) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_prok_leaves, permutations=9999) 
adonis(t(otu_prok_leaves) ~ Growth_Stage + Management + Growth_Stage : Management, data=metadata_prok_leaves, permutations=9999)


# leaves betadisper
vegan::vegdist(t(otu_prok_leaves), method="bray") -> dist_otu_prok_leaves

permdisp_otu_prok_leaves_M <- betadisper(dist_otu_prok_leaves, metadata_prok_leaves$Management)
permdisp_otu_prok_leaves_GS<- betadisper(dist_otu_prok_leaves, metadata_prok_leaves$Growth_Stage)

anova(permdisp_otu_prok_leaves_M, permutations = 9999)
permutest(permdisp_otu_prok_leaves_M, permutations = 9999, pairwise = T)
plot(permdisp_otu_prok_leaves_M)
plot(TukeyHSD(permdisp_otu_prok_leaves_M), las=1)
boxplot(permdisp_otu_prok_leaves_M)

anova(permdisp_otu_prok_leaves_GS, permutations = 9999)
permutest(permdisp_otu_prok_leaves_GS, permutations = 9999, pairwise = T)
plot(permdisp_otu_prok_leaves_GS)
plot(TukeyHSD(permdisp_otu_prok_leaves_GS), las=1)
boxplot(permdisp_otu_prok_leaves_GS)

### Split NMDS for supplemental figure 3

#soil V2
otu_table(ps.noncontam_soil_obj1)
ps.noncontam_soil_obj1_V2 <- subset_samples(ps.noncontam_soil_obj1, Growth_Stage%in%c("V2"))
ord_soil_v2 = ordinate(ps.noncontam_soil_obj1_V2 , method ="NMDS", distance="bray", try=200)
ord_soil_v2 



# shape is management
NMDS_soil_V2= plot_ordination(ps.noncontam_soil_obj1_V2, ord_soil_v2 , color="Management") + 
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
otu_fungi_soil_V2<- as.data.frame(otu_table(ps.noncontam_soil_obj1_V2))
taxa_fungi_soil_V2<- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_obj1_V2)))
metadata_fungi_soil_V2 <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_obj1_V2)))
metadata_fungi_soil_V2
library("vegan")
library("RVAideMemoire")

adonis(t(otu_fungi_soil_V2) ~ Management, data=metadata_fungi_soil_V2, permutations=9999)
vegan::vegdist(t(otu_fungi_soil_V2), method="bray") -> dist_otu_fungi_soil_V2

permdisp_otu_fungi_soil_V2<- betadisper(dist_otu_fungi_soil_V2, metadata_fungi_soil_V2$Management)
anova(permdisp_otu_fungi_soil_V2, permutations = 9999)
#soil R2
otu_table(ps.noncontam_soil_obj1) <- subset(otu_table(ps.noncontam_soil_obj1),
                                            select = -c(T4R6AR2S))
ps.noncontam_soil_obj1_R2 <- subset_samples(ps.noncontam_soil_obj1, Growth_Stage%in%c("R2"))
ord_soil_R2 = ordinate(ps.noncontam_soil_obj1_R2 , method ="NMDS", distance="bray", try=200)
ord_soil_R2 



# shape is management

NMDS_soil_R2= plot_ordination(ps.noncontam_soil_obj1_R2, ord_soil_R2 , color="Management") + 
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
otu_fungi_soil_R2<- as.data.frame(otu_table(ps.noncontam_soil_obj1_R2))
taxa_fungi_soil_R2<- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_obj1_R2)))
metadata_fungi_soil_R2 <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_obj1_R2)))
metadata_fungi_soil_R2


adonis(t(otu_fungi_soil_R2) ~ Management, data=metadata_fungi_soil_R2, permutations=9999)
vegan::vegdist(t(otu_fungi_soil_R2), method="bray") -> dist_otu_fungi_soil_R2

permdisp_otu_fungi_soil_R2<- betadisper(dist_otu_fungi_soil_R2, metadata_fungi_soil_R2$Management)
anova(permdisp_otu_fungi_soil_R2, permutations = 9999)

#soil R6
ps.noncontam_soil_obj1_R6 <- subset_samples(ps.noncontam_soil_obj1, Growth_Stage%in%c("R6"))
ord_soil_R6 = ordinate(ps.noncontam_soil_obj1_R6 , method ="NMDS", distance="bray", try=200)
ord_soil_R6 



# shape is management
NMDS_soil_R6= plot_ordination(ps.noncontam_soil_obj1_R6, ord_soil_R6 , color="Management") + 
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
otu_fungi_soil_R6<- as.data.frame(otu_table(ps.noncontam_soil_obj1_R6))
taxa_fungi_soil_R6<- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_obj1_R6)))
metadata_fungi_soil_R6 <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_obj1_R6)))
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
ps.noncontam_soil_obj1_Conventional <- subset_samples(ps.noncontam_soil_obj1, Management%in%c("Conventional"))
ord_soil_Conventional = ordinate(ps.noncontam_soil_obj1_Conventional , method ="NMDS", distance="bray", try=200)
ord_soil_Conventional 



# shape is management
NMDS_soil_Conventional= plot_ordination(ps.noncontam_soil_obj1_Conventional, ord_soil_Conventional , color="Growth_Stage") + 
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
otu_fungi_soil_Conventional<- as.data.frame(otu_table(ps.noncontam_soil_obj1_Conventional))
taxa_fungi_soil_Conventional<- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_obj1_Conventional)))
metadata_fungi_soil_Conventional <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_obj1_Conventional)))
metadata_fungi_soil_Conventional

adonis(t(otu_fungi_soil_Conventional) ~ Growth_Stage, data=metadata_fungi_soil_Conventional, permutations=9999)
vegan::vegdist(t(otu_fungi_soil_Conventional), method="bray") -> dist_otu_fungi_soil_Conventional

permdisp_otu_fungi_soil_Conventional<- betadisper(dist_otu_fungi_soil_Conventional, metadata_fungi_soil_Conventional$Growth_Stage)
anova(permdisp_otu_fungi_soil_Conventional, permutations = 9999)
#soil No-Till
ps.noncontam_soil_obj1_No_Till <- subset_samples(ps.noncontam_soil_obj1, Management%in%c("No-Till"))
ord_soil_No_Till = ordinate(ps.noncontam_soil_obj1_No_Till , method ="NMDS", distance="bray", try=200)
ord_soil_No_Till 



# shape is Growth_Stage
NMDS_soil_No_Till= plot_ordination(ps.noncontam_soil_obj1_No_Till, ord_soil_No_Till , color="Growth_Stage") + 
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
otu_fungi_soil_No_Till<- as.data.frame(otu_table(ps.noncontam_soil_obj1_No_Till))
taxa_fungi_soil_No_Till<- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_obj1_No_Till)))
metadata_fungi_soil_No_Till <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_obj1_No_Till)))
metadata_fungi_soil_No_Till


adonis(t(otu_fungi_soil_No_Till) ~ Growth_Stage, data=metadata_fungi_soil_No_Till, permutations=9999)
vegan::vegdist(t(otu_fungi_soil_No_Till), method="bray") -> dist_otu_fungi_soil_No_Till

permdisp_otu_fungi_soil_No_Till<- betadisper(dist_otu_fungi_soil_No_Till, metadata_fungi_soil_No_Till$Growth_Stage)
anova(permdisp_otu_fungi_soil_No_Till, permutations = 9999)

#soil Organic

ps.noncontam_soil_obj1_Organic <- subset_samples(ps.noncontam_soil_obj1, Management%in%c("Organic"))
ps.noncontam_soil_obj1_Organic 
# remove 3 outliers to make NMDS visible
otu_table(ps.noncontam_soil_obj1_Organic) <- subset(otu_table(ps.noncontam_soil_obj1_Organic),
                                                    select = -c(T4R1CR6S,T4R2AR6S,T4R2AR2S))


ord_soil_Organic = ordinate(ps.noncontam_soil_obj1_Organic , method ="NMDS", distance="bray", try=200)
ord_soil_Organic 




# shape is management

NMDS_soil_Organic= plot_ordination(ps.noncontam_soil_obj1_Organic, ord_soil_Organic  , color="Growth_Stage") + 
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
otu_fungi_soil_Organic<- as.data.frame(otu_table(ps.noncontam_soil_obj1_Organic))
taxa_fungi_soil_Organic<- as.data.frame(as.matrix(tax_table(ps.noncontam_soil_obj1_Organic)))
metadata_fungi_soil_Organic <- as.data.frame(as.matrix(sample_data(ps.noncontam_soil_obj1_Organic)))
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
# remove 3 outliers for nmds to work
otu_table(ps.noncontam_roots_R6) <- subset(otu_table(ps.noncontam_roots_R6),
                                           select = -c(T1R1BR6R,T1R6AR6R,T4R5AR6R))

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
#remove two outliers 
otu_table(ps.noncontam_roots_No_Till) <- subset(otu_table(ps.noncontam_roots_No_Till),
                                                select = -c(T2R1AR2R,T2R5AR2R))

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
otu_table(ps.noncontam_leaves_obj1_Conventional) <- subset(otu_table(ps.noncontam_leaves_obj1_Conventional),
                                                           select = -c(T1R5AR2L))
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


# Note on supplemental Figure - Each set of 3 was arranged using ggarrange and then they
# were combined outside of R

         

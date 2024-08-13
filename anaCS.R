library(phyloseq)
library(bendiR)
library(ggplot2)
library(ggpubr)
library(readxl)
library(vegan)
library(viridis)
library(dplyr)
library(performance)
library(RColorBrewer)
library(phylogeo)
library(ggordiplots)
library(lme4)
library(piecewiseSEM)
library(semPlot)
library(lavaan)

setwd("/Users/u2061312/Documents/PhD/Countryside survey/ISME submisson/R")
load("/Users/u2061312/Documents/PhD/Countryside survey/ISME submisson/R/CountrySS_master.R")
CountrySS_master

load("/Users/u2061312/Documents/PhD/Countryside survey/ISME submisson/R/Tetracladium_RA_group.R")
Tetracladium_master <- subset_taxa(CountrySS_master, Genus=="Tetracladium")
CountrySS_master_RA <- transform_sample_counts(CountrySS_master, function(x) x/sum(x)*100)
Tetracladium_master_RA <- subset_taxa(CountrySS_master_RA, Genus=="Tetracladium")

#OTU table all samples
data_otu = as(otu_table(phylo2_RA), "matrix")
sum(data_otu)

#OTU table Tetracladium samples
stat_phylo_tetra <- prune_samples(sample_sums(phylo2)>=1,phylo2)
data_otuT = as(otu_table(stat_phylo_tetra), "matrix")
sum(data_otuT)


sample_sums(CountrySS_master_RA)

save(Tetracladium_geo_E_RA,file="/Users/u2061312/Documents/PhD/Countryside survey/R/Tetracladium_geo_E_RA")

load("/Users/u2061312/Documents/PhD/Countryside survey/ISME submisson/R/old2/Tetracladium_geo_E_RA_grouped.R")
load("/Users/u2061312/Documents/PhD/Countryside survey/ISME submisson/R/old2/Tetracladium_geo_E_grouped.R")

#Remove samples with no Broad habitat
sample_variables(Tetracladium_master)
phylo <- subset_samples(Tetracladium_master, Broad.habitat !="")

#Remove samples with no Habitat
phylo2 <- subset_samples(Tetracladium_master, Habitat !="")
phylo2_RA <- subset_samples(Tetracladium_master_RA, Habitat !="")
phylo_group <- subset_samples(Tetracladium_RA_group, Habitat !="")

#Remove samples with no metadata
phylo_R <- subset_samples(Tetracladium_master, Broad.habitat !="")
phylo_R <- subset_samples(phylo_R, N !="")
phylo_R <- subset_samples(phylo_R, pH !="")
phylo_R <- subset_samples(phylo_R, Organic.C !="")
phylo_R <- subset_samples(phylo_R, Moisture !="")
phylo_R <- subset_samples(phylo_R, P !="")
phylo_R <- subset_samples(phylo_R, C !="")
phylo_R <- subset_samples(phylo_R, Ellenburg.N !="")
phylo_R <- subset_samples(phylo_R, avc07 !="")
phylo_R <- subset_samples(phylo_R, Habitat !="")
phylo_R <- subset_samples(phylo_R, Longitude !="")
phylo_R <- subset_samples(phylo_R, Latitude !="")

#REMOVE taxa and samples with 0 no reads  
phylo_R <- prune_taxa(taxa_sums(phylo_R) > 0, phylo_R)
phylo_R <- prune_samples(sample_sums(phylo_R)>=0.000000000000000001,phylo_R)


#Make sure this phyloseq matches the phyloseq you are using.
metadata <- as(sample_data(phylo_R), "data.frame")

#Remove samples with no metadata grouped
sample_variables(Tetracladium_master_RA)
phylo_RA <- subset_samples(Tetracladium_master_RA, Broad.habitat !="")
phylo_RA <- subset_samples(phylo_RA, N !="")
phylo_RA <- subset_samples(phylo_RA, pH !="")
phylo_RA <- subset_samples(phylo_RA, Organic.C !="")
phylo_RA <- subset_samples(phylo_RA, Moisture !="")
phylo_RA <- subset_samples(phylo_RA, P !="")
phylo_RA <- subset_samples(phylo_RA, C !="")
phylo_RA <- subset_samples(phylo_RA, Ellenburg.N !="")
phylo_RA <- subset_samples(phylo_RA, avc07 !="")
phylo_RA <- subset_samples(phylo_RA, Habitat !="")
phylo_RA <- subset_samples(phylo_RA, Longitude !="")
phylo_RA <- subset_samples(phylo_RA, Latitude !="")

#REMOVE taxa and samples with 0 no reads  
phylo_RA <- prune_taxa(taxa_sums(phylo_RA) > 0, phylo_RA)
phylo_RA <- prune_samples(sample_sums(phylo_RA)>=0.000000000000000001,phylo_RA)


#Make sure this phyloseq matches the phyloseq you are using.
metadata_RA <- as(sample_data(phylo_RA), "data.frame")

#rarefy
## Rarefaction curve, ggplot style
ggrare <- function(physeq, step = 10, label = NULL, color = NULL, plot = TRUE, parallel = FALSE, se = TRUE) {
  ## Args:
  ## - physeq: phyloseq class object, from which abundance data are extracted
  ## - step: Step size for sample size in rarefaction curves
  ## - label: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - color: (Optional). Default ‘NULL’. Character string. The name of the
  ##          variable to map to colors in the plot. This can be a sample
  ##          variable (among the set returned by
  ##          ‘sample_variables(physeq)’ ) or taxonomic rank (among the set
  ##          returned by ‘rank_names(physeq)’).
  ##
  ##          Finally, The color scheme is chosen automatically by
  ##          ‘link{ggplot}’, but it can be modified afterward with an
  ##          additional layer using ‘scale_color_manual’.
  ## - color: Default `NULL`. Character string. The name of the variable
  ##          to map to text labels on the plot. Similar to color option
  ##          but for plotting text.
  ## - plot:  Logical, should the graphic be plotted.
  ## - parallel: should rarefaction be parallelized (using parallel framework)
  ## - se:    Default TRUE. Logical. Should standard errors be computed.
  ## require vegan
  x <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) { x <- t(x) }
  
  ## This script is adapted from vegan `rarecurve` function
  tot <- rowSums(x)
  S <- rowSums(x > 0)
  nr <- nrow(x)
  
  rarefun <- function(i) {
    cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
    n <- seq(1, tot[i], by = step)
    if (n[length(n)] != tot[i]) {
      n <- c(n, tot[i])
    }
    y <- rarefy(x[i, ,drop = FALSE], n, se = se)
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
    }
  }
  if (parallel) {
    out <- mclapply(seq_len(nr), rarefun, mc.preschedule = FALSE)
  } else {
    out <- lapply(seq_len(nr), rarefun)
  }
  df <- do.call(rbind, out)
  
  ## Get sample data
  if (!is.null(sample_data(physeq, FALSE))) {
    sdf <- as(sample_data(physeq), "data.frame")
    sdf$Sample <- rownames(sdf)
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x))
    labels <- merge(labels, sdf, by = "Sample")
  }
  
  ## Add, any custom-supplied plot-mapped variables
  if( length(color) > 1 ){
    data$color <- color
    names(data)[names(data)=="color"] <- deparse(substitute(color))
    color <- deparse(substitute(color))
  }
  if( length(label) > 1 ){
    labels$label <- label
    names(labels)[names(labels)=="label"] <- deparse(substitute(label))
    label <- deparse(substitute(label))
  }
  
  p <- ggplot(data = data, aes_string(x = "Size", y = ".S", group = "Sample", color = color))
  p <- p + labs(x = "Sample Size", y = "Species richness")
  if (!is.null(label)) {
    p <- p + geom_text(data = labels, aes_string(x = "x", y = "y", label = label, color = color),
                       size = 4, hjust = 0)
  }
  p <- p + geom_line()
  if (se) { ## add standard error if available
    p <- p + geom_ribbon(aes_string(ymin = ".S - .se", ymax = ".S + .se", color = NULL, fill = color), alpha = 0.2)
  }
  if (plot) {
    plot(p)
  }
  invisible(p)
}


plot_rare <- ggrare(CountrySS_master, step = 100, color = "Habitat", label = NULL, se = FALSE)

ggsave(filename = "RareCurve.pdf", plot = plot_rare,
       scale = 1,
       dpi = 300)


#write OTU ra table
OTU1 = as(otu_table(Tetracladium_master_RA), "matrix")
if(taxa_are_rows(Tetracladium_master_RA)){OTU1 <- t(OTU1)}
OTUdf = as.data.frame(OTU1)
Sample <- rownames(OTUdf)
rownames(OTUdf) <- NULL
OTUdf2 <- cbind(Sample,OTUdf)

write.csv(OTUdf2, "RA_OTUtable_master_R.csv")


#Richness
diversity_plot4<-function(physeq,grouping,diversity=c("fishers","shannon","simpson","observed"),test=c("anova","kruskal"),correction_method=stats::p.adjust.methods){
  if(!phyloseq::taxa_are_rows(physeq)){
    phyloseq::otu_table(physeq)<-t(phyloseq::otu_table(physeq))
  }
  test <- base::match.arg(test,c("anova","kruskal"))
  diversity <- base::match.arg(diversity,c("fishers","shannon","simpson","observed"))
  if(length(correction_method)>1){
    corr_specified<-FALSE
  } else{
    corr_specified<-TRUE
  }
  if(length(test)==2){
    warning("Defaulting to ANOVA")
    test<-"anova"
  }
  if(!length(grep(grouping,colnames(phyloseq::sample_data(physeq))))==0){
    if(diversity=="fishers"){
      div<-vegan::fisher.alpha(t(phyloseq::otu_table(physeq)))
    }  else if(diversity=="observed"){
      div<-colSums(otu_table(physeq) != 0)
    } else {
      div<-vegan::diversity(t(phyloseq::otu_table(physeq)),diversity)
    } 
    if(test=="anova"){
      anova<-stats::aov(div ~ as.character(phyloseq::sample_data(physeq)[[grouping]]))
      post_hoc<-stats::TukeyHSD(anova)
      means<-stats::aggregate(anova$model[, 1], list(anova$model[,2]), mean)
      letters<-compact_letters(phyloseq::sample_data(physeq)[[grouping]],post_hoc)
    } else{
      correction_method <- base::match.arg(correction_method)
      if(corr_specified == FALSE){
        warning("Multiple testing correction post-hoc p-values with FDR")
        correction_method<-"fdr"
      }
      kruskal<-stats::kruskal.test(div ~ as.factor(phyloseq::sample_data(physeq)[[grouping]]))
      post_hoc<-DescTools::DunnTest(div~as.factor(phyloseq::sample_data(physeq)[[grouping]]),method=correction_method)
      diversity_table<-cbind.data.frame(div,group=as.character(phyloseq::sample_data(physeq)[[grouping]]))
      means<-stats::aggregate(diversity_table[,"div"], list(diversity_table[,"group"]), mean)
      letters<-compact_letters(phyloseq::sample_data(physeq)[[grouping]],post_hoc)
    }
    
    names(means)<-c("Group","mean")
    
  } else{
    stop(paste("Grouping factor entered does not exist in phyloseq::sample_data. Check phyloseq::sample_data(",substitute(physeq),").",sep=""))
  }
  if(!exists("diversity_table")){
    diversity_table<-cbind.data.frame(div,group=as.character(phyloseq::sample_data(physeq)[[grouping]]))
  }
  
  if(typeof(letters)=="character"){
    letters<-cbind.data.frame(sample=names(letters),letters)
  } else{
    letters<-cbind.data.frame(sample=names(letters$Letters),letters=letters$Letters)
    print(letters)
  }
  plot<-ggplot2::ggplot(data=diversity_table)+ggplot2::geom_boxplot(aes(x=group,y=div,fill=group, colour=group),alpha=0.3,outlier.shape=NA)+ggplot2::geom_jitter(aes(x=group,y=div, colour=group),alpha=0.5)+ggplot2::geom_text(data=letters,aes(x=sample,y=Inf,label=letters,vjust = 1,fontface="bold"),size=4)+theme_journal()+scale_y_continuous("Observed species")+ theme(legend.title = element_blank())#+
  #theme(axis.title.x=element_blank(),
  #axis.text.x=element_blank(),
  #axis.ticks.x=element_blank())
  return(list(POSTHOC=post_hoc,PLOT=plot))
}

alpha <- diversity_plot4(phylo,"Broad.habitat",diversity="observed",test="kruskal", correction_method = "none") 
plot1 <- alpha$PLOT
plot2 <- plot1  +
  theme_pubr() +
  theme(legend.position="right", axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank(), legend.title = element_blank())

alpha$POSTHOC
ggsave(filename = "Richness.pdf", plot = plot2,
       scale = 1.5,
       dpi = 300)

beta <- diversity_plot4(phylo2,"Habitat",diversity="observed",test="kruskal", correction_method = "none") 
plot1b <- beta$PLOT
plot2b <- plot1b  +
  theme_pubr() +
  theme(legend.position="right", axis.text.x.bottom = element_blank(),
        axis.ticks.x = element_blank(), legend.title = element_blank())

beta$POSTHOC
ggsave(filename = "Richness_agg.pdf", plot = plot2b,
       scale = 1.5,
       dpi = 300)



#PCOA instead of NMDS ph(#E09B54","#604225), mc("#CC8DBD","#563B51")

GP.ord <- ordinate(phylo_R, distance="bray", method="PCoA")
p2 = plot_ordination(phylo_R, GP.ord, type="samples", color="pH",shape="Habitat",label=NULL,axes=c(1, 2))+
  scale_shape_manual(values=c(15,4,17,18,19,2,8,11))
                             
pcoa_plot <- p2 + geom_point(aes(color=pH),size=2) +
  ggtitle("PCoA")+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +labs(fill="pH", shape="Habitat")+
  scale_color_gradientn(colours = c("#E09B54","#604225"))


ggsave(filename = "pcoa_plot_PH.pdf", plot = pcoa_plot,
       scale = 1,
       dpi = 300)

#raupcrick instead of PCOA ph(#E09B54","#604225), mc("#CC8DBD","#563B51")
raupcrick <- raupcrick(otu_table(phylo_R))
p_raupcrick = plot_ordination(phylo_R, raupcrick, type="samples", color="pH",shape="Habitat",label=NULL,axes=c(1, 2))+
  scale_shape_manual(values=c(15,4,17,18,19,2,8,11))

#convert phyloseq object to vegan matrix
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

#print the sample variables
sample_variables(phylo_R)

#get vegan OTU table
votu <- veganotu(phylo_R)

#make sure vegan OTU is binary (this is necessary for RaupCrick)
votu[votu>0] <-1
#recreate binary table as phyloseq object
OTU = otu_table(t(votu), taxa_are_rows = TRUE)
physeq2 = phyloseq(OTU, tax_table(phylo_R), sample_data(phylo_R))
#loop through treatments, keep the  treatment options in biom that are relic of 
keep <- c("pH","Habitat")
treatments <- intersect(sample_variables(phylo_R), keep)
print(treatments)
special_distances <- c("raupcrick")

  print(treatments[y])
  variables = get_variable(phylo_R, treatments[y])
  num_variables = length(unique(variables))
  
  if ( is.element('no_data',variables) ) {
    num_variables <- num_variables - 1
  }
  
  if ( num_variables > 1 ) {
    remove_idx = as.character(get_variable(phylo_R, treatments[y])) != "no_data"
    PhyBray <- filter_taxa(prune_samples(remove_idx, phylo_R), function(x) sum(x) > 1, TRUE)
    PhyRaup <- filter_taxa(prune_samples(remove_idx, physeq2), function(x) sum(x) > 1, TRUE)
    print(PhyBray) }
    
    #get metadata again for updated phyloseq object
    metadata <- as(sample_data(PhyBray), "data.frame")
    
    #convert to vegan OTU matrix for hypothesis test
    votuRaup <- veganotu(PhyRaup)
    votuBray <- veganotu(PhyBray)
    
    #run distance for raupcrick
     {
      d <- raupcrick(votuRaup, null="r1", nsimul=999, chase=FALSE)
      PhyORD = PhyRaup
      veganData = votuRaup
    } 
    
    #hypothesis test for significance
    testA <- adonis2(d~metadata[[treatments[y]]], permutations=9999)
    betaA <- betadisper(d, metadata[[treatments[y]]])
    pA <- permutest(betaA, permutations=9999)
    
    #create dataframe of stats from tests
    Df <- c(testA$aov.tab$Df[1], pA$tab$Df[1])
    Fstat <- c(testA$aov.tab$F.Model[1], pA$tab$F[1])
    R2 <- c(testA$aov.tab$R2[1], 'NA')
    Pvalue <- c(testA$aov.tab$`Pr(>F)`[1], pA$tab$`Pr(>F)`[1])
    xTab <- data.frame(Df,Fstat,R2,Pvalue)
    
    
    #setup color palette, if more than 15 don't even draw it as it is dumb
    num_colors = length(unique(get_variable(PhyBray, treatments[y])))
    print(num_colors)
    getPalette = colorRampPalette(brewer.pal(9, "Set1"))
    if (num_colors < 10 ) {
      colors = brewer.pal(9, "Set1")
    } else if ( num_colors > 15) {
      next
    } else {
      colors = getPalette(num_colors)
    }
    
      
    #run ordination on previously computed distances, then save in a 2x2 image with shared legend.  
    ord <- ordinate(PhyORD,method="PCoA", distance=d)

p2b = plot_ordination(PhyORD, ord, type="samples", color="Habitat",shape="Habitat",label=NULL,axes=c(1, 2))+
      scale_shape_manual(values=c(15,4,17,18,19,2,8,11))  
oedi_plot <- p2b + geom_point(aes(color=Habitat),size=2) +
      ggtitle("NMDS")+
      theme_bw()+
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +labs(fill="Habitat", shape="Habitat")

ggsave(filename = "NMDS_new.pdf", plot = oedi_plot,
       scale = 1,
       dpi = 300)

#PERMANOVA use cleaned RA dataset (NO NAs) this case phylo

#Run this function to convert physeq to vegan format
veganotu = function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

sample_variables(phylo_RA)


#Test each variable (categorical or continuous) separately,

adonis2(veganotu(phylo_RA) ~ C, data = metadata_RA, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA) ~ N, data = metadata_RA, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA) ~ pH, data = metadata_RA, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA) ~ Organic.C, data = metadata_RA, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA) ~ Moisture, data = metadata_RA, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA) ~ P, data = metadata_RA, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA) ~ Habitat, data = metadata_RA, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA) ~ Longitude, data = metadata_RA, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA) ~ Latitude, data = metadata_RA, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA) ~ Ellenburg.N, data = metadata_RA, method="bray", permutations = 100)

#Then put together in order of highest R first.

set.seed(12345)

adonis <- adonis2(veganotu(phylo_RA) ~ Habitat+pH+Longitude+Ellenburg.N+
                   Moisture+N+Organic.C+C+Latitude+P,
                 data = metadata_RA, method="bray", permutations = 100)


adonis

#Adjusted p values. Not sure if needed. n = the number of variables added. Can change to holm, fdr etc. Canoco uses FDR. 
adonis$aov.tab$`Adj_P` <- p.adjust(adonis$aov.tab$`Pr(>F)`, method = 'fdr', n = 11)
#See table with adjusted p values for multiple comparisons included. 

adonis

#Write output to text.

out <- capture.output(adonis)
cat(out, file="adonis.txt", sep="\n")


#adonis for arable only

phylo_RA_Agri <- subset_samples(phylo_RA, Habitat == "Crops and weeds" )
unique(get_variable(phylo_RA_Agri, "Habitat"))


#Make sure this phyloseq matches the phyloseq you are using.
metadataAgri <- as(sample_data(phylo_RA_Agri), "data.frame")


#Test each variable (categorical or continuous) separately,

adonis2(veganotu(phylo_RA_Agri) ~ C, data = metadataAgri, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA_Agri) ~ N, data = metadataAgri, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA_Agri) ~ pH, data = metadataAgri, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA_Agri) ~ Organic.C, data = metadataAgri, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA_Agri) ~ Moisture, data = metadataAgri, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA_Agri) ~ P, data = metadataAgri, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA_Agri) ~ Longitude, data = metadataAgri, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA_Agri) ~ Latitude, data = metadataAgri, method="bray", permutations = 100)
adonis2(veganotu(phylo_RA_Agri) ~ Ellenburg.N, data = metadataAgri, method="bray", permutations = 100)

#Then put together in order of highest R first.

set.seed(12345)

adonisAgri <- adonis2(veganotu(phylo_RA_Agri) ~ pH+Moisture+Latitude+N+Latitude+Organic.C+C+
                        Ellenburg.N+P,
                 data = metadataAgri, method="bray", permutations = 100)


adonisAgri

#Adjusted p values. Not sure if needed. n = the number of variables added. Can change to holm, fdr etc. Canoco uses FDR. 
adonisAgri$aov.tab$`Adj_P` <- p.adjust(adonisAgri$aov.tab$`Pr(>F)`, method = 'fdr', n = 9)
#See table with adjusted p values for multiple comparisons included. 

adonisAgri

#Write output to text.

outAgri <- capture.output(adonisAgri)
cat(outAgri, file="adonisAgri.txt", sep="\n")

##############################RDA
#RDA analyses source: http://www.hiercourse.com/docs/microbial/04_betaDiversity_multTables.html
# extract the OTU table
Y <- veganotu(phylo_RA)

# extract the sample data from the 'phyloseq' object
sample_variables(phylo_RA)
# then remove the non needed columns
dat <- data.frame(sample_data(phylo_RA))
dat <- select(dat, -our.ref)
dat <- select(dat, -REP.ID)
dat <- select(dat, -REP.group)
dat <- select(dat, -Rep)
dat <- select(dat, -Eastings )
dat <- select(dat, -Northings)
dat <- select(dat, -BH)
dat <- select(dat, -H)
dat <- select(dat, -Broad.habitat)
dat <- select(dat, -seq.lookup.bac)
dat <- select(dat, -seq.lookup.fun)
dat <- select(dat, -seq.lookup.euks)
dat <- select(dat, -D.plate)
dat <- select(dat, -Well.row)
dat <- select(dat, -Well.col)
dat <- select(dat, -Band)
dat <- select(dat, -Sample)
dat <- select(dat, -PH07)
dat <- select(dat, -avc07)
dat <- select(dat, -TOT.RICH07)
dat <- select(dat, -NATIVE07)
dat <- select(dat, -NONNATIVE07)
dat <- select(dat, -Ellenburg.pH)
dat <- select(dat, -Ellenburg.wetness)
dat <- select(dat, -Ellenburg.light)
dat <- select(dat, -Competitor.strategy)
dat <- select(dat, -Stress.tolerator.strategy)
dat <- select(dat, -Ruderal.strateg)
dat <- select(dat, -BUTTRICH07)
dat <- select(dat, -CVS07)
dat <- select(dat, -BIRDRICH07)
dat <- select(dat, -AVC071)
dat <- select(dat, -SCPTDATA.ID.07AN)
dat <- select(dat, -OS.2.FIG.10KM)
dat <- select(dat, -Habitat.name)
dat <- select(dat, -HABT.CODE1)


# select the continuous variables to use in the constraint (Y)
# then standardise
X <- select_if(dat, is.numeric)
X <- decostand(X, method='standardize')

# RDA
res <- rda(Y ~ ., data=X)
resplot1 <- plot(res, xlab=ord_labels(res)[1], ylab=ord_labels(res)[2])
summary(res, display=NULL)



#There are multiple approaches to select only the variables that are explaining variation efficiently.
# 1.calculate variance inflation
# variables with scores >10 are redundant
sort(vif.cca(res))

# 2. calculate fit - Fitting environmental vectors/factors onto an ordination 
envfit(Y ~ ., data=X)

# 3. stepwise selection - set up full and null models for 'ordistep'
# full model
rda1 <- rda(Y ~ ., data=X)

# intercept-only (null) model
rda0 <- rda(Y ~ 1, data=X)

# perform forward and backward selection of explanatory variables
# output not shown
step.env <- ordistep(rda0, scope=formula(rda1), direction='both')
# look at the significant variables 
step.env$anova
# code to get variable names from 'ordistep' and 'envfit' results
vars.ordistep <- gsub('^. ', '', rownames(step.env$anova))

vars.envfit <- names(which(vif.cca(res) <= 10))

vars <- unique(c(vars.ordistep, vars.envfit))

# select variables to keep from table 'Y'
X1 <- X[, vars]
str(X1)
#RDA REDUCED VARIABLES
res <- rda(Y ~ ., data=X1)

# summary of the results
summary(res, display=NULL)
anova(res)


# set up dataframes for plotting the results
sit <- cbind(dat, scores(res, display='sites'))
spp <- cbind(data.frame(tax_table(phylo_RA)), scores(res, display='species'))
vec <- data.frame(scores(res, display='bp'))

# use these to adjust length of arrows and position of arrow labels
adj.vec <- 2
adj.txt <- 2.5

# 'site' ordination
p1 <- ggplot(sit, aes(x=RDA1, y=RDA2, color=Habitat, shape=Habitat)) +
  geom_point(size=3) +
  scale_shape_manual(values=c(15,4,17,18,19,2,8,11))+
  geom_segment(data=vec, inherit.aes=F, 
               mapping=aes(x=0, y=0, xend=adj.vec*RDA1, yend=adj.vec*RDA2), 
               arrow=arrow(length=unit(0.2, 'cm'))) + 
  geom_text(data=vec, inherit.aes=F, 
            mapping=aes(x=adj.txt*RDA1, y=adj.txt*RDA2, 
                        label=c('pH', 'Moisture',"Latitude", 'N',"C", "Ellenburg.N","P","Longitude"))) +
  theme_pubr()
p1$labels$x <- ord_labels(res)[1]
p1$labels$y <- ord_labels(res)[2]

ggsave("finalRDA_ordi.pdf", p1, dpi = 300)

#Variation partitioning - more than two tables
# extract the OTU table
#Y <- veganotu(phylo)

# extract the sample data from the 'phyloseq' object
# then remove the 'SampleID' column
#remove cat variable to standardize (AVC add back later)
#dat <- data.frame(sample_data(phylo))
#dat <- select(dat, -our.ref)
#dat <- select(dat, -AVC07)
#dat <- select(dat, -BH)
#dat <- select(dat, -REP.ID07)
#dat <- select(dat, -BROAD.HABITAT)

dat2 <- select(dat, -Habitat)


X = dat2

X <- decostand(X, method='standardize')

# define the partitions
# soil variables from ordistep and envfit analyses
X1 <- X[, vars]
X1 <- select(X1, -Longitude)
X1 <- select(X1, -Latitude)
# habitat variables
X2 <- select(dat, Habitat)
# spatial variables
X3 <- X[, grep('tude', names(X), value=T)]

# partition variation and plot the result
vp <- varpart(Y, X1, X2, X3)
plot(vp, Xnames=c('Soil prop', 'Habitat', 'Location'))

#We can test the significance of each partition by performing RDA and conditioning the OTU matrix with all of the other partitions than the one that we are focussing on in that analysis. For example:
# Test the 'Habitat' partition (X2)
# first argument contains the community data 
# second argument is the partition of interest (constraint)
# third argument is a column-bound dataframe of remaining partitions (condition)
rda.x2 <- rda(Y, X2, cbind(X1, X3))
rda.x1 <- rda(Y, X1, cbind(X2, X3))
rda.x3 <- rda(Y, X3, cbind(X2, X1))

# test significance
anova(rda.x2)
anova(rda.x1)
anova(rda.x3)


#seperate ph and moisture vs habitat
X4 <- select(X, pH)
X5 <- select(X, Moisture)


vp2 <- varpart(Y, X4, X2)
plot(vp2, Xnames=c('pH', "Habitat"))

#seperate all vs all
X4 <- select(X, PH)
X5 <- select(X, MOISTURE.CONTENT)
X6 <- select(X, N..PERCENT)
X7 <- select(X, Latitude)
X8 <- select(X, C.PERCENT)
X9 <- select(X, Ellenburg.N)
X10 <- select(X, PO4.P.OLSEN.MG.KG)
X11 <- select(X, Longitude)


vp2 <- varpart(Y, X4, X9, X5)
plot(vp2, Xnames=c('pH', "EllenN", "Moisture"))

# #plot pH response to OTUOTU62642 abundance faceted in different AVCs
# trybb<-cbind(Y,dat)
# responseplot1 <- ggplot(trybb, aes(x=OTU62642, y=PH))+
#   geom_point(size=2) +
#   facet_wrap(~AGG.CLASS.DESCRIPTION)
# #plot pH response to OTUOTU62642 abundance in different AVCs
# responseplot2 <- ggplot(trybb, aes(x=OTU62642, y=PH, color = AGG.CLASS.DESCRIPTION))+
#   geom_point(size=2) +
#   geom_smooth(method=lm)+
#   theme_pubr()

#plot pH response to OTU diversity in different AVCs
observed_D_samples <- estimate_richness(phylo_R, split = TRUE, measures = "Observed")
observed_D_samples <- filter(observed_D_samples, Observed > 0)
tryaa<-cbind(observed_D_samples,dat)

responseplot2 <- ggplot(tryaa, aes(x=Observed, y=pH, color = Habitat))+
  geom_point(size=2, alpha = 0.3) +
  geom_smooth(method = "lm")+
  theme_pubr()+
  facet_wrap(~Habitat)+
  stat_poly_line() +
  stat_poly_eq(use_label(c("p", "R2")))

ggsave("ph_observed_hab.pdf",responseplot2,  dpi = 300)

responseplot_RA <- ggplot(mod_data, aes(x=x, y=pH, color = Habitat))+
  geom_point(size=2, alpha = 0.3) +
  geom_smooth(method = "lm")+
  theme_pubr()+
  facet_wrap(~Habitat)+
  stat_poly_line() +
  stat_poly_eq(use_label(c("p", "R2")))


ggsave("ph_RA_hab.pdf",responseplot_RA,  dpi = 300)


responseplot2b <- ggplot(tryaa, aes(x=Observed, y=pH, color = Habitat))+
   geom_point(size=2) +
   geom_smooth(method = "lm")+
   theme_pubr()+
   facet_wrap(~Habitat)+
  stat_poly_eq(use_label(c("p", "R2")))

#plot MC response to OTUOTU62642 abundance faceted in different AVCs
# responseplot1MC <- ggplot(trybb, aes(x=OTU62642, y=MOISTURE.CONTENT))+
#   geom_point(size=2) +
#   facet_wrap(~AGG.CLASS.DESCRIPTION)
#plot pH response to OTU62642 abundance in different AVCs
# responseplotMC2 <- ggplot(trybb, aes(x=OTU62642, y=PH, color = AGG.CLASS.DESCRIPTION))+
#   geom_point(size=2) +
#   geom_smooth(method=lm)+
#   theme_pubr()

#plot MC response to OTU diversity in different AVCs
# responseplot2MC <- ggplot(tryaa, aes(x=Observed, y=Moisture, color = Habitat))+
#   geom_point(size=2) +
#   geom_smooth(method = "lm")+
#   theme_pubr()
# 
# responseplot2bMC <- ggplot(tryaa, aes(x=Observed, y=Moisture, color = Habitat))+
#   geom_point(size=2) +
#   geom_smooth(method = "gam")+
#   theme_pubr()+
#   facet_wrap(~Habitat)


###########################MODELING
#grab data (has to be data frame)
#get data OTU RA
Y
#add OTU RA
Ysum <- rowSums(Y,na.rm = FALSE)
#get meta data
dat

mod_data <- merge(Ysum, dat, by ="row.names", all = TRUE )
rownames(mod_data) <- mod_data[,1]
mod_data[,1] <- NULL

#normalise continuous variables
#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

mod_data2 <- mod_data[, c(1,2,3,4,5,6,7,9,10,11,8)]
mod_data2
#apply Min-Max normalization
mod_data_norm <- as.data.frame(lapply(mod_data2[1:10], min_max_norm))

#add cat variable back
mod_data_norm$Habitat <- mod_data2$Habitat

#selected cont variables from before
vars


#AGG.CLASS.DESCRIPTION as random variable, fixed variables based on RDA
mod1 <- lmer(x ~  pH + Ellenburg.N + Moisture  + N + C + P + (1|Habitat), data = mod_data_norm)
tab_model(mod1, show.est=TRUE)
plot_mod1 <- plot_model(mod1, sort.est=FALSE, type = "std", show.values=TRUE, value.size=2, value.offset=0.4)
plot_mod1

mod2 <- lmer(x ~  pH + Moisture * C *N + Ellenburg.N *P + (1|Habitat), data = mod_data_norm)
tab_model(mod2, show.est=TRUE)
plot_mod2 <- plot_model(mod2, sort.est=FALSE, type = "std", show.values=TRUE, value.size=2, value.offset=0.4)
plot_mod2

mod3 <- lmer(x ~  pH  + (1|Habitat), data = mod_data_norm)
tab_model(mod3, show.est=TRUE)
plot_mod3 <- plot_model(mod3, sort.est=FALSE, type = "std", show.values=TRUE, value.size=2, value.offset=0.4)
plot_mod3

comp <- compare_performance(mod1, mod2,mod3,
                                 metrics = "all",
                                 rank = TRUE)
comp.p <- ggtexttable(comp, rows = NULL, 
                        theme = ttheme("lBlack"))

comppolt <- plot(comp)

performance_accuracy(
  mod3,
  method = "boot",
  k = 5,
  n = 1000,
  verbose = TRUE)

CMplot <- check_model(
  mod3,
  dot_size = 2,
  line_size = 0.8,
  panel = TRUE,
  check = "all",
  alpha = 0.2,
  dot_alpha = 0.8,
  colors = c("#3aaf85", "#1b6ca8", "#cd201f"),
  theme = "see::theme_lucid",
  detrend = FALSE,
  verbose = TRUE)

# Use the effects package --> effect function. term= the fixed effect you want to get data on, mod= name of your model.

effects <- effects::effect(term= "pH", mod= mod3)
summary(effects) #output of what the values are
x_pH <- as.data.frame(effects)

pH_plot <- ggplot() + 
  geom_point(data=mod_data_norm, aes(pH, x), colour= "#C46854", alpha = 0.3, size = 3 ) + 
  geom_point(data=x_pH, aes(x=pH, y=fit), color="#C46854") +
  geom_line(data=x_pH, aes(x=pH, y=fit), color="#C46854" ) +
  geom_ribbon(data= x_pH, aes(x=pH, ymin=lower, ymax=upper), alpha= 0.3, fill="#C46854") +
  labs(x="pH", y="Tetracladium total relative abundance") +
  theme_pubr()

pH_plot

ggsave("pH_plot_abund.pdf", pH_plot, dpi = 300)

model_coefs <- coef(mod3)$Habitat %>% 
  rename(Intercept = `(Intercept)`, Slope = pH) %>% 
  rownames_to_column("Habitat")

sleep_groups_rani <- left_join(mod_data_norm, model_coefs, by = "Habitat")

model_coef_plot <- ggplot(data = sleep_groups_rani, 
                          mapping = aes(x = pH, 
                                        y = x, 
                                        colour = Habitat)) +
  geom_point(na.rm = F, alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, 
                  slope = Slope,
                  colour = Habitat),linewidth = 1.5) +
  theme_pubr()




#STACKED BAR PLOTS

#Make barplot

#Merge by phyla, change Phylum to Class (or different ta level) if needed
phyla<-(tax_glom(phylo_group,"Species"))

#Merge by Sample type (column heading of your metadata). Don't worry about NAs.
merged_phyla<-merge_samples(phyla, "Habitat")


#otu_table(merged_phyla)<-otu_table((otu_table(merged_phyla))/(data.frame(sample_data(phyloG)) %>% group_by(AGG.CLASS.DESCRIPTION) %>% tally())$n,taxa_are_rows=FALSE)

#In this one you have the choice of leaving all taxa that are above x% in at least one sample in, or just merging everything low abundance on a sample by sample basis, keep_all=TRUE

#merged_phyla<-merge_low_abun_taxa(merged_phyla,1,keep_all=TRUE)

#Melt for use with ggplot
melted_phyla<-psmelt(merged_phyla)

plot_phyla<-ggplot(data=melted_phyla, 
                   aes(fill=Species,x=Sample, y=Abundance))+
  geom_bar(stat="identity")+
  labs(y="Total abundance",x="") +
  theme_pubr()+
  theme(legend.title=element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))
  

plot_phyla


ggsave("Stacked_bar_total_abundance_habitat.pdf", width = 20, height = 15, units = "cm") 

#STACKED BAR PLOTS manual merge from new phyloseq Tetracladium_geo_group

#Merge by Sample type (column heading of your metadata). Don't worry about NAs.
merged_group_bh <-merge_samples(phylo_group, "Habitat")


otu_table(merged_group_bh)<-otu_table((otu_table(merged_group_bh))/(data.frame(sample_data(phylo_group)) %>% group_by(Habitat) %>% tally())$n,taxa_are_rows=FALSE)

#In this one you have the choice of leaving all taxa that are above x% in at least one sample in, or just merging everything low abundance on a sample by sample basis, keep_all=TRUE

#merged_phyla<-merge_low_abun_taxa(merged_phyla,1,keep_all=TRUE)

#Melt for use with ggplot
melted_group_bh<-psmelt(merged_group_bh)


plot_group_bh<-ggplot(data=melted_group_bh, 
                   aes(fill=Species,x=Sample, y=Abundance))+
  geom_bar(stat="identity")+
  labs(y="Number of reads",x="") +
  theme_pubr()+
  theme(legend.title=element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) 


plot_group_bh


ggsave(plot =plot_group_bh,"Stacked_bar_actual_habitat_new.pdf", width = 20, height = 15, units = "cm") 

#Merge by AGG.CLASS.DESCRIPTION (column heading of your metadata). Don't worry about NAs.
merged_group_acd <-merge_samples(Tetracladium_geo_E_grouped, "AGG.CLASS.DESCRIPTION")


otu_table(merged_group_acd)<-otu_table((otu_table(merged_group_acd))/(data.frame(sample_data(Tetracladium_geo_E_grouped)) %>% group_by(AGG.CLASS.DESCRIPTION) %>% tally())$n,taxa_are_rows=FALSE)

#In this one you have the choice of leaving all taxa that are above x% in at least one sample in, or just merging everything low abundance on a sample by sample basis, keep_all=TRUE

#merged_phyla<-merge_low_abun_taxa(merged_phyla,1,keep_all=TRUE)

#Melt for use with ggplot
melted_group_acd<-psmelt(merged_group_acd)


plot_group_acd<-ggplot(data=melted_group_acd, 
                      aes(fill=Species,x=Sample, y=Abundance))+
  geom_bar(stat="identity")+
  labs(y="Number of reads",x="") +
  theme_pubr()+
  theme(legend.title=element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_cvi_d("CSS_colours") +
  scale_color_cvi_d("CSS_colours")


plot_group_acd


ggsave(plot =plot_group_acd,"Stacked_bar_actual_aggclass_new.pdf", width = 20, height = 15, units = "cm") 
#FACATED bar plots from same melted data as above

facated_bar_tax <- ggplot(melted_group_bh, aes(fill=Species, y=Abundance, x=Sample)) + 
  geom_bar(position="fill", stat="identity")+
  labs(y="Percent of Tetracladium reads",x="") +
  theme_pubr()+
  theme(legend.title=element_blank(),
        legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = c("#B36458","#AC7A40","#7A733B","#4B7544","#4D7E6E","#5995AE","#567CA2","#666990","#765F84","#875B7E","#A05D80","#B15E7D"))

facated_bar_tax
ggsave(plot =facated_bar_tax,"facated_bar_tax.pdf", width = 20, height = 15, units = "cm") 


#Heat map for finding unique OTUs use phyloseq with RA
#REMOVE taxa and samples with 0 no reads 
phylo_RA_removed <- prune_taxa(taxa_sums(phylo2_RA) > 0, phylo2_RA)
phylo_RA_removed <- prune_samples(sample_sums(phylo_RA_removed)>=0.000000000000000001,phylo_RA_removed)


mergedT_RA_all <- merge_samples(phylo_group, "Habitat")

HM_OTUsg <- plot_heatmap(mergedT_RA_all, low="#94D2D9", high="#233A39",na.value="white", taxa.order = "Species",taxa.label="Species",
                         sample.order = "Habitat")+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
HM_OTUsg

            
ggsave(filename = "2Heatmap_OTU_RA_habitat.pdf", plot = HM_OTUsg,
       scale = 1,
       dpi = 300)

HM_Samples2 <- plot_heatmap(phylo2_RA, "PCoA", "unifrac",low="#94D2D9", high="#233A39",na.value="white" )+
  theme_pubr()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
HM_Samples2

ggsave(filename = "HM_Samples.pdf", plot = HM_Samples,
       scale = 1,
       dpi = 300)

heatdata <- tryaa2 %>% select(1,11)
heatdata2 <- tibble::rownames_to_column(heatdata, "row_names") # Apply rownames_to_column
heatdata2

heatmap_data <- veganotu(Tetracladium_master_RA)
mheatdata <- data.matrix(heatmap_data, rownames.force = NA)

library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "Blues"))(25)
li2 <- append(coul,'#FFFFFF', after = 0)

heatmap(data_otu, col = li2, Colv = NA, Rowv = NA, scale="column", labCol = F)
heat_legend(data_otu, cols =li2)




#plotly
plot_ly(x=colnames(data_otu), y=rownames(data_otu), z = data_otu, type = "heatmap") %>%
  layout(margin = list(l=120))



#geo spacial analyses

UK <- map_data(map = "world", region = "UK")

MAP_PH <-map_phyloseq(phylo2_RA,region = "UK", projection="mercator", color = "pH", size="Abundance")+
  scale_color_gradientn(colours = c("#ff1b6b", "#45caff"))

MAP_BH <-map_phyloseq(phylo_RA,region = "UK", projection="mercator", color = "Broad.habitat", size="Abundance")

MAP_AC <-map_phyloseq(phylo2_RA,region = "UK", projection="mercator", color = "Habitat", size="Abundance")

MAP_MC <-map_phyloseq(phylo2_RA,region = "UK", projection="mercator", color = "Moisture", size="Abundance")+
  scale_color_gradientn(colours = c("#ffcb6b", "#3d8bff"))
  

ggsave(filename = "MAP_MC.pdf", plot = MAP_MC,
       scale = 1,
       dpi = 300)




 #model with richnes
###########################MODELING
#grab data (has to be data frame)
datall <- data.frame(sample_data(phylo2))
observed_D_samples <- estimate_richness(phylo2, split = TRUE, measures = "Observed")
#observed_D_samples <- filter(observed_D_samples, Observed > 0)
tryaa<-cbind(observed_D_samples,datall)
#tryaa<-na.omit(tryaa)


#normalise continous variables
#define Min-Max normalization function
min_max_norm <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

tryaa2 <- tryaa[, c(1,6,7,8,9,10,11,18,44,45,17)]
tryaa2
tryaa2 <- na.omit(tryaa2)

#apply Min-Max normalization
tryaa2_norm <- as.data.frame(lapply(tryaa2[1:10], min_max_norm))

#add cat variable back
tryaa2_norm$Habitat <- tryaa2$Habitat

#selected cont variables from before
vars



#AGG.CLASS.DESCRIPTION as random variable, fixed variables based on RDA

mod3R <- lmer(Observed ~  pH  + (1|Habitat), data = tryaa2_norm)
tab_model(mod3R, show.est=TRUE)
plot_mod3R <- plot_model(mod3R, sort.est=FALSE, type = "std", show.values=TRUE, value.size=2, value.offset=0.4)
plot_mod3R


performance_accuracy(
  mod3R,
  method = "boot",
  k = 5,
  n = 1000,
  verbose = TRUE)

CMplotR <- check_model(
  mod3R,
  dot_size = 2,
  line_size = 0.8,
  panel = TRUE,
  check = "all",
  alpha = 0.2,
  dot_alpha = 0.8,
  colors = c("#3aaf85", "#1b6ca8", "#cd201f"),
  theme = "see::theme_lucid",
  detrend = FALSE,
  verbose = TRUE)

# Use the effects package --> effect function. term= the fixed effect you want to get data on, mod= name of your model.

effectsR <- effects::effect(term= "pH", mod= mod3R)
summary(effectsR) #output of what the values are
x_PH_R <- as.data.frame(effectsR)

pH_plot_R <- ggplot() + 
  geom_point(data=tryaa2_norm, aes(pH, Observed), colour= "#C46854", fill ="#C46854",size = 3, alpha = 0.3 ) + 
  geom_point(data=x_PH_R, aes(x=pH, y=fit), color="#C46854") +
  geom_line(data=x_PH_R, aes(x=pH, y=fit), color="#C46854") +
  geom_ribbon(data= x_PH_R, aes(x=pH, ymin=lower, ymax=upper), alpha= 0.3, fill="#C46854") +
  labs(x="pH", y="Tetracladium OTU Richness") +
  theme_pubr()

pH_plot_R

ggsave("pH_plot_R.pdf", pH_plot_R, dpi = 300)

model_coefsR <- coef(mod3R)$ Habitat %>% 
  rename(Intercept = `(Intercept)`, Slope = pH) %>% 
  rownames_to_column("Habitat")

sleep_groups_raniR <- left_join(tryaa2_norm, model_coefsR, by = "Habitat")

model_coef_plotR <- ggplot(data = sleep_groups_raniR, 
                          mapping = aes(x = pH, 
                                        y = Observed, 
                                        colour = Habitat)) +
  geom_point(na.rm = T, alpha = 0.5) +
  geom_abline(aes(intercept = Intercept, 
                  slope = Slope,
                  colour = Habitat),size = 1.5) +
  theme_pubr()

###ANOSIM
unique(get_variable(phylo_RA, "Habitat"))
group = get_variable(phylo_RA, "Habitat")
ano = anosim(phyloseq::distance(phylo_RA, "bray"), group, permutations = 999)
# p value
ano$signif
#R value
ano$statistic
#Matrix of pairwise results

#Run this function to line 95 as it copes with hyphens in names (Adriana data). 
make_anosim_table<-function(otus,grouping_factor,meta=NULL){
  require(vegan)
  if(is(otus)=="phyloseq"){
    groupings<-as.character(sample_data(otus)[[grouping_factor]])
    if(class(groupings)=="character"){
      ind_groups<-unique(groupings)
    } else {
      ind_groups<-levels(sample_data(otus)[[grouping_factor]])
    }
    
    
    names(ind_groups)<-ind_groups
    ## Call it
    test<-lapply(ind_groups, function(x) lapply(ind_groups, function(y) do_anosim(otus,groupings=groupings,grp1=x,grp2=y,grouping_factor=grouping_factor)))
  } else {
    if(is.null(meta)){
      stop("You must specify a metadata table.")
    }
    ord =match(colnames(otus),rownames(metadata)) #Get matching order
    metadata = metadata[ord, ] #rearrange order
    
    groupings<-meta[,grouping_factor]
    ind_groups<-levels(droplevels(groupings))
    names(ind_groups)<-ind_groups
    
    ## Call it
    test<-lapply(ind_groups, function(x) lapply(ind_groups, function(y) do_anosim(otus,groupings=groupings,grp1=x,grp2=y,grouping_factor=grouping_factor,meta=meta)))
  } 
  # Format results
  R<-matrix()
  P<-matrix()
  pvals<-list()
  rvals<-list()
  for(i in ind_groups){
    pvals<-c(pvals,sapply(ind_groups,function(x) test[[i]][x][[1]]["P"]))
    rvals<-c(rvals,sapply(ind_groups,function(x) test[[i]][x][[1]]["R"]))
  }
  pvals<-matrix(unlist(pvals),ncol=length(ind_groups))
  rownames(pvals)<-ind_groups
  colnames(pvals)<-ind_groups
  rvals<-matrix(unlist(rvals),ncol=length(ind_groups))
  rownames(rvals)<-ind_groups
  colnames(rvals)<-ind_groups
  return(list("all_res"=test,"p"=pvals,"r"=rvals))
}

anosim_table <- make_anosim_table(phylo_RA,"Habitat")

anosim_table

#To plot a grid. Change pvals to TRUE to see actual values on chart. If get "Panic" then this means that one or more of the R values are <0. This doesn't matter. The plot will still draw.
plotANO <- plot_anosim_table(anosim_table, pvals = T) +
  scale_fill_gradient2(low="#CB8FBF",
                         high="#4F858A", mid = "white", limits = c(-0.2,1)) 

plotANO

#Original P values
anosim_table$p

#Remove one side of the tables.
anosim_table$p[lower.tri(anosim_table$p,diag=TRUE)] <- NA

anosim_table$p

##SEM
#data RA
mod_data_norm

model <- ' x ~ pH + Habitat +Longitude +Moisture +Ellenburg.N'
fit <- sem(model, data=mod_data_norm)
summary(fit)
semPaths(fit, 
         what = "std",
         whatLabels = "std",
        layout = "tree",
        intercepts = T,
        residuals = F)


Model_pSEM <- psem(
  lmer(x ~ pH +Longitude +Moisture +P+ Latitude+Organic.C  +C+ N+Ellenburg.N +(1|Habitat), data =mod_data_norm),
  lmer(Longitude ~pH +Moisture +P+Organic.C  +C+ N+Ellenburg.N +(1|Habitat), data =mod_data_norm ))


PSEM <- plot(
  Model_pSEM,
  return = F,
  node_attrs = data.frame(shape = "rectangle", style = "rounded", fixedsize = F, color = "#2F4F4F", fillcolor = "white", fontsize = 8, penwidth = 2),
  edge_attrs = data.frame(style = "solid", color = "black", arrowhead = "none"),
  ns_dashed = T,
  alpha = 0.05,
  show = "std",
  digits = 3,
  add_edge_label_spaces = T,
  title = "RA PSEM")

PSEM
 summary(Model_pSEM, test.statistic = "T")
#data richness
tryaa2_norm

Model_pSEM_D <- psem(
  lmer(Observed ~ pH +Longitude +Moisture +P+ Latitude+Organic.C  +C+ N+Ellenburg.N +(1|Habitat), data =tryaa2_norm),
  lmer (Longitude ~pH +Moisture +P+Organic.C  +C+ N+Ellenburg.N +(1|Habitat), data =tryaa2_norm ))


PSEM_D <- plot(
  Model_pSEM_D,
  return = F,
  node_attrs = data.frame(shape = "rectangle", style = "rounded", fixedsize = F, color = "#2F4F4F", fillcolor = "white", fontsize = 8, penwidth = 2),
  edge_attrs = data.frame(style = "solid", color = "black", arrowhead = "none"),
  ns_dashed = T,
  alpha = 0.05,
  show = "std",
  digits = 3,
  add_edge_label_spaces = T,
  title = "RA PSEM")

PSEM_D
summary(Model_pSEM_D, test.statistic = "T")
sink("lm_D.txt")
print(summary(Model_pSEM_D, test.statistic = "T"))
sink()  


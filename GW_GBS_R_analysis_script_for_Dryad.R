# Introduction ----
# GW_GBS_R_analysis_script_for_Dryad.R
#
# This script contains the R processing code used in this manuscript:
#
# Irwin DE, Alcaide M, Delmore KE, Irwin JH, Owens GL. Recurrent selection explains parallel evolution of genomic regions of high relative but low absolute differentiation in a ring species. In revision, Molecular Ecology.
#
# also posted on bioRxiv:
#
# Irwin DE, Alcaide M, Delmore KE, Irwin JH, Owens GL. 2016. Recurrent selection explains parallel evolution of genomic regions of high relative but low absolute differentiation in greenish warblers. bioRxiv, doi: http://dx.doi.org/10.1101/041467
#
# The functions are in a file called genomics_R_functions.R
# This script will give general info and use those functions.
# Please cite the above paper if you use these scripts.

# Set directories and functions ----

# Set directory where files stored:
setwd("~/Dropbox/Darren's current work/genomics training 2014/GW_genomics_Darren/")

# Load functions:
source ("~/Dropbox/Darren's current work/YRW genomics/YRWA_genomics/genomics_R_functions.R")

# Setup for main processing ----
# Note that the main processing of big chromosomes can take a long time (e.g. overnight);
# For the actual paper I ran the chromosomes in batches as defined by "chromosomes.to.analyze" below.
# Here I set it just to analyze chromosome 1A, which will reproduce Figure 3.
# This takes some time; If you want to skip to producing the other figures,
# skip below down to heading "GENOME-WIDE plots"

# choose the chromosomes to analyze in this run:
chromosomes.to.analyze <- c("1A")
# to process other chromosomes, choose from among these:
# chromosomes.to.analyze <- c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","Z")

# Options to calculate the per-site and windowed stats, or to load already calculated stats:
calculate_or_load_stats <- 1  # 1) calculate site stats;  
# 2) load previously calculated per-site stats; 
# 3) skip per-site stats (and load windowed data from file--set "load.rolling.means" to T below)
saveSiteInfo <- F   # If TRUE, will save a file for per-site stats
saveWindowedStats <- F   # If TRUE, will save a file for per-window stats
load.rolling.means <- T   # If TRUE, will load rolling mean data (rather than calculate)

set <- 1  # 1 is 45 samples troch_vir_plumb;  
          # 2 is nine_taxa one.per.taxon (the "Outgroup analysis"--see that heading below)

if (set == 1) {
  # choose path and filename for the 012NA files
  base.name <- "GW_Lane5plusLiz.GWref.012NA_files/GW_Lane5plusLiz.GWref.45samples"
  filename.text.middle <- ".max2allele_noindel.maxmiss"
  # indicate percent threshold for missing genotypes for each SNP--
  # this was set by earlier filtering, and is just a record-keeper for the filenames:
  missing.genotype.threshold <- 60 
  filename.text.end <-".MQ20.lowHet.tab"
  # choose a tag name for this analysis
  tag.name <- ".troch_vir_plumb."    # to avoid writing over existing files, change name, e.g.: ".TEST.troch_vir_plumb."
  # indicate name of metadata file, a text file with these column headings:
  # ID	location	group	Fst_group	plot_order
  metadata.file <- "GW_Lane5_plus_Liz.45samples.Fst_groups.txt"
  # specify number of individuals in file (good for error checking)
  num.individuals <- 45
  # specify window size (number of bp with info) and step size
  window_size <- 5000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("troch", "vir", "plumb")
  group_count <- length(groups)
  # specify groups for plotting, and their colors
  groups.to.plot.WC84_Fst <- c("troch_vir", "troch_plumb", "vir_plumb")
  group.colors.WC84_Fst <- c("green3", "orange", "purple")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
  group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
  groups.to.plot.pi <- c("vir", "plumb", "troch") 
  group.colors.pi <- c("blue", "red", "gold1")
  groups.to.plot.PCA <- c("vir","vir_S","nit","lud","troch_MN","troch_LN","troch_EM","obs","plumb")
  group.colors.PCA <- c("blue","turquoise1","grey","green","greenyellow","yellow","gold","orange","red")
} else if (set == 2) {   # This is the set for the "Outgroup analysis"--see that heading below
  # choose path and filename for the 012NA files
  base.name <- "GW_Lane5plusLiz.GWref.012NA_files/GW_Lane5plusLiz.GWref.nine_taxa"
  filename.text.middle <- ".max2allele_noindel.maxmiss"
  # indicate percent threshold for missing genotypes for each SNP--
  # this was set by earlier filtering, and is just a record-keeper for the filenames:
  missing.genotype.threshold <- 60 
  filename.text.end <-".MQ20.lowHet.tab"
  # choose a tag name for this analysis
  tag.name <- ".one_per_taxon."  # to avoid writing over existing files, change name, e.g.: ".TEST.one_per_taxon."
  # indicate name of metadata file, a text file with these column headings:
  # ID	location	group	Fst_group	plot_order
  metadata.file <- "GW_Lane5_plus_Liz.nine_taxa.Fst_groups.txt"
  # specify number of individuals in file (good for error checking)
  num.individuals <- 9
  # specify window size (number of bp with info) and step size
  window_size <- 10000
  step_size <- window_size  # could change if wanting overlapping windows
  # specify groups for calculation of statistics
  groups <- c("nit", "vir", "troch", "obs", "plumb", "inor", "hume", "fusc", "burk")
  group.colors <- c("purple", "blue", "yellow", "orange", "red", "grey20", "grey80", "salmon", "turquoise1")
  groups.to.plot.WC84_Fst <- c("vir_plumb", "inor_hume")
  group.colors.WC84_Fst <- c("purple", "grey50")
  groups.to.plot.Dxy <- groups.to.plot.WC84_Fst 
  group.colors.Dxy <- group.colors.WC84_Fst 
  groups.to.plot.pi <- c("vir", "plumb", "inor", "hume")  #, "fusc", "burk")
  group.colors.pi <- c("blue", "red", "grey", "orange")  #, "salmon", "turquoise1")
  groups.to.plot.PCA <- groups
  group.colors.PCA <- group.colors
} else {
  stop("No set chosen")
}


# Option to focus on a region of chromosome ----
# (not used in paper)
focus.region <- F  # choose T for a subset of the chromosome, F for the whole thing)
if (focus.region==T){
  position.min <-  1500000
  position.max <- 1750000
}



# MAIN LOOP ----
# -----
for (i in 1:length(chromosomes.to.analyze)) {
  chr <- chromosomes.to.analyze[i] 
  
  # Get chr data ---- 
  # read in position data for this chromosome
  position.file.name <- paste0(base.name, ".chr", chr, filename.text.middle, missing.genotype.threshold, filename.text.end, ".012.pos")
  pos.whole.chr <- read.table(position.file.name, col.names = c("chrom", "position"))
  # read in genotype data for this chromosome
  column_names <- c("null", paste("c", pos.whole.chr$chrom, pos.whole.chr$position, sep="."))
  genotype.file.name <- paste0(base.name, ".chr", chr, filename.text.middle, missing.genotype.threshold, filename.text.end, ".012NA")
  geno<-read.table(genotype.file.name, nrows = num.individuals, colClasses = "integer", col.names = column_names)
  loci_count <- length(geno[1,]) -1   # because the first column is not a SNP (just a count from zero)
  # read in individual names for this chromosome dataset
  individuals.file.name <- paste0(base.name, ".chr", chr, filename.text.middle, missing.genotype.threshold, filename.text.end, ".012.indv")
  ind<-read.table(individuals.file.name)
  
  # Get metadata ----
  locations <- read.table(paste0(metadata.file), header=TRUE)
  num_loc_cols <- length(locations[1,])
  ind_with_locations <- cbind(ind,locations) 
  print(ind_with_locations)    
  print("check first two columns to make sure the same")
  # combine metadata with genotype data
  combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])
  # If need to filter out individuals, based on low read number:
  filter <- F
  if (filter==T){
    # Specify individuals to filter out:
    combo.NApass.whole.chr <- combo[c(-20,-163),]
  } else if (1==1) {
    combo.NApass.whole.chr <- combo
  }
  
  # Get region text ----
  if (focus.region==F) {
    position.min <- 1
    position.max <- pos.whole.chr$position[length(pos.whole.chr$position)]
    pos <- pos.whole.chr
    combo.NApass <- combo.NApass.whole.chr
    region.text <- paste0("Chr",chr,"_whole")
  } else if (focus.region==T) {
    selection <- pos.whole.chr$position >= position.min & pos.whole.chr$position <= position.max
    pos <- pos.whole.chr[selection,]
    selection <- c(rep(TRUE, times=num_loc_cols), selection)
    combo.NApass <- combo.NApass.whole.chr[, selection]
    region.text <- paste0("Chr",chr,"_from_",position.min,"_to_",position.max)
  }
  
  # Make site stats ---- 
  if (calculate_or_load_stats==1) {
    # Calculate allele freqs and sample sizes (use column Fst_group)
    temp.list <- getFreqsAndSampleSizes(combo.NApass, num_loc_cols, groups)
    freqs <- temp.list$freqs
    sample_size <- temp.list$sample_size
    rm(temp.list)
    print("Calculated population allele frequencies and sample sizes")
    
    # calculate nucleotide diversity (pi) at each site for each population
    site_pi <- getSitePi(freqs) 
    print("Calculated population pi values")
    
    # calculate rownames for pairwise comparisons, for use in Dxy and Fst matrices:
    # rownames <- getPairwiseNames(groups)   # NOT NEEDED HERE since called from getDxy
    
    # calculate Dxy at each site, between pairs of groups
    Dxy <- getDxy(freqs, groups)
    print("Calculated Dxy values")
    
    # calculate Fst (and numerator and denominator) for each site, 
    # between pairs of groups (so pops (r) is 2), 
    # using the Weir&Cockerham 1984 approach to correct for sample size and number of pops
    temp.list <- getWC84Fst(freqs, sample_size, groups, among=TRUE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
    WC84_Fst <- temp.list$WC84_Fst
    WC84_Fst_numerator <- temp.list$WC84_Fst_numerator
    WC84_Fst_denominator <- temp.list$WC84_Fst_denominator
    rm(temp.list)
    print("Calculated WC84_Fst values")
    
    if (saveSiteInfo == TRUE) {  # save the per-site stats, if chosen to in Intro section
      save(pos, freqs, sample_size, site_pi, WC84_Fst, WC84_Fst_numerator, WC84_Fst_denominator, Dxy, file=paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
      print("Saved summary site stats")
      WC84_Fst_among <- WC84_Fst[rownames(WC84_Fst)=="Fst_among",]
      print(paste0(length(WC84_Fst_among)," markers in total (",sum(!is.na(WC84_Fst_among))," variable and ",sum(is.na(WC84_Fst_among)), " invariant)"))
    } else print("Site stats not saved")
    
  } else if (calculate_or_load_stats==2) {
    load(paste0(base.name,tag.name,region.text,"_SiteStats_from_R.R"))
    print("Loaded saved summary stats")
  }
  
  
  # Make windowed stats ---- 
  if (calculate_or_load_stats==1 | calculate_or_load_stats==2) {
    
    # calculate windowed pi, in whole windows starting on left side of chromosome
    temp.list <- getWindowedPi(site_pi, pos, window_size, step_size)
    rolling.mean.pos.pi <- temp.list$rolling.mean.pos.pi
    rolling.mean.pi <- temp.list$rolling.mean.pi
    rm(temp.list)
    
    # calculate windowed Dxy, in whole windows starting on left side of chromosome
    temp.list <- getWindowedDxy(Dxy, pos, window_size, step_size)
    rolling.mean.pos.Dxy <- temp.list$rolling.mean.pos.Dxy
    rolling.mean.Dxy <- temp.list$rolling.mean.Dxy
    rm(temp.list)
    
    # calculate windowed Fst according to according to Weir&Cockerham1984 
    # (with sample size and pop number correction),
    # calculated as windowed numerator over windowed denominator.
    temp.list <- getWindowedWC84_Fst(WC84_Fst_numerator, WC84_Fst_denominator, pos, window_size, step_size)
    rolling.mean.pos.WC84_Fst <- temp.list$rolling.mean.pos.WC84_Fst
    rolling.mean.WC84_Fst <- temp.list$rolling.mean.WC84_Fst
    rm(temp.list)
  }
  
  # Save or load ----
  # save the rolling mean data, if chosen in Intro section
  if (saveWindowedStats == TRUE) {
    save(rolling.mean.pos.pi, rolling.mean.pi, rolling.mean.pos.Dxy, rolling.mean.Dxy, rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst, file=paste0(base.name,tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R"))
    print("Saved rolling mean stats")
  }
  
  # load the rolling mean data, if chosen in Intro section:
  if (load.rolling.means == TRUE) {
    load(paste0(base.name,tag.name,region.text,"_window",window_size,"_WindowStats_from_R.R"))
  }
  
  # Make plot for chr ----
  # make ggplots for quick inspection of rolling mean results (not shown in paper)
  makeRollingMeanPlots(rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst,
                       rolling.mean.pos.Dxy, rolling.mean.Dxy,
                       rolling.mean.pos.pi, rolling.mean.pi, 
                       group.colors.pi, groups.to.plot.pi,
                       group.colors.Dxy, groups.to.plot.Dxy,
                       group.colors.WC84_Fst, groups.to.plot.WC84_Fst,
                       region.text)
  
}   
# End main loop ----


# FIGURE 3 ----
# Detailed chr plot ----
# Plots of Fst, Dxy, and pi for three groups across chromosome
# Note there is a small amount of vertical jitter in point locations
# (depends on running this file from the top)
makeDetailedChrPlots(pos, 
                     WC84_Fst, rolling.mean.pos.WC84_Fst, rolling.mean.WC84_Fst,
                     Dxy, rolling.mean.pos.Dxy, rolling.mean.Dxy,
                     site_pi, rolling.mean.pos.pi, rolling.mean.pi, 
                     group.colors.pi, groups.to.plot.pi,
                     group.colors.Dxy, groups.to.plot.Dxy,
                     group.colors.WC84_Fst, groups.to.plot.WC84_Fst,
                     region.text)


# GENOME-WIDE plots -------------------------------------------

# FIGURE 4 ----
# Read windowed data files for each chromosome, and plot Fst, Dxy, and pi for each chromosome in a single window.
base.name <- "~/Dropbox/Darren's current work/genomics training 2014/GW_genomics_Darren/GW_Lane5plusLiz.GWref.012NA_files/GW_Lane5plusLiz.GWref.45samples"
tag.name <- ".troch_vir_plumb."
window_size <- 5000  # this must be the same as was used to generate the windowed data
groups.to.plot.WC84_Fst <- c("troch_vir", "troch_plumb", "vir_plumb")
group.colors.WC84_Fst <- c("green3", "orange", "purple")
groups.to.plot.Dxy <- groups.to.plot.WC84_Fst   # or specify differences if necessary
group.colors.Dxy <- group.colors.WC84_Fst  # or specify differences if necessary
groups.to.plot.pi <- c("vir", "plumb", "troch") 
group.colors.pi <- c("blue", "red", "gold1")
chromosomes.to.plot <- c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","Z")
max_Dxy_axis <- 0.008  # height of Dxy axis
max_pi_axis <- 0.008  # height of pi axis
transparency <- 0.2  # transparency of the colored area under the lines
line_transparency <- 0.8  # transparency of the colored lines
# plot Genome-wide plot of windowed Fst, Dxy, pi
plotGenomeFstDxyPi(base.name, tag.name, window_size, chromosomes.to.plot, 
                   max_Dxy_axis, max_pi_axis, transparency, line_transparency,
                   groups.to.plot.WC84_Fst, group.colors.WC84_Fst,
                   groups.to.plot.Dxy, group.colors.Dxy,
                   groups.to.plot.pi, group.colors.pi)


# COMPILE genome-wide info  -----------------
# compile windowed WC84_Fst, Dxy, pi for a bunch of chromosomes
# NOTE: This is needed for most remaining figures in the paper.
# compile autosomal-genome-wide info:
chromosomes.to.combine <- c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28")
autosome.genome.rolling.stats <- compileWindowedStats(base.name, tag.name, chromosomes.to.combine, window.size)


# Use the compiled info for the following:

# Fst boundary analysis ----
# make a histogram of Fst for one group
# and compare pi and Dxy above and below an Fst cutoff
group1 <- "troch"         
group2 <- "vir" 
Fst_boundary <- 0.6
Fst_boundary_stats_0.6 <- Fst_hist_and_boundary_stats(group1, group2, 
                                                  autosome.genome.rolling.stats, 
                                                  Fst_boundary)
# Some statistics presented in Discussion of paper:
Fst_boundary_stats_0.6$mean_pi_below_Fst_boundary / Fst_boundary_stats_0.6$mean_pi_above_Fst_boundary
# 5.158131
Fst_boundary <- 0.9
Fst_boundary_stats_0.9 <- Fst_hist_and_boundary_stats(group1, group2, 
                                                  autosome.genome.rolling.stats, 
                                                  Fst_boundary)
Fst_boundary_stats_0.6$mean_pi_below_Fst_boundary / Fst_boundary_stats_0.9$mean_pi_above_Fst_boundary
# 14.83475


# FIGURE 5 ----
# Fst 3 plots ----
# make scatterplots of Fst vs. Fst for three comparisons of group pairs;
# each of the "comp_" objects should have two group names and a color name;
# returns the statistical correlation tests, using the cor.method as specified
comp1 <- c("troch", "vir", "green3")  
comp2 <- c("troch", "plumb", "orange")
comp3 <- c("vir", "plumb", "purple")
cor.method <- "spearman"
cor.tests <- plots3Fst_Fst(comp1, comp2, comp3,
                           autosome.genome.rolling.stats, cor.method)
cor.tests  # prints correlation tests for 1x2, 1x3, 2x3



# FIGURE 6 ----
# Fst_Dxy 3 plots ----
# make scatterplots of Fst vs. Dxy for three group comparisons
# each of the "comp_" objects should be have two group names and a color name, e.g.:
# returns the statistical correlation tests, using the cor.method as specified
comp1 <- c("troch", "vir", "green3")  
comp2 <- c("troch", "plumb", "orange")
comp3 <- c("vir", "plumb", "purple")
cor.method <- "spearman"
cor.tests <- plots3Fst_Dxy(comp1, comp2, comp3,
                           autosome.genome.rolling.stats, cor.method)
cor.tests  # prints correlation tests for comp1, 2, 3


# FIGURE 7 ----
# Dxy_MeanPi plot ----
# Scatterplot of Dxy vs. mean_pi using compiled genome-wide info .
# Color the points based on Fst.
group1 <- "vir"
group2 <- "plumb"
max.x <- 0.008
max.y <- 0.006  
cor.method <- "spearman"
test <- plotDxy_MeanPi(group1, group2, cor.method,
                       autosome.genome.rolling.stats,
                       max.x, max.y, color_bar=T)
test   # print test of correlation


# FIGURE 8A ----
# pi_pi 3 plots ----
# Make 3 scatterplots of pi in a line
groups <- c("troch", "vir", "plumb")  
colors <- c("gold1", "blue", "red")
cor.method <- "spearman"
cor.tests <- plots3pi_pi(groups, colors, cor.method, 
                         autosome.genome.rolling.stats)
cor.tests  # prints correlation tests for 1x2, 1x3, 2x3

# FIGURE 8B ----
# Pi/DxyMax 3 plots ----
# Make 3 scatterplots of pi/Dxy_max in a line:
# Dxy_max is the highest Dxy value for all three pairwise comparisons
groups <- c("troch", "vir", "plumb")  
colors <- c("gold1", "blue", "red")
cor.method <- "spearman"
cor.tests <- plots3PiOverDxyMax_corr(groups, colors, cor.method,
                                     autosome.genome.rolling.stats)
cor.tests  # prints correlation tests for 1x2, 1x3, 2x3


# Compile variant numbers ----
# Reads files containing per-bp info for each chromosome;
# requires that "Fst_among" is in the WC84_Fst matrix;
# returns a summary of variant and invariant bp per chromosome
base.name <- "~/Dropbox/Darren's current work/genomics training 2014/GW_genomics_Darren/GW_Lane5plusLiz.GWref.012NA_files/GW_Lane5plusLiz.GWref.45samples"
tag.name <- ".troch_vir_plumb."
chromosomes.to.summarize <- c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28","Z")
num_loci_summary <- compileVariantNumbers(base.name, tag.name, chromosomes.to.summarize) 
num_loci_summary   # print table of variant and invariant bp per chromosome
# loci  SNPs invariant
# 1  1073594 33166   1040428
# 1A  738167 18145    720022
# 1B   20850   460     20390
# 2  1282248 39394   1242854
# 3  1153399 34745   1118654
# 4   633419 18157    615262
# 4A  408083 13061    395022
# 5   755285 22966    732319
# 6   511309 16406    494903
# 7   509346 15304    494042
# 8   390698 12175    378523
# 9   480076 15118    464958
# 10  361535 10787    350748
# 11  367564 10567    356997
# 12  450185 12873    437312
# 13  365141 10022    355119
# 14  380334 10482    369852
# 15  342955  9568    333387
# 17  337854  8706    329148
# 18  256295  6852    249443
# 19  306760  9280    297480
# 20  397720 11045    386675
# 21  117309  2928    114381
# 22   70718  1535     69183
# 23  179519  4544    174975
# 24  187501  5118    182383
# 25   26574   582     25992
# 26  145916  3405    142511
# 27  112795  3074    109721
# 28  128320  2986    125334
# Z   522481 11388    511093
sapply(num_loci_summary, FUN=sum)   # print totals of loci, SNPs, and invariant sites
#      loci      SNPs invariant 
#  13013950    374839  12639111 


# compare Z & autosome ----
# script to compare autosomal and Z-chromosome differentiation statistics
base.name <- "~/Dropbox/Darren's current work/genomics training 2014/GW_genomics_Darren/GW_Lane5plusLiz.GWref.012NA_files/GW_Lane5plusLiz.GWref.45samples"
tag.name <- ".troch_vir_plumb."
window_size <- 5000
# compile autosomal-genome-wide info on WC84_Fst, Dxy, pi:
chromosomes.to.combine <- c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28")
autosome.genome.rolling.stats <- compileWindowedStats(base.name, tag.name, chromosomes.to.combine, window.size)

# compile Z-chromosome info on WC84_Fst, Dxy, pi:
chromosomes.to.combine <- c("Z")
Z.rolling.stats <- compileWindowedStats(base.name, tag.name, chromosomes.to.combine, window.size)

#Stats for 45-sample analysis:
rowMeans(autosome.genome.rolling.stats$pi)
#       troch         vir       plumb 
# 0.003007980 0.002072852 0.002875735 
rowMeans(Z.rolling.stats$pi)
#       troch         vir       plumb 
# 0.002026023 0.001175741 0.001818695 
rowMeans(autosome.genome.rolling.stats$Dxy)
# troch_vir troch_plumb   vir_plumb 
# 0.003859455 0.003737568 0.003822110 
rowMeans(Z.rolling.stats$Dxy)
# troch_vir troch_plumb   vir_plumb 
# 0.003218566 0.002898285 0.003298915 
rowMeans(autosome.genome.rolling.stats$WC84_Fst)
# troch_vir troch_plumb   vir_plumb   Fst_among 
# 0.3245106   0.1999270   0.3335126   0.2830489 
rowMeans(Z.rolling.stats$WC84_Fst)
# troch_vir troch_plumb   vir_plumb   Fst_among 
# 0.4772868   0.3122621   0.5183402   0.4308442 

# t-test comparing Fst across autosomes and across Z
row.choice.1 <- which(rownames(autosome.genome.rolling.stats$WC84_Fst) == "vir_plumb")
row.choice.2 <- which(rownames(Z.rolling.stats$WC84_Fst) == "vir_plumb")
t.test(autosome.genome.rolling.stats$WC84_Fst[row.choice.1,], Z.rolling.stats$WC84_Fst[row.choice.2,])
# Welch Two Sample t-test
# data:  WC84_Fst.vector and Z.WC84_Fst.vector
# t = -12.0768, df = 115.944, p-value < 2.2e-16
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
# -0.2151399 -0.1545153
# sample estimates:
# mean of x mean of y 
# 0.3335126 0.5183402 

# FIGURE 10
# MeanPi_Dxy autosome vs. Z plot ----
# uses compiled data above
group1 <- "vir"
group2 <- "plumb"
tests <- plotMeanPi_Dxy.autosome_Z(group1, group2,
                                   autosome.genome.rolling.stats, Z.rolling.stats) 
tests  # show t-tests comparing autosomes and Z in both pi and Dxy



# FIGURE S1 ----
# PCA whole-genome ----
# Load vcf file containing only variable sites throughout genome;
# construct a PCA based on all sites passing an Fst threshold between the "groups" below;
# and all individuals in "groups.to.plot.PCA" according to colors in "group.colors.PCA"
groups <- c("troch", "vir", "plumb")  # for purpose of calculating pairwise Fst and Fst_group (to determine SNPs)
groups.to.plot.PCA <- c("vir","vir_S","nit","lud","troch_MN","troch_LN","troch_EM","obs","plumb")
group.colors.PCA <- c("blue","turquoise1","grey","green","greenyellow","yellow","gold","orange","red")
base.file.name <- "~/Dropbox/Darren's current work/genomics training 2014/GW_genomics_Darren/GW_Lane5plusLiz.GWref.012NA_files/GW_Lane5plusLiz.GWref.genotypes.SNPs_only.whole_genome.GW_only.max2allele_noindel.maxmiss60.MQ20.lowHet.tab"
pos <- read.table(paste0(base.file.name, ".012.pos"))
column_names <- c("null", paste("c", pos$V1, pos$V2, sep="."))
geno <- read.table(paste0(base.file.name, ".012NA"), colClasses = "integer", col.names = column_names)
SNPnum <- length(geno[1,]) -1   # because the first column is not a SNP (just a count from zero)
ind <- read.table(paste0(base.file.name, ".012.indv"))
locations <- read.table("~/Dropbox/Darren's current work/genomics training 2014/GW_genomics_Darren/GW_Lane5_plus_Liz.GW_only.adults_only.Fst_groups.txt", header=TRUE)
num_loc_cols <- length(locations[1,])
ind_with_locations <- cbind(ind,locations)
combo <- cbind(ind_with_locations[,2:(num_loc_cols+1)],geno[,2:length(geno[1,])])
# determine number of missing SNPs per bird, and filter out those with more than X% missing SNPs
X <- 25   # this is the percentage threshold
threshold_NA <- SNPnum * X/100
numNAs <- rowSums(is.na(combo))
selection <- (numNAs < threshold_NA)
combo.NApass.all <- combo[selection,]
ind[which(selection==F),]  # prints out the individuals being left out.
# [1] GW_Lane5_AA8       GW_Lane5_DA6       GW_Liz_GBS_Liz5118 GW_Liz_GBS_Liz5195
# option to filter out all but selected chromosome (or set of them):
choose.chrom <- FALSE
if (choose.chrom == TRUE) {
  chrom <- "5"
  loci.selection <- (pos$V1 == chrom)
  loci.selection <- c(TRUE, TRUE, TRUE, loci.selection)  # add placeholders for info columns
  combo.NApass <- combo.NApass.all[,loci.selection]
  region.text <- paste0("chr", chrom)	
}	else {
  region.text <- "whole_genome"
  combo.NApass <- combo.NApass.all
}
# Calculate allele freqs and sample sizes (use column Fst_group)
temp.list <- getFreqsAndSampleSizes(combo.NApass, num_loc_cols, groups)
freqs <- temp.list$freqs
sample_size <- temp.list$sample_size
rm(temp.list)
# calculate WC84_Fst 
temp.list <- getWC84Fst(freqs, sample_size, groups, among=TRUE)  # set among to FALSE if no among Fst wanted (some things won't work without it)
WC84_Fst <- temp.list$WC84_Fst
rm(temp.list)
# make the figure:
Fst.filter <- F   # option to filter to high-Fst markers only, using cutoff below
WC84_Fst.cutoff <- 0.75  # has no effect if Fst.filter is FALSE
# choose whether to filter by Fst between pair of populations, or by Fst_among (as defined above)
groups.to.compare <- "Fst_among"
axes <- 3
PCA_results <- plotPCA(Fst.filter, WC84_Fst.cutoff, groups.to.compare, WC84_Fst, combo.NApass, num_loc_cols, region.text,
                       groups.to.plot.PCA, group.colors.PCA, axes)
# Re-plot with axis 2 flipped to make more parallel to geography (it is legit to flip PCA axes):
quartz(title="PCA of whole genome",8, 7)
plot(PCA_results$scores$PC1,-1*PCA_results$scores$PC2, asp=1, pch=23, xlab="PC1", ylab="PC2", main="PCA of whole genome")
for (i in 1:length(groups.to.plot.PCA)) {
  selection <- PCA_results$data$group == groups.to.plot.PCA[i]
  points(PCA_results$scores$PC1[selection], -1*PCA_results$scores$PC2[selection], pch=23, bg = group.colors.PCA[i])
}
# The above produces FIGURE S1 as shown in the paper.

#####################################


# Outgroup Analysis ----
# Used for FIGURE 9, FIGURE S2, FIGURE S3, FIGURE S4, FIGURE S5. 
# For analysis of patterns among just nine taxa, 
# one individual in each (burk, fusc, inor, humei, nit, vir, troch, obs, plumb). 

# Setup: 
base.name <- "~/Dropbox/Darren's current work/genomics training 2014/GW_genomics_Darren/GW_Lane5plusLiz.GWref.012NA_files/GW_Lane5plusLiz.GWref.nine_taxa"
tag.name <- ".one_per_taxon."        
window_size <- 10000  # this must be the same as was used to generate the windowed data
groups <- c("nit", "vir", "troch", "obs", "plumb", "inor", "hume", "fusc", "burk")
group.colors <- c("purple", "blue", "yellow", "orange", "red", "grey20", "grey80", "salmon", "turquoise1")
groups.to.plot.WC84_Fst <- c("vir_plumb", "inor_hume")
group.colors.WC84_Fst <- c("purple", "grey50")
groups.to.plot.Dxy <- groups.to.plot.WC84_Fst 
group.colors.Dxy <- group.colors.WC84_Fst 
groups.to.plot.pi <- c("vir", "plumb", "inor", "hume")  #, "fusc", "burk")
group.colors.pi <- c("blue", "red", "grey", "orange")  #, "salmon", "turquoise1")
groups.to.plot.PCA <- groups
group.colors.PCA <- group.colors

# COMPILE autosomal genome-wide info:
chromosomes.to.combine <- c("1","1A","1B","2","3","4","4A","5","6","7","8","9","10","11","12","13","14","15","17","18","19","20","21","22","23","24","25","26","27","28")
autosome.genome.rolling.stats <- compileWindowedStats(base.name, tag.name, chromosomes.to.combine, window.size)

# FIGURE 9 ----
# Pi scatterplot matrix 4groups
group1 <- "vir"
group2 <- "plumb"
group3 <- "inor"
group4 <- "hume"
cor.method <-  "spearman"   # or "pearson" 
cor.tests <- PiMatrixPlot4(group1, group2, group3, group4, cor.method, 
                           autosome.genome.rolling.stats,
                           xlim=c(0,0.003), ylim=c(0,0.003), xaxp=c(0,0.002,1), yaxp=c(0,0.002,1))
cor.tests  # prints correlation tests for 1x2, 1x3, 1x4, 2x3, 2x4, 3x4
# The above produces the data plots for FIGURE 9; 
# additional label processing was done using Adobe Illustrator.


# FIGURE S2 ----
# Dxy scatterplot matrix
# Make a scatterplot matrix for Dxy using above compiled genome-wide info.
comp1 <- "vir_troch"			
comp2 <- "obs_plumb" 			
comp3 <- "fusc_burk" 
cor.method <- "pearson"
cor.tests <- DxyMatrixPlot(comp1, comp2, comp3, cor.method, 
                           autosome.genome.rolling.stats, cex.labels = 2)
cor.tests  # prints correlation tests for 1x2, 1x3, 2x3


# FIGURE S3 ----
# MeanPi/Dxy corr plot
# Make a scatterplot of mean_pi / Dxy in two taxa vs. mean_pi / Dxy in two other taxa.
group1A <- "vir"
group1B <- "plumb"
color1 <- "purple"
group2A <- "inor"
group2B <- "hume"
color2 <- "black"
cor.method <- "spearman"
test <- plotMeanPiOverDxyCorr(group1A, group1B, color1,
                              group2A, group2B, color2,
                              autosome.genome.rolling.stats, cor.method)
test   # print test of correlation

# FIGURE S4 left ----
# Dxy_MeanPi plot
# Scatterplot of Dxy vs. mean_pi using compiled genome-wide info;
# Color the points based on Fst.
group1 <- "inor"
group2 <- "hume"
max.x <- 0.010
max.y <- 0.004
cor.method <- "spearman"
test <- plotDxy_MeanPi(group1, group2, cor.method,
                       autosome.genome.rolling.stats,
                       max.x, max.y, color_bar=F)
test   # print test of correlation

# FIGURE S4 right ----
# Dxy_MeanPi plot
# Scatterplot of Dxy vs. mean_pi using compiled genome-wide info;
# Color the points based on Fst.
group1 <- "vir"
group2 <- "plumb"
max.x <- 0.010
max.y <- 0.004
cor.method <- "spearman"
test <- plotDxy_MeanPi(group1, group2, cor.method,
                       autosome.genome.rolling.stats,
                       max.x, max.y, color_bar=F)
test   # print test of correlation


# FIGURE S5 ----
# script to compare autosomal and Z-chromosome differentiation statistics
# compile autosomal-genome-wide info on WC84_Fst, Dxy, pi: done above under "Outgroup Analysis" heading
# compile Z-chromosome info on WC84_Fst, Dxy, pi:
chromosomes.to.combine <- c("Z")
Z.rolling.stats <- compileWindowedStats(base.name, tag.name, chromosomes.to.combine, window.size)

# Left plot in figure S5:
group1 <- "fusc"  
group2 <- "burk"
limits <- c(0, 0.06)
groups.for.graph <- paste0(group1, "_", group2)
Dxy.row.choice.autosome <- which(rownames(autosome.genome.rolling.stats$Dxy) == groups.for.graph)
autosome.Dxy <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice.autosome,])
Dxy.row.choice.Z <- which(rownames(Z.rolling.stats$Dxy) == groups.for.graph)
Z.Dxy <- as.vector(Z.rolling.stats$Dxy[Dxy.row.choice.Z,])
hist.autosome.Dxy <- hist(autosome.Dxy, breaks=seq(limits[1], limits[2]+0.005, by=0.005), plot=FALSE)
hist.Z.Dxy <- hist(Z.Dxy, breaks=seq(limits[1], limits[2]+0.005, by=0.005), plot=FALSE)
quartz(title=paste0("Autosomes vs. Z: Proportions of Dxy between ", group1,"_", group2, sep=""), width=4, height=3)
plot(NA, xlim=limits, ylim=c(0,0.4), xaxp=c(0, 0.06, 3), yaxp=c(0, 0.4, 2),  cex.axis=0.8, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
lines(hist.autosome.Dxy$mids, hist.autosome.Dxy$counts/sum(hist.autosome.Dxy$counts), type="o", pch=16, col="grey50", lwd=2)
lines(hist.Z.Dxy$mids, hist.Z.Dxy$counts/sum(hist.Z.Dxy$counts), type="o", pch=16, col="blue", lwd=2)
label.size <- 1.2
title(xlab=expression(italic('D')[xy]), line=2.5, cex.lab=label.size)
title(ylab="Proportion", line=2, cex.lab=label.size)
t.test(autosome.Dxy, Z.Dxy, var.equal=T)

# Right plot in figure S5:
group1 <- "troch"  
group2 <- "burk"
limits <- c(0, 0.03)
groups.for.graph <- paste0(group1, "_", group2)
Dxy.row.choice.autosome <- which(rownames(autosome.genome.rolling.stats$Dxy) == groups.for.graph)
autosome.Dxy <- as.vector(autosome.genome.rolling.stats$Dxy[Dxy.row.choice.autosome,])
Dxy.row.choice.Z <- which(rownames(Z.rolling.stats$Dxy) == groups.for.graph)
Z.Dxy <- as.vector(Z.rolling.stats$Dxy[Dxy.row.choice.Z,])
hist.autosome.Dxy <- hist(autosome.Dxy, breaks=seq(limits[1], limits[2]+0.0025, by=0.0025), plot=FALSE)
hist.Z.Dxy <- hist(Z.Dxy, breaks=seq(limits[1], limits[2]+0.0025, by=0.0025), plot=FALSE)
quartz(title=paste0("Autosomes vs. Z: Proportions of Dxy between ", group1,"_", group2, sep=""), width=4, height=3)
plot(NA, xlim=limits, ylim=c(0,0.4), xaxp=c(0, 0.03, 3), yaxp=c(0, 0.4, 2),  cex.axis=0.8, xlab=NA, ylab=NA, mgp=c(3,0.5,0))
lines(hist.autosome.Dxy$mids, hist.autosome.Dxy$counts/sum(hist.autosome.Dxy$counts), type="o", pch=16, col="grey50", lwd=2)
lines(hist.Z.Dxy$mids, hist.Z.Dxy$counts/sum(hist.Z.Dxy$counts), type="o", pch=16, col="blue", lwd=2)
label.size <- 1.2
title(xlab=expression(italic('D')[xy]), line=2.5, cex.lab=label.size)
title(ylab="Proportion", line=2, cex.lab=label.size)
t.test(autosome.Dxy, Z.Dxy, var.equal=T)


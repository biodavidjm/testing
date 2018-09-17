

library(data.table)
library(bit64)
library(reshape)
library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(gplots)
library(ggalt)
library(ggrepel)
library(UpSetR)

library(d3heatmap)
hmcol <- colorRampPalette(brewer.pal(10, "RdBu"))(256)
cols6 <- brewer.pal(9, "YlGnBu")[1:6]
cbbPalette <- c("#0072B2", "#D55E00", "#009E73", "#000000", "#E69F00", "#56B4E9", "#F0E442", "#CC79A7")


qcut <- function(x, n) {
  quantiles = seq(0, 1, length.out = n+1)
  cutpoints = unname(quantile(x, quantiles, na.rm = TRUE))
  as.character(cut(x, cutpoints, labels = 1:n, include.lowest = TRUE))
}


sampleOverlap <- function(data,sampleID='Experiment',referenceID='Sequence'){
  pept <- by(data,data[,sampleID],function(x) unique(as.character(x[,referenceID])))
  pepmatch <- lapply(pept, function(x,y) { lapply(y, function(y) length(intersect(x,y))/length(union(x,y))) }, pept)
  m <- matrix(unlist(pepmatch),nrow=length(pept),ncol=length(pept),byrow=TRUE)
  colnames(m) <- rownames(m) <- names(pept)
  ans <- list(Peptides=pept,M=m)
  return(ans)
}

panel.cor <- function(x, y, digits=2, prefix="", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="pairwise.complete.obs", method="pearson"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.cor.log <- function(x, y, digits=2, prefix="", cex.cor, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(log2(x), log2(y), use="pairwise.complete.obs", method="pearson"))
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  txt <- paste(prefix, txt, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

panel.hist <- function(x, ...){
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col="cyan", ...)
}

panel.ma <- function (x, y,  pch = par("pch"), col.lm = "red",  ...) {   
  ymin <- min(y)
  ymax <- max(y)
  xmin <- min(x)
  xmax <- max(x)
  ylim <- c(min(ymin,xmin),max(ymax,xmax))
  xlim <- ylim
  plot(mean(x,y), log2(x)-log2(y), pch = pch,...)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) 
    abline(0,1, col = col.lm, ...)
}


setwd("~/Box Sync/apps/artMS/1D")
evidence_file <- "evidence.txt"
keys_file <- "keys.txt"
summary_file <- "summary.txt"



mergeEvidenceKeys <- function(evidence_file, keys_file) {
  
  #suppressMessages(library(reshape2))
  #suppressMessages(library(data.table)) 
  # data.table here is just to use the "setnames" function 
  #(lazy motherfucker mode)
  
  data <- read.delim(evidence_file, sep='\t', quote = "", header = T, stringsAsFactors = F)
  keys <- read.delim(keys_file, sep='\t', quote = "", header = T, stringsAsFactors = F)
  
  if( !('RawFile' %in% colnames(data)) ) {
    tryCatch(setnames(data, 'Raw.file', 'RawFile'), error=function(e) stop('\nRaw.file not found in EVIDENCE FILE\n'))
  }
  
  if('Raw.file' %in% colnames(keys)){
    setnames(keys, 'Raw.file', 'RawFile')
  }
  
  # Check that the keys file is correct
  if(any(!c('RawFile','IsotopeLabelType','Condition','BioReplicate','Run') %in% colnames(keys))){ #,'SAINT','BioReplicaSaint'
    cat('\nERROR: COLUMN NAMES IN KEYS NOT CONFORM TO SCHEMA. One of these is lost\n
        \tRawFile\n\tIsotopeLabelType\n\tCondition\n\tBioReplicate\n\tRun\n\n') # \tSAINT\n\tBioReplicaSaint\n\n
    stop('Please, try again once revised\n\n')
  }
  
  # MERGING THE DATA
  # Checking that the keys make sense
  unique_data <- unique(data$RawFile)
  unique_keys <- unique(keys$RawFile)
  
  # Rawfiles on Keys not found on the data
  keys_not_found <- setdiff(unique_keys, unique_data)
  # Rawfiles on Data not found on the keys
  data_not_found <- setdiff(unique_data, unique_keys)
  
  if ( (length(keys_not_found) != 0) & ( length(data_not_found) != 0) ) {
    cat(sprintf("keys found: %s \t keys not in data file:\n%s\n", length(unique_keys)-length(keys_not_found), paste(keys_not_found,collapse='\t')))
    cat(sprintf("data found: %s \t data not in keys file:\n%s\n", length(unique_data)-length(data_not_found), paste(data_not_found, collapse='\t')))
    stop('\nThis script is sorry, but it needs to stop this because something is going on between your keys and evidence files so you better check\n')
  }
  
  ## select only required attributes from MQ format
  datamerged <- merge(data, keys, by='RawFile')
  
  return(datamerged)
}


mergeSummaryKeys <- function(summary_file, keys_file) {
  
  summary <- read.delim(summary_file, sep='\t', quote = "", header = T, stringsAsFactors = F)
  summary <- subset(summary, Experiment != "")
  
  keys <- read.delim(keys_file, sep='\t', quote = "", header = T, stringsAsFactors = F)
  
  if( !('RawFile' %in% colnames(summary)) ) {
    tryCatch(setnames(summary, 'Raw.file', 'RawFile'), error=function(e) stop('\nRaw.file not found in SUMMARY FILE\n'))
  }
  
  if('Raw.file' %in% colnames(keys)){
    setnames(keys, 'Raw.file', 'RawFile')
  }
  
  # Check that the keys file is correct
  if(any(!c('RawFile','IsotopeLabelType','Condition','BioReplicate','Run') %in% colnames(keys))){ #,'SAINT','BioReplicaSaint'
    cat('\nERROR: COLUMN NAMES IN KEYS NOT CONFORM TO SCHEMA. One of these is lost\n
        \tRawFile\n\tIsotopeLabelType\n\tCondition\n\tBioReplicate\n\tRun\n\n') # \tSAINT\n\tBioReplicaSaint\n\n
    stop('Please, try again once revised\n\n')
  }
  
  # MERGING THE DATA
  # Checking that the keys make sense
  unique_summary <- unique(summary$RawFile)
  unique_keys <- unique(keys$RawFile)
  
  # Rawfiles on Keys not found on the data
  keys_not_found <- setdiff(unique_keys, unique_summary)
  # Rawfiles on Data not found on the keys
  summary_not_found <- setdiff(unique_summary, unique_keys)
  
  if ( (length(keys_not_found) != 0) & ( length(summary_not_found) != 0) ) {
    cat(sprintf("keys found: %s \t keys not in data file:\n%s\n", length(unique_keys)-length(keys_not_found), paste(keys_not_found,collapse='\t')))
    cat(sprintf("data found: %s \t data not in keys file:\n%s\n", length(unique_summary)-length(summary_not_found), paste(summary_not_found, collapse='\t')))
    stop('\nThis script is sorry, but it needs to stop this because something is going on between your keys and evidence files so you better check\n')
  }
  
  ## select only required attributes from MQ format
  datamerged <- merge(summary, keys, by='RawFile')
  
  return(datamerged)
}




evidencekeys <- mergeEvidenceKeys(evidence_file, keys_file)
summarykeys <- mergeSummaryKeys(summary_file, keys_file)



evidence1 <- evidencekeys
colnames(evidence1) <- gsub("Experiment.x", "Experiment", colnames(evidence1))
colnames(evidence1) <- tolower(colnames(evidence1))
evidence1 <- subset(evidence1, reverse != "+")
# Add potential.contaminant column in case search was performed without it
`%ni%` <- Negate(`%in%`)
if("potential.contaminant" %ni% colnames(evidence1)){
  evidence1$potential.contaminant <- ""
}


evidence1.dt <- data.table(evidence1)

evidence2 <- ddply(evidence1, c("experiment", "condition", "potential.contaminant"), summarise, 
                   PSMs = sum(ms.ms.count), Ions = length(sequence[ms.ms.count>0]), Peptides = length(unique(sequence[ms.ms.count>0])), 
                   Proteins = length(unique(proteins[ms.ms.count>0])), Intensity = as.numeric(format(sum(intensity, na.rm=TRUE), scientific=TRUE, digits = 2)), 
                   PeakWidth.mean = mean(retention.length, na.rm=TRUE))

if("fraction" %in% colnames(evidence1)){
  
  evidence2fx <- ddply(evidence1, c("experiment", "condition", "potential.contaminant", "fraction"), summarise, 
                       PSMs = sum(ms.ms.count), Ions = length(sequence[ms.ms.count>0]), Peptides = length(unique(sequence[ms.ms.count>0])), 
                       Proteins = length(unique(proteins[ms.ms.count>0])), Intensity = as.numeric(format(sum(intensity, na.rm=TRUE), scientific=TRUE, digits = 2)), 
                       PeakWidth.mean = mean(retention.length, na.rm=TRUE))
  
}


evidence3 <- ddply(evidence2, c("condition", "potential.contaminant"), summarise, 
                   PSMs.mean = mean(PSMs),PSMs.max = max(PSMs),PSMs.min = min(PSMs),PSMs.sem = sd(PSMs)/sqrt(length(PSMs)),
                   Ions.mean = mean(Ions), Ions.max = max(Ions), Ions.min = min(Ions), Ions.sem = sd(Ions)/sqrt(length(Ions)),
                   Peptides.mean = mean(Peptides), Peptides.max = max(Peptides),Peptides.min = min(Peptides), Peptides.sem = sd(Peptides)/sqrt(length(Peptides)),
                   Proteins.mean = mean(Proteins), Proteins.max = max(Proteins),Proteins.min = min(Proteins),Proteins.sem = sd(Proteins)/sqrt(length(Proteins)),
                   Intensity.mean = mean(Intensity),Intensity.max = max(Intensity),Intensity.min = min(Intensity),Intensity.sem = sd(Intensity)/sqrt(length(Intensity)),
                   PeakWidth.mean.mean = mean(PeakWidth.mean),PeakWidth.mean.max = max(PeakWidth.mean),PeakWidth.mean.min = min(PeakWidth.mean),PeakWidth.mean.sem = sd(PeakWidth.mean)/sqrt(length(PeakWidth.mean)))




summary1 <- summarykeys
colnames(summary1) <- gsub("Experiment.x", "Experiment", colnames(summary1))
colnames(summary1) <- tolower(colnames(summary1))

if("fraction" %in% colnames(summary1)){
  
  summary2fx <- summary1
}

if("fraction" %in% colnames(summary1)){
  
  summary1 <- data.table(summary1[,c("condition", "experiment", "ms", "ms.ms", "ms.ms.identified....", "isotope.patterns", "isotope.patterns.sequenced..z.1.")])
  #summary1 <- aggregate.data.frame(summary1[,-1:-2], by=list(condition=summary1$condition, experiment=summary1$experiment), FUN="sum")
  summary1 <- summary1[,list(ms=sum(ms), ms.ms=sum(ms.ms), ms.ms.identified....=mean(ms.ms.identified....), isotope.patterns=sum(isotope.patterns), isotope.patterns.sequenced..z.1.=sum(isotope.patterns.sequenced..z.1.)),
                       by=list(condition, experiment)]
  
} else {
  summary1 <- data.table(summary1[,c("condition", "experiment", "ms", "ms.ms", "ms.ms.identified....", "isotope.patterns", "isotope.patterns.sequenced..z.1.")])
}

summary2 <- ddply(summary1, c("condition"), summarise, 
                  num.MS1.mean = mean(ms), num.MS1.max = max(ms), num.MS1.min = min(ms), num.MS1.sem = sd(ms)/sqrt(length(ms)),
                  num.MS2.mean = mean(ms.ms), num.MS2.max = max(ms.ms), num.MS2.min = min(ms.ms), num.MS2.sem = sd(ms.ms)/sqrt(length(ms.ms)),
                  num.IsotopePatterns.mean = mean(isotope.patterns), num.IsotopePatterns.max = max(isotope.patterns),
                  num.IsotopePatterns.min = min(isotope.patterns), num.IsotopePatterns.sem = sd(isotope.patterns)/sqrt(length(isotope.patterns)),
                  num.IsotopePatternsSeq.mean = mean(isotope.patterns.sequenced..z.1.), num.IsotopePatternsSeq.max = max(isotope.patterns.sequenced..z.1.),
                  num.IsotopePatternsSeq.min = min(isotope.patterns.sequenced..z.1.), num.IsotopePatternsSeq.sem = sd(isotope.patterns.sequenced..z.1.)/sqrt(length(isotope.patterns.sequenced..z.1.)),
                  pct.MS2Id.mean = mean(ms.ms.identified....), pct.MS2Id.max = max(ms.ms.identified....),
                  pct.MS2Id.min = min(ms.ms.identified....), pct.MS2Id.sem = sd(ms.ms.identified....)/sqrt(length(ms.ms.identified....))
                  ) 








# Oversampling
## Oversampling is calculated for peptides identified by MS/MS and detected in MS1 (i.e., peptides with intensity calculated)
oversampling0 <- evidence1.dt[, .N, by=list(ms.ms.count, experiment)]
oversampling0$MSMS.counts <- findInterval(oversampling0$ms.ms.count, seq(1, 3, by=1)) 
oversampling0 <- oversampling0[, list(N=sum(N)), by=list(MSMS.counts, experiment)]
oversampling0$MSMS.counts <- paste("n=", oversampling0$MSMS.counts, sep="")
oversampling0$MSMS.counts <- gsub("n=3", "n=3+", oversampling0$MSMS.counts)
oversampling0.total <- evidence1.dt[, .N, by=list(experiment)]
oversampling0 <- merge(oversampling0, oversampling0.total, by="experiment", all = T)
oversampling0$FxOverSamp <- as.numeric(format(100*(oversampling0$N.x/oversampling0$N.y), digits = 3))

oversampling1 <- evidence1.dt[intensity > 0, .N, by=list(ms.ms.count, experiment)]
oversampling1$MSMS.counts <- findInterval(oversampling1$ms.ms.count, seq(1, 3, by=1)) 
oversampling1 <- oversampling1[, list(N=sum(N)), by=list(MSMS.counts, experiment)]
oversampling1$MSMS.counts <- paste("n=", oversampling1$MSMS.counts, sep="")
oversampling1$MSMS.counts <- gsub("n=3", "n=3+", oversampling1$MSMS.counts)
oversampling1.total <- evidence1.dt[intensity > 0, .N, by=list(experiment)]
oversampling1 <- merge(oversampling1, oversampling1.total, by="experiment", all = T)
oversampling1$FxOverSamp <- as.numeric(format(100*(oversampling1$N.x/oversampling1$N.y), digits = 3))

oversampling2 <- evidence1.dt[intensity > 0 & ms.ms.count > 0, .N, by=list(ms.ms.count, experiment)]
oversampling2$MSMS.counts <- findInterval(oversampling2$ms.ms.count, seq(1, 3, by=1)) 
oversampling2 <- oversampling2[, list(N=sum(N)), by=list(MSMS.counts, experiment)]
oversampling2$MSMS.counts <- paste("n=", oversampling2$MSMS.counts, sep="")
oversampling2$MSMS.counts <- gsub("n=3", "n=3+", oversampling2$MSMS.counts)
oversampling2.total <- evidence1.dt[intensity > 0 & ms.ms.count > 0, .N, by=list(experiment)]
oversampling2 <- merge(oversampling2, oversampling2.total, by="experiment", all = T)
oversampling2$FxOverSamp <- as.numeric(format(100*(oversampling2$N.x/oversampling2$N.y), digits = 3))



# Charge state
chargeState <- evidence1.dt[, .N, by=list(charge, experiment)]
chargeState$charge <- paste("z=", chargeState$charge, sep="")
chargeState.total <- evidence1.dt[, .N, by=list(experiment)]
chargeState <- merge(chargeState, chargeState.total, by="experiment", all = T)
chargeState$FxOverSamp <- as.numeric(format(100*(chargeState$N.x/chargeState$N.y), digits = 1))

chargeStateCond <- evidence1.dt[, .N, by=list(charge, condition)]
chargeStateCond$charge <- paste("z=", chargeStateCond$charge, sep="")
chargeStateCond.total <- evidence1.dt[, .N, by=list(condition)]
chargeStateCond <- merge(chargeStateCond, chargeStateCond.total, by="condition", all = T)
chargeStateCond$FxOverSamp <- as.numeric(format(100*(chargeStateCond$N.x/chargeStateCond$N.y), digits = 1))

# Type
mstype <- evidence1.dt[, .N, by=list(type, experiment, condition)]
mstype.total <- evidence1.dt[, .N, by=list(experiment)]
mstype <- merge(mstype, mstype.total, by="experiment", all = T)
mstype$FxOverSamp <- as.numeric(format(100*(mstype$N.x/mstype$N.y), digits = 2))




############# PLOTS ###############


nsamples <- length(unique(evidence1$experiment))
nconditions <- length(unique(evidence1$condition))

# Check the largest number so that width and height of pdf is not too large
if(nsamples > 20){
  nsamples <- 20
}

if(nconditions > 7){
  nconditions <- 7
}

pdf('QC-Plots%03d.pdf', width=nsamples*3, height=20, onefile = FALSE)


### PSM ###

ggplot(evidence2, aes(x = experiment, y = PSMs, fill = factor(condition))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_text(aes(label=round(PSMs, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
  facet_wrap(~potential.contaminant, ncol=1) +
  xlab("Experiment") + ylab("Counts") +
  ggtitle("Number of PSMs: bottom = Potential contaminants; top = non-contaminants") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")

ggplot(evidence3, aes(x = condition, y = PSMs.mean, fill = factor(potential.contaminant))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_errorbar(aes(ymin=PSMs.mean-PSMs.sem, ymax=PSMs.mean+PSMs.sem), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=round(PSMs.mean, digits=0)), hjust=0.5, vjust=-1.5, size = 10, position = position_dodge(width = 1)) +
  xlab("Condition") + ylab("Counts") +
  ggtitle("Mean number of PSMs per condition for contaminants and non-contaminants, error bar= std error of the mean") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")

if("fraction" %in% colnames(evidence1)){
  
  ggplot(evidence2fx, aes(x = experiment, y = PSMs, fill = factor(fraction))) +
    geom_bar(stat="identity", alpha=0.7) +
    geom_text(aes(label=round(PSMs, digits=0)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
    facet_wrap(~potential.contaminant, ncol=1, scales = "free") +
    xlab("Experiment") + ylab("Counts") +
    ggtitle("Number of PSMs per Fraction: bottom = Potential contaminants; top = non-contaminants") +
    theme(legend.text = element_text(size=20)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
    theme(axis.text.y = element_text(size=20)) + 
    theme(axis.title.x = element_text(size=30)) + 
    theme(axis.title.y = element_text(size=30)) + 
    theme(plot.title = element_text(size = 40)) 
  
}


### IONS ###

ggplot(evidence2, aes(x = experiment, y = Ions, fill = factor(condition))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_text(aes(label=round(Ions, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
  facet_wrap(~potential.contaminant, ncol=1) +
  xlab("Experiment") + ylab("Counts") +
  ggtitle("Number of unique Peptide Ions: bottom = Potential contaminants; top = non-contaminants") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")

ggplot(evidence3, aes(x = condition, y = Ions.mean, fill = factor(potential.contaminant))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_errorbar(aes(ymin=Ions.mean-Ions.sem, ymax=Ions.mean+Ions.sem), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=round(Ions.mean, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
  xlab("Condition") + ylab("Counts") +
  ggtitle("Mean number of unique Peptide Ions per condition for contaminants (blue) and non-contaminants (red), error bar= std error of the mean") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Set1")

if("fraction" %in% colnames(evidence1)){
  
  ggplot(evidence2fx, aes(x = experiment, y = Ions, fill = factor(fraction))) +
    geom_bar(stat="identity", alpha=0.7) +
    geom_text(aes(label=round(Ions, digits=0)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
    facet_wrap(~potential.contaminant, ncol=1, scales = "free") +
    xlab("Experiment") + ylab("Counts") +
    ggtitle("Number of unique Peptide Ions in each Fraction: bottom = Potential contaminants; top = non-contaminants") +
    theme(legend.text = element_text(size=20)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
    theme(axis.text.y = element_text(size=20)) + 
    theme(axis.title.x = element_text(size=30)) + 
    theme(axis.title.y = element_text(size=30)) + 
    theme(plot.title = element_text(size = 40)) 
}



#### TYPE ######
ggplot(mstype, aes(x = experiment, y = FxOverSamp, fill = factor(type))) +
  geom_bar(stat="identity", position = position_stack(), alpha=0.7) +
  geom_text_repel(aes(label=round(FxOverSamp, digits=1)), vjust=1, size = 10, position = position_stack()) +
  xlab("Experiment") + ylab("Fraction") +
  ggtitle("Type of identification (MaxQuant type column)") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Set1")



### PEPTIDES ###

ggplot(evidence2, aes(x = experiment, y = Peptides, fill = factor(condition))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_text(aes(label=round(Peptides, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
  facet_wrap(~potential.contaminant, ncol=1) +
  xlab("Experiment") + ylab("Counts") +
  ggtitle("Number of unique Peptides: bottom = Potential contaminants; top = non-contaminants") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")

ggplot(evidence3, aes(x = condition, y = Peptides.mean, fill = factor(potential.contaminant))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_errorbar(aes(ymin=Peptides.mean-Peptides.sem, ymax=Peptides.mean+Peptides.sem), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=round(Peptides.mean, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
  xlab("Condition") + ylab("Counts") +
  ggtitle("Mean number of unique Peptides per condition for contaminants (blue) and non-contaminants (red), error bar= std error of the mean") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Set1")

if("fraction" %in% colnames(evidence1)){
  
  ggplot(evidence2fx, aes(x = experiment, y = Peptides, fill = factor(fraction))) +
    geom_bar(stat="identity", alpha=0.7) +
    geom_text(aes(label=round(Peptides, digits=0)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
    facet_wrap(~potential.contaminant, ncol=1, scales = "free") +
    xlab("Experiment") + ylab("Counts") +
    ggtitle("Number of unique Peptides in each Fraction: bottom = Potential contaminants; top = non-contaminants") +
    theme(legend.text = element_text(size=20)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
    theme(axis.text.y = element_text(size=20)) + 
    theme(axis.title.x = element_text(size=30)) + 
    theme(axis.title.y = element_text(size=30)) + 
    theme(plot.title = element_text(size = 40)) 
}

lst.pepExp <- list()
for(i in unique(evidence1$experiment)) {
  lst.pepExp[[i]] <- unique(subset(evidence1, ms.ms.count > 0 & experiment == i)[,c("sequence")])
}

upset(fromList(lst.pepExp), nsets=length(unique(evidence1$experiment)), nintersects = 50, order.by = "degree")

upset(fromList(lst.pepExp), nsets=length(unique(evidence1$experiment)), nintersects = 50, order.by = "freq")


lst.pepCond <- list()
for(i in unique(evidence1$condition)) {
  lst.pepCond[[i]] <- unique(subset(evidence1, ms.ms.count > 0 & condition == i)[,c("sequence")])
}

#upset(fromList(lst.pepCond), nsets=length(unique(evidence1$condition)), nintersects = 50, order.by = "degree")

upset(fromList(lst.pepCond), nsets=length(unique(evidence1$condition)), nintersects = 50, order.by = "freq")



# Proteins

ggplot(evidence2, aes(x = experiment, y = Proteins, fill = factor(condition))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_text(aes(label=round(Proteins, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
  facet_wrap(~potential.contaminant, ncol=1) +
  xlab("Experiment") + ylab("Counts") +
  ggtitle("Number of unique Protein Groups: bottom = Potential contaminants; top = non-contaminants") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")

ggplot(evidence3, aes(x = condition, y = Proteins.mean, fill = factor(potential.contaminant))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_errorbar(aes(ymin=Proteins.mean-Proteins.sem, ymax=Proteins.mean+Proteins.sem), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=round(Proteins.mean, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
  xlab("Condition") + ylab("Counts") +
  ggtitle("Mean number of unique Proteins per condition for contaminants (blue) and non-contaminants (red), error bar= std error of the mean") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Set1")



lst.protExp <- list()
for(i in unique(evidence1$experiment)) {
  lst.protExp[[i]] <- unique(subset(evidence1, ms.ms.count > 0 & experiment == i)[,c("proteins")])
}

upset(fromList(lst.protExp), nsets=length(unique(evidence1$experiment)), nintersects = 50, order.by = "degree")

upset(fromList(lst.protExp), nsets=length(unique(evidence1$experiment)), nintersects = 50, order.by = "freq")


lst.protCond <- list()
for(i in unique(evidence1$condition)) {
  lst.protCond[[i]] <- unique(subset(evidence1, ms.ms.count > 0 & condition == i)[,c("proteins")])
}

#upset(fromList(lst.protCond), nsets=length(unique(evidence1$condition)), nintersects = 50, order.by = "degree")

upset(fromList(lst.protCond), nsets=length(unique(evidence1$condition)), nintersects = 50, order.by = "freq")




# Number of MS1 scans
ggplot(summary1, aes(x = experiment, y = ms, fill = condition)) +
  geom_bar(stat="identity", alpha=0.7) +
  geom_text(aes(label=round(ms, digits=0)), vjust=1 , size = 15) +
  xlab("Experiment") + ylab("Counts") +
  ggtitle("Number of MS1 scans") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")



ggplot(summary2, aes(x = condition, y = num.MS1.mean, fill = factor(condition))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_errorbar(aes(ymin=num.MS1.mean-num.MS1.sem, ymax=num.MS1.mean+num.MS1.sem), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=round(num.MS1.mean, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
  xlab("Condition") + ylab("Counts") +
  ggtitle("Mean number of MS1 scans per condition, error bar= std error of the mean") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")



if("fraction" %in% colnames(evidence1)){
  
  ggplot(summary2fx, aes(x = experiment, y = ms, fill = factor(fraction))) +
    geom_bar(stat="identity", alpha=0.7) +
    geom_text(aes(label=round(ms, digits=0)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
    xlab("Experiment") + ylab("Counts") +
    ggtitle("Number of MS1 scans per Fraction") +
    theme(legend.text = element_text(size=20)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
    theme(axis.text.y = element_text(size=20)) + 
    theme(axis.title.x = element_text(size=30)) + 
    theme(axis.title.y = element_text(size=30)) + 
    theme(plot.title = element_text(size = 40)) 
  
}


#Number of MS2 scans
ggplot(summary1, aes(x = experiment, y = ms.ms, fill = condition)) +
  geom_bar(stat="identity", alpha=0.7) +
  geom_text(aes(label=round(ms.ms, digits=0)), vjust=1 , size = 15) +
  xlab("Experiment") + ylab("Counts") +
  ggtitle("Number of MS2 scans") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")


ggplot(summary2, aes(x = condition, y = num.MS2.mean, fill = factor(condition))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_errorbar(aes(ymin=num.MS2.mean-num.MS2.sem, ymax=num.MS2.mean+num.MS2.sem), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=round(num.MS2.mean, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
  xlab("Condition") + ylab("Counts") +
  ggtitle("Mean number of MS2 scans per condition, error bar= std error of the mean") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")

if("fraction" %in% colnames(evidence1)){
  
  ggplot(summary2fx, aes(x = experiment, y = ms.ms, fill = factor(fraction))) +
    geom_bar(stat="identity", alpha=0.7) +
    geom_text(aes(label=round(ms.ms, digits=0)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
    xlab("Experiment") + ylab("Counts") +
    ggtitle("Number of MS2 per Fraction") +
    theme(legend.text = element_text(size=20)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
    theme(axis.text.y = element_text(size=20)) + 
    theme(axis.title.x = element_text(size=30)) + 
    theme(axis.title.y = element_text(size=30)) + 
    theme(plot.title = element_text(size = 40)) 
  
}



summary1.scans <- melt(summary1[,1:4], id.vars=1:2)
ggplot(summary1.scans, aes(x = interaction(variable, experiment), y = value, fill = condition)) +
  geom_bar(stat="identity", alpha=0.7) +
  geom_text(aes(label=round(value, digits=0)), angle = 90, vjust=0.5 , size = 10) +
  xlab("Experiment") + ylab("Counts") +
  ggtitle("Number of MS1 and MS2 scans") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")



# Number of msms.identification rate
ggplot(summary1, aes(x = experiment, y = ms.ms.identified...., fill = condition)) +
  geom_bar(stat="identity", alpha=0.7) +
  geom_text(aes(label=round(ms.ms.identified...., digits=2)), vjust=1 , size = 15) +
  xlab("Experiment") + ylab("Rate") +
  ggtitle("MS2 Identification rate") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")


ggplot(summary2, aes(x = condition, y = pct.MS2Id.mean, fill = factor(condition))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_errorbar(aes(ymin=pct.MS2Id.mean-pct.MS2Id.sem, ymax=pct.MS2Id.mean+pct.MS2Id.sem), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=round(pct.MS2Id.mean, digits=2)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
  xlab("Condition") + ylab("Rate") +
  ggtitle("Mean MS2 Identification rate across experiments and fractions") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")

if("fraction" %in% colnames(evidence1)){
  
  ggplot(summary2fx, aes(x = experiment, y = ms.ms.identified...., fill = factor(fraction))) +
    geom_bar(stat="identity", alpha=0.7) +
    geom_text(aes(label=round(ms.ms.identified...., digits=1)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
    xlab("Experiment") + ylab("Counts") +
    ggtitle("MS2 Identification rate per Fraction") +
    theme(legend.text = element_text(size=20)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
    theme(axis.text.y = element_text(size=20)) + 
    theme(axis.title.x = element_text(size=30)) + 
    theme(axis.title.y = element_text(size=30)) + 
    theme(plot.title = element_text(size = 40)) 
  
}



# Number of isotope patterns
ggplot(summary1, aes(x = experiment, y = isotope.patterns, fill = condition)) +
  geom_bar(stat="identity", alpha=0.7) +
  geom_text(aes(label=round(isotope.patterns, digits=0)), vjust=1 , size = 15) +
  xlab("Experiment") + ylab("Counts") +
  ggtitle("Number of detected Isotope Patterns") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")


ggplot(summary2, aes(x = condition, y = num.IsotopePatterns.mean, fill = factor(condition))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_errorbar(aes(ymin=num.IsotopePatterns.mean-num.IsotopePatterns.sem, ymax=num.IsotopePatterns.mean+num.IsotopePatterns.sem), width=.2, position=position_dodge(.9)) +
  geom_text(aes(label=round(num.IsotopePatterns.mean, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
  xlab("Condition") + ylab("Counts") +
  ggtitle("Mean number of detected Isotope Patterns") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")

if("fraction" %in% colnames(evidence1)){
  
  ggplot(summary2fx, aes(x = experiment, y = isotope.patterns, fill = factor(fraction))) +
    geom_bar(stat="identity", alpha=0.7) +
    geom_text(aes(label=round(isotope.patterns, digits=0)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
    xlab("Experiment") + ylab("Counts") +
    ggtitle("Number of detected Isotope Patterns per Fraction") +
    theme(legend.text = element_text(size=20)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
    theme(axis.text.y = element_text(size=20)) + 
    theme(axis.title.x = element_text(size=30)) + 
    theme(axis.title.y = element_text(size=30)) + 
    theme(plot.title = element_text(size = 40)) 
  
}




# Number of sequenced isotope patterns with charge = 2 or more
ggplot(summary1, aes(x = experiment, y = isotope.patterns.sequenced..z.1., fill = condition)) +
  geom_bar(stat="identity", alpha=0.7) +
  geom_text(aes(label=round(isotope.patterns.sequenced..z.1., digits=0)), vjust=1 , size = 15) +
  xlab("Experiment") + ylab("Counts") +
  ggtitle("Number of sequenced Isotope Patterns with charge state greater than 1") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")


ggplot(summary2, aes(x = condition, y = num.IsotopePatternsSeq.mean, fill = factor(condition))) +
  geom_bar(stat="identity", position = position_dodge(width = 1), alpha=0.7) +
  geom_errorbar(aes(ymin=num.IsotopePatternsSeq.mean-num.IsotopePatternsSeq.sem, ymax=num.IsotopePatternsSeq.mean+num.IsotopePatternsSeq.sem), width=.2, position=position_dodge(.9)) +
  geom_text_repel(aes(label=round(num.IsotopePatternsSeq.mean, digits=0)), hjust=0.5, vjust=-0.5, size = 10, position = position_dodge(width = 1)) +
  xlab("Condition") + ylab("Counts") +
  ggtitle("Mean number of sequenced Isotope Patterns with charge state greater than 1") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")

if("fraction" %in% colnames(evidence1)){
  
  ggplot(summary2fx, aes(x = experiment, y = isotope.patterns.sequenced..z.1., fill = factor(fraction))) +
    geom_bar(stat="identity", alpha=0.7) +
    geom_text(aes(label=round(isotope.patterns.sequenced..z.1., digits=0)), hjust=0.5, vjust=1.5, size = 7, position = position_stack()) +
    xlab("Experiment") + ylab("Counts") +
    ggtitle("Number of sequenced Isotope Patterns  with charge state greater than 1 per Fraction") +
    theme(legend.text = element_text(size=20)) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, size=20)) + 
    theme(axis.text.y = element_text(size=20)) + 
    theme(axis.title.x = element_text(size=30)) + 
    theme(axis.title.y = element_text(size=30)) + 
    theme(plot.title = element_text(size = 40)) 
  
}


# Peptide ion oversampling
ggplot(oversampling0[with(oversampling0, order(MSMS.counts)),], aes(x = experiment, y = FxOverSamp, fill = MSMS.counts, label=FxOverSamp)) + 
  geom_bar(stat="identity", alpha=0.7) + 
  geom_text(aes(label=round(FxOverSamp, digits=1)), vjust=1 , size = 10, position = position_stack()) +
  xlab("Experiment") + ylab("Fraction (percentage)") +
  ggtitle("Peptide ion oversampling, all peptides reported by MaxQuant") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Set1")

ggplot(oversampling1[with(oversampling1, order(MSMS.counts)),], aes(x = experiment, y = FxOverSamp, fill = MSMS.counts, label=FxOverSamp)) + 
  geom_bar(stat="identity", alpha=0.7) + 
  geom_text(aes(label=round(FxOverSamp, digits=1)), vjust=1 , size = 10, position = position_stack()) +
  xlab("Experiment") + ylab("Fraction (percentage)") +
  ggtitle("Peptide ion oversampling, only peptides detected (MS1 AUC calculated)") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Set1")


ggplot(oversampling2[with(oversampling2, order(MSMS.counts)),], aes(x = experiment, y = FxOverSamp, fill = MSMS.counts, label=FxOverSamp)) + 
  geom_bar(stat="identity", alpha=0.7) + 
  geom_text(aes(label=round(FxOverSamp, digits=1)), vjust=1 , size = 10, position = position_stack()) +
  xlab("Experiment") + ylab("Fraction (percentage)") +
  ggtitle("Peptide ion oversampling, only peptides detected (MS1 AUC calculated) and identified (confidence MS/MS)") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Set1")

# Charge state distribution
ggplot(chargeState[with(chargeState, order(charge)),], aes(x = experiment, y = FxOverSamp, fill = charge)) + 
  geom_bar(stat="identity", alpha=0.7) + 
  geom_text_repel(aes(label=round(FxOverSamp, digits=1)), vjust=-0.5, size = 5, position = position_stack()) +
  xlab("Experiment") + ylab("Fraction (percentage)") +
  ggtitle("Precursor charge state distribution") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")


# Mass error
ggplot(evidence1, aes(x=experiment, y=uncalibrated.mass.error..ppm.)) + 
  geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha=0.7) +
  geom_text(data = aggregate(uncalibrated.mass.error..ppm. ~ experiment, evidence1, median), aes(label = round(uncalibrated.mass.error..ppm., digits=1), y = max(evidence1$uncalibrated.mass.error..ppm., na.rm=T) + 2), size=20) +
  xlab("Experiment") + ylab("mass error") +
  ggtitle("Precursor mass error (in ppm) distribution, global median mass error on the top") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")

# mass-over-charge distribution
ggplot(evidence1, aes(x=experiment, y=m.z)) +  
  geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha=0.7) +
  geom_text(data = aggregate(m.z ~ experiment, evidence1, median), aes(label = round(m.z, digits=1), y = max(evidence1$m.z, na.rm=T) + 30), size=20) +
  xlab("Experiment") + ylab("m/z") +
  ggtitle("Precursor mass-over-charge distribution, global median m/z on the top") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")




#m Peptide Intensity CV
peptCV <- data.table(subset(evidence1, !is.na(intensity)))
peptCV <- peptCV[, list(intensity=sum(intensity, na.rm=T)), by=list(condition, experiment, modified.sequence)]
peptCV <- peptCV[, list(pCV = 100*(sd(intensity)/mean(intensity)), sumInt=sum(intensity), pDX = length(unique(experiment))), by=list(condition, modified.sequence)]
peptCV <- peptCV[, bin.all := qcut(sumInt, 4)]
peptCV <- peptCV[, bin.condition := qcut(sumInt, 4), by = condition]

ggplot(subset(peptCV, !is.na(pCV)), aes(condition, pCV)) + 
  geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha=0.7) + 
  geom_text(data = aggregate(pCV ~ condition, subset(peptCV, !is.na(pCV)), median), aes(label = round(pCV, digits=1), y = max(peptCV$pCV, na.rm=T) + 1), size = 10) +
  geom_text(data = aggregate(pCV ~ condition, subset(peptCV, !is.na(pCV)), length), aes(label = round(pCV, digits=1), y = 0), size = 10) +
  xlab("Condition") + ylab("Coefficient of variance (%)") +
  ggtitle("Distribution of peptide feature intensity CV within each condition \n 
          Overall median CV for each condition is given on the top and number of features used to calculate CVs is shown on the bottom") +
  
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")


ggplot(subset(peptCV, !is.na(pCV)), aes(interaction(condition, bin.condition), pCV)) + 
  geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha=0.7) + 
  geom_text(data = aggregate(pCV ~ condition + bin.condition, subset(peptCV, !is.na(pCV)), median), aes(label = round(pCV, digits=1), y = max(peptCV$pCV, na.rm=T) + 1), size = 7, angle=90) +
  geom_text(data = aggregate(pCV ~ condition + bin.condition, subset(peptCV, !is.na(pCV)), length), aes(label = round(pCV, digits=1), y = 0), size = 7, angle=90) +
  xlab("Condition") + ylab("Coefficient of variance (%)") +
  ggtitle("Distribution of peptide feature intensity CV within each condition \n 
          For each condition, peptides were ranked by summed intensity and the CV for each peptide was calculated, \n 
          therefore each condition shows 4 distribution (box) for low (1) to high (4) intensity peptides \n 
          Overall median CV within each bin/condition is shown on the top and number of features used to calculate CV is given on the bottom") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 20, hjust = 0.5)) + 
  scale_fill_brewer(palette="Spectral")



# peptide detection (using modified.sequence)
peptDX <- peptCV[, .N, by=list(condition, pDX)]
peptDXT <- peptDX[, list(N.total = sum(N)), by=list(condition)]
peptDX <- merge(peptDX, peptDXT, by="condition")
peptDX$fxDx <- 100*(peptDX$N/peptDX$N.total)
ggplot(peptDX, aes(x=condition, y=fxDx, fill=factor(pDX))) + 
  geom_bar(stat="identity", alpha=0.7) + 
  geom_text_repel(aes(label=round(fxDx, digits=1)), vjust=1, size = 10, position = position_stack()) +
  xlab("Condition") + ylab("Counts") +
  ggtitle("Frequency of peptides detection across replicates within the same condition") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(plot.title = element_text(size = 40)) +
  scale_fill_brewer(palette="Set1")



#m Protein Intensity CV
protCV <- data.table(subset(evidence1, !is.na(intensity)))
protCV <- protCV[, list(intensity=sum(intensity, na.rm=T)), by=list(condition, experiment, proteins)]
protCV <- protCV[, list(pCV = 100*(sd(intensity)/mean(intensity)), sumInt=sum(intensity), pDX = length(unique(experiment))), by=list(condition, proteins)]
protCV <- protCV[, bin.all := qcut(sumInt, 4)]
protCV <- protCV[, bin.condition := qcut(sumInt, 4), by = condition]

ggplot(subset(protCV, !is.na(pCV)), aes(condition, pCV)) + 
  geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha=0.7) + 
  geom_text(data = aggregate(pCV ~ condition, subset(protCV, !is.na(pCV)), median), aes(label = round(pCV, digits=1), y = max(protCV$pCV, na.rm=T) + 1), size = 15) +
  xlab("Condition") + ylab("Coefficient of variance (%)") +
  ggtitle("Distribution of Protein intensity CV within each condition \n 
          Overall median CV for each condition is given on the top and number of proteins used to calculate CVs is shown on the bottom") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 20, hjust = 0.5)) + 
  scale_fill_brewer(palette="Spectral")


ggplot(subset(protCV, !is.na(pCV)), aes(interaction(condition, bin.condition), pCV)) + 
  geom_boxplot(varwidth = TRUE, aes(fill = factor(condition)), alpha=0.7) + 
  geom_text(data = aggregate(pCV ~ condition + bin.condition, subset(protCV, !is.na(pCV)), median), aes(label = round(pCV, digits=1), y = max(protCV$pCV, na.rm=T) + 1), size = 7, angle=90) +
  geom_text(data = aggregate(pCV ~ condition + bin.condition, subset(protCV, !is.na(pCV)), length), aes(label = round(pCV, digits=1), y = 0), size = 7, angle=90) +
  xlab("Condition") + ylab("Coefficient of variance (%)") +
  ggtitle("Distribution of Protein (summed) intensity CV within each condition \n 
          For each condition, Proteins were ranked by summed intensity and the CV for each protein was calculated, \n 
          therefore each condition shows 4 distribution (box) for low (1) to high (4) intensity proteins \n 
          Overall median CV within each bin/condition is shown on the top and number of protein groups used to calculate CV is given on the bottom") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 20, hjust = 0.5)) + 
  scale_fill_brewer(palette="Spectral")



# Protein detection
protDX <- protCV[, .N, by=list(condition, pDX)]
protDXT <- protDX[, list(N.total = sum(N)), by=list(condition)]
protDX <- merge(protDX, protDXT, by="condition")
protDX$fxDx <- 100*(protDX$N/protDX$N.total)
ggplot(protDX, aes(x=condition, y=fxDx, fill=factor(pDX))) + 
  geom_bar(stat="identity", alpha=0.7) + 
  geom_text_repel(aes(label=round(fxDx, digits=1)), vjust=1, size = 10, position = position_stack()) +
  xlab("Condition") + ylab("Counts") +
  ggtitle("Detection of Protein across replicates of the same condition") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 0, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  scale_fill_brewer(palette="Set1")



ggplot(subset(evidence1, !is.na(intensity)), aes(experiment, log2(intensity))) + 
  geom_boxplot(varwidth = TRUE, aes(fill = potential.contaminant), alpha=0.7) + 
  #geom_text_repel(data = aggregate(intensity ~ experiment + potential.contaminant, subset(evidence1, !is.na(intensity)), median), aes(label = round(log2(intensity), digits=1), y = log2(max(evidence1$intensity, na.rm=T))+0.5 ), size = 15) +
  xlab("Experiment") + ylab("Log2 Intensity") +
  ggtitle("Peptide feature intensity distribution") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Set1")

if("fraction" %in% colnames(evidence1)){
  
ggplot(subset(evidence1, !is.na(intensity)), aes(experiment, log2(intensity))) + 
  geom_boxplot(varwidth = TRUE, aes(fill = potential.contaminant), alpha=0.7) + 
  #geom_text_repel(data = aggregate(intensity ~ experiment + potential.contaminant, subset(evidence1, !is.na(intensity)), median), aes(label = round(log2(intensity), digits=1), y = log2(max(evidence1$intensity, na.rm=T))+0.5 ), size = 15) +
  xlab("Experiment") + ylab("Log2 Intensity") +
  facet_wrap(~fraction,ncol=5) +
  ggtitle("Peptide feature intensity distribution by Fraction") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Set1")

}



ggplot(evidence1, aes(x=log2(intensity))) + 
  geom_density(alpha=.5, aes(fill = type)) + 
  ggtitle("Peptide feature intensity distribution by ID Type") +
  facet_wrap(~experiment,ncol=5) +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Set1") 

dev.off()






pdf('ID-Overlap%03d.pdf', width=nsamples*4, height=nsamples*4, onefile = FALSE)
ovlSeq <- subset(evidence1, !is.na(intensity) & !is.na(ms.ms.count), select=c("experiment", "sequence"))
ovlSeq <- unique(ovlSeq)
ovlSeq <- sampleOverlap(ovlSeq, sampleID='experiment', referenceID='sequence')
ovlSeqM <- ovlSeq$M
ovlSeqM[ovlSeqM == 1] <- NA
heatmap.2(ovlSeqM*100, 
          col=hmcol, Rowv=FALSE, Colv="Rowv", 
          cexRow=5, cexCol=5, 
          margins=c(40,40), 
          cellnote=round(ovlSeqM*100,1), 
          notecex=8,
          main="Pairwise peptide identification overlap (only peptides with at least 1 peptide detected and identified)", 
          notecol="black", 
          trace="none", 
          key = TRUE,
          keysize=1, 
          scale="none",
          symm=TRUE, 
          colsep=c(1:length(colnames(ovlSeqM))), 
          rowsep=c(1:length(rownames(ovlSeqM))), 
          sepcolor="white", 
          sepwidth=c(0.01,0.01),
          dendrogram="none")

if("fraction" %in% colnames(evidence1)){
  
  ovlSeqFx <- subset(evidence1, !is.na(intensity) & !is.na(ms.ms.count), select=c("experiment", "fraction", "sequence"))
  ovlSeqFx <- unique(ovlSeqFx)
  ovlSeqFx$experiment <- paste(ovlSeqFx$fraction, ovlSeqFx$experiment, sep=".")
  ovlSeqFx <- sampleOverlap(ovlSeqFx, sampleID='experiment', referenceID='sequence')
  ovlSeqFxM <- ovlSeqFx$M
  ovlSeqFxM[ovlSeqFxM == 1] <- NA
  heatmap.2(ovlSeqFxM*100, 
            col=hmcol, Rowv=FALSE, Colv="Rowv", 
            cexRow=0.5, cexCol=0.5, 
            margins=c(40,40), 
            #cellnote=round(ovlSeqFxM*100,1), 
            #notecex=3,
            main="Pairwise peptide identification overlap by Fraction (only peptides with at least 1 peptide detected and identified)", 
            #notecol="black", 
            trace="none", 
            key = TRUE,
            keysize=1, 
            scale="none",
            symm=TRUE, 
            colsep=c(1:length(colnames(ovlSeqFxM))), 
            rowsep=c(1:length(rownames(ovlSeqFxM))), 
            sepcolor="white", 
            sepwidth=c(0.0001,0.0001),
            dendrogram="none")
  
}


ovlProt <- subset(evidence1, !is.na(intensity) & !is.na(ms.ms.count), select=c("experiment", "proteins"))
ovlProt <- unique(ovlProt)
ovlProt <- sampleOverlap(ovlProt, sampleID='experiment', referenceID='proteins')
ovlProtM <- ovlProt$M
ovlProtM[ovlProtM == 1] <- NA
heatmap.2(ovlProtM*100, 
          col=hmcol, Rowv=FALSE, Colv="Rowv", 
          cexRow=5, cexCol=5, 
          margins=c(40,40), 
          cellnote=round(ovlProtM*100,1), 
          notecex=8,
          main="Pairwise protein identification overlap (only proteins with at least 1 peptide detected and identified)", 
          notecol="black", 
          trace="none", 
          key = TRUE,
          keysize=1, 
          scale="none",
          symm=TRUE, 
          colsep=c(1:length(colnames(ovlProtM))), 
          rowsep=c(1:length(rownames(ovlProtM))), 
          sepcolor="white", 
          sepwidth=c(0.01,0.01),
          dendrogram="none")

dev.off()


pdf('IntCorrelation%03d.pdf', width=nsamples*4, height=nsamples*4, onefile = FALSE)

peptintmtx <- subset(evidence1, !is.na(intensity), select=c("experiment", "sequence", "intensity"))
peptintmtx <- dcast(peptintmtx, sequence ~ experiment, sum)
peptintmtx <- as.matrix(peptintmtx[,-1])
peptintmtx <- log2(peptintmtx)
peptintmtx[!is.finite(peptintmtx)] <- NA
pairs(peptintmtx, upper.panel=panel.cor, diag.panel=panel.hist, lower.panel=panel.smooth, pch=".")

mpept <- peptintmtx
mpept[is.na(mpept)] <- 0
pept.pca <- prcomp(t(mpept))
pept.pcax <- as.data.frame(pept.pca$x)
pept.pcax <- merge(unique(evidence1[,c("experiment", "condition")]), pept.pcax, by.x="experiment", by.y="row.names")

ggplot(pept.pcax, aes(PC1, PC2)) + #, color=condition, shape=condition
  geom_point(alpha=.8, size=30) + 
  geom_label_repel(aes(label = experiment),
                   #label.size = 10,
                   box.padding   = 0.1, 
                   point.padding = 0.1,
                   segment.color = 'grey50') +
  ggtitle("Peptide Intensity Principal Component Analysis") +
  geom_encircle(aes(group=condition, fill=condition),alpha=0.4)




protintmtx <- subset(evidence1, !is.na(intensity), select=c("experiment", "proteins", "intensity"))
protintmtx <- dcast(protintmtx, proteins ~ experiment, sum)
protintmtx <- as.matrix(protintmtx[,-1])
protintmtx <- log2(protintmtx)
protintmtx[!is.finite(protintmtx)] <- NA
pairs(protintmtx, upper.panel=panel.cor, diag.panel=panel.hist, lower.panel=panel.smooth, pch=".")


mprot <- protintmtx
mprot[is.na(mprot)] <- 0
prot.pca <- prcomp(t(mprot))
prot.pcax <- as.data.frame(prot.pca$x)
prot.pcax <- merge(unique(evidence1[,c("experiment", "condition")]), prot.pcax, by.x="experiment", by.y="row.names")

ggplot(prot.pcax, aes(PC1, PC2)) + #, color=condition, shape=condition
  geom_point(alpha=.8, size=40) + 
  geom_label_repel(aes(label = experiment),
                   #label.size = 10,
                   box.padding   = 0.1, 
                   point.padding = 0.1,
                   segment.color = 'grey50') +
  ggtitle("Protein Intensity Principal Component Analysis") +
  geom_encircle(aes(group=condition, fill=condition),alpha=0.4)


dev.off()







pdf('SamplePrep%03d.pdf', width=nsamples*4, height=20, onefile = FALSE)

evidence.misscleavages <- data.table(subset(evidence1, ms.ms.count > 0))
misscleavages.tot <- evidence.misscleavages[, .N, by=list(experiment)]
misscleavages.dt <- evidence.misscleavages[, .N, by=list(experiment, missed.cleavages)]
misscleavages.dt <- merge(misscleavages.dt, misscleavages.tot, by="experiment", all=T)
misscleavages.dt$mc <- as.numeric(format(100*(misscleavages.dt$N.x/misscleavages.dt$N.y), digits = 5))

ggplot(misscleavages.dt, aes(x = experiment, y = mc, fill = as.factor(missed.cleavages), label=mc)) + 
  geom_bar(stat="identity", alpha=0.7) +
  geom_text(aes(label=round(mc, digits=1)), vjust=0, size = 10, position = position_stack()) +
  xlab("Experiment") + ylab("Fraction (percentage)") +
  ggtitle("Missing cleavage stats") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Set1")


oxMet.df <- ddply(evidence1, c("experiment", "condition"), summarise, 
                  pct.OxM = (length(oxidation..m.[oxidation..m.>0])/length(sequence))*100) 
ggplot(oxMet.df, aes(x = experiment, y = pct.OxM, fill = condition, label=pct.OxM)) + 
  geom_bar(stat="identity", alpha=0.7) +
  geom_text(aes(label=round(pct.OxM, digits=1)), vjust=0, size = 15, position = position_stack()) +
  xlab("Experiment") + ylab("Fraction (percentage)") +
  ggtitle("Percentage of peptides where at least 1 Methionine is oxidized") +
  theme(legend.text = element_text(size=20)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=20)) + 
  theme(axis.text.y = element_text(size=20)) + 
  theme(axis.title.x = element_text(size=30)) + 
  theme(axis.title.y = element_text(size=30)) + 
  theme(plot.title = element_text(size = 40)) + 
  scale_fill_brewer(palette="Spectral")




dev.off()










             
             
             



# Retention length
if(input$groupby == "Experiment"){
  
  g <- ggplot(evidence1(), aes(interaction(replicate, condition), retention.length)) + 
    geom_boxplot(varwidth = TRUE, aes(fill = condition), alpha=0.7) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.margin = unit(c(1, 1, 1, 0.5), "cm")) + 
    scale_fill_brewer(palette="Spectral") #+ theme(axis.title.x =element_blank())
  
} else if(input$groupby == "Condition"){
  
  
  g <- ggplot(evidence1(), aes(condition, retention.length, fill = condition)) + 
    geom_boxplot(varwidth = TRUE, aes(fill = condition), alpha=0.7) +  
    theme(axis.text.x = element_text(angle = 0, hjust = 1), plot.margin = unit(c(1, 1, 1, 0.5), "cm")) + 
    scale_fill_brewer(palette="Spectral") #+ theme(axis.title.x =element_blank())

}

######################################################

if(input$groupby == "Experiment"){
  
  g <- ggplot(evidence1(), aes(x=interaction(replicate, condition), y=retention.time)) + 
    geom_boxplot(varwidth = TRUE, aes(fill = condition), alpha=0.7) + 
    theme(axis.text.x = element_text(angle = 90), plot.margin = unit(c(1, 1, 1, 0.5), "cm")) + 
    scale_fill_brewer(palette="Spectral") 
  
} else if(input$groupby == "Condition"){
  
  g <- ggplot(evidence1(), aes(x=condition, y=retention.time)) + 
    geom_boxplot(varwidth = TRUE, aes(fill = condition), alpha=0.7) + 
    theme(axis.text.x = element_text(angle = 90), plot.margin = unit(c(1, 1, 1, 0.5), "cm")) + 
    scale_fill_brewer(palette="Spectral") 
  
} 

######################################################

# RT monitor
if(input$rt.monitor == "Gene.names"){
  
  s <- which(evidence1()$gene.names %in% input$RTMonitor)
  subsetData <- evidence1()[c(s),]
  
} else if(input$rt.monitor == "Protein.IDs"){
  
  s <- which(evidence1()$proteins %in% input$RTMonitor)
  subsetData <- evidence1()[c(s),]
  
} 


g <- ggplot(subsetData, aes(x = interaction(replicate, condition), y = as.numeric(retention.time), group=interaction(modified.sequence, charge), colour=modified.sequence)) + 
  geom_line() + 
  geom_point() +
  theme(axis.text.x = element_text(angle = 90), plot.margin = unit(c(1, 1, 1, 0.5), "cm")) #+
#facet_wrap(~gene.names, ncol=2) + guides(fill=FALSE)

######################################################

# 2D RT correlation
if(input$corrRT2DParam == "Calibrated RT") {
  
  q <- evidence1()[, c("modified.sequence", "experiment", "calibrated.retention.time")]
  
  
} else if(input$corrRT2DParam == "RT") {
  
  q <- evidence1()[, c("modified.sequence", "experiment", "retention.time")]
  
  
}

q <- dcast(q, modified.sequence ~ experiment, min)
qx <- as.matrix(q[,-1])
rownames(qx) <- q[,1]
qx <- as.matrix(qx)
qx[!is.finite(qx)] <- NA

pairs(qx, upper.panel=panel.cor, diag.panel=panel.hist, lower.panel=panel.smooth, pch=".")


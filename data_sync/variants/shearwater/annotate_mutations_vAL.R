# RUN like this:
# cd /lustre/scratch116/casm/cgp/users/fa8/Katerina/SW_RESULTS
# nohup /software/R-3.3.0/bin/Rscript katerina.R >& katerina.log


#for a in PD*; do echo $a; cd $a; Rscript /lustre/scratch116/casm/cgp/users/fa8/SHEARWATER_OESOPHAGUS/esopipe_perpatient_annotatemuts_Oesophagus.R $a >& LOG.$a.err; cd ..; done;

# Inigo Martincorena - Feb 2017
# Annotating mutations from the output of esopipe_perpatient.R
#
# fa8: I've done several modifications:
#       1. After FDR correction I do not remove all non-significant p-values, I set q<0.1
#          This is done to allow for the posterior rescue of variants: 
#             when a variant is reliable in a sample we can relax the threshold in another sample
#          In this version, variants that have q>=0.01 and q <0.1 are not removed but labelled as no-fdr
#          This rescuing worked very well for some samples but not for others.
#          All the variants "rescued" are labelled as OK-rescued.
#       2. I have included a piece of code to estimate contamination based on alternative homozygotes (AA)
#          Any read matching the reference is considered a contaminant (not necessarily true but works)
#          The code also prints the list of AA sites and read counts, allowing to identify sites that are
#          not good contamination estimators (i.e. sites that show a higher "contamination" proportion for
#          all the samples)
#       3. As with FDR, variants filtered based on indels, strandness, etc are not filtered now, they
#          are just labelled accordingly. In the end, we will only be interested in OK variants or OK + OK-rescued
#          (depending on whether we want to include the rescued variants or not)
#       4. There is a initial step, before FDR and any filtering in which consecutive indels (- or INS, 
#          but not mixed) are merged together if they have consistent VAFs (consistent meaning not considered
#          different in a Fisher's exact test). Later, if one of the indels is considered reliable (it passed
#          all the filters), all the associated indels will be included in the output.
#          ** Possible improvement: instead of comparing the VAF of indels, we could check that the indels
#             happen in the same reads
#       5. There are some variables that have to be hard coded at present (not provided as arguments), like the 
#          sample_table and the path to the bam files.
        
# Input arguments:
# 1. Dataset name (e.g. PD prefix for the patient)
#
# Example:
# nohup Rscript annotate_mutations.R dataset_name [bait_regions.bed] &



####################################################################################################
## 1. Environment
####args = commandArgs(TRUE)
####dataset_name = args[1]; #"18003";
####dataset_name = gsub(" ","",dataset_name);
####if (length(args)>1) { 
####	baits_bed = args[2];
####} else { 
####	baits_bed = "/lustre/scratch116/casm/cgp/users/fa8/Katerina/v2_panel_new.bed";
####}
####outdir = "Mutation_calls"; system(sprintf("mkdir %s",outdir));
#####bams_path = "/lustre/scratch116/casm/cgp/users/fa8/BLADDER_SHEARWATER/BAMS_REALIGNED/"; #change accordingly
#####bams_path = "/lustre/scratch112/sanger/casm/cgp/im3/Projects/Oesophagus/20170814/BAM_files/"
####bams_path = "/lustre/scratch116/casm/cgp/users/fa8/Katerina/BAM_files/";
#####sample_table = read.table(sprintf("samples_%s.txt",dataset_name), header=1, sep="\t", stringsAsFactors=F) # changed by fede to:
####sample_table = read.table("/lustre/scratch116/casm/cgp/users/fa8/Katerina/1484_samples.txt", header=1, sep="\t", stringsAsFactors=F)
#####sample_table = read.table("/lustre/scratch116/casm/cgp/users/fa8/BLADDER_SHEARWATER/1385.list", header=1, sep="\t", stringsAsFactors=F)
#####The include_REF_REF variable can be modified below to change the way contamination is estimated
#####sample_table = read.table("/lustre/scratch112/sanger/casm/cgp/im3/Projects/Oesophagus/20170814/all_samples.txt", header=0, sep="\t", stringsAsFactors=F)
#####colnames(sample_table) = c("pid","sampleID")
#####sample_table = sample_table[substr(sample_table$sampleID,1,7)==dataset_name,]


# Libraries
library("GenomicRanges")
library("rtracklayer")
library("deepSNV")

baits_bed = "/lustre/scratch116/casm/cgp/users/fa8/TARGETED_PANCREAS/normal_bait_all_covered.bed";
outdir = "Mutation_calls"; system(sprintf("mkdir %s",outdir));
bams_path = "/lustre/scratch116/casm/cgp/users/fa8/TARGETED_PANCREAS/BAMS/";
#sample_table = read.table("/lustre/scratch116/casm/cgp/users/fa8/Katerina/SCORT_katerina/tumour_samples.tsv", header=1, sep="\t", stringsAsFactors=F)
#sample_table = read.table("/lustre/scratch116/casm/cgp/users/fa8/EXOMES_PANCREAS/41j/41j_calls/tumor.list.tsv", header=F, sep="\t", stringsAsFactors=F)
sample_table = read.table("/lustre/scratch116/casm/cgp/users/fa8/TARGETED_PANCREAS/tumor.list.tsv", header=F, sep="\t", stringsAsFactors=F)
if(colnames(sample_table)[1] != "sampleID") {
	colnames(sample_table) <- c("sampleID")
}
#mutations_file = "/lustre/scratch116/casm/cgp/users/fa8/Katerina/CALLING_NEW_NORMAL_PANEL.NEWBAM2R/NEW/ALL_Mismatches.tsv";




####################################################################################################
## 2. Calling mutations from the Shearwater output files

baits = read.table(baits_bed, header=0, sep="\t", stringsAsFactors=F)
baits = GRanges(baits[,1], IRanges(baits[,2],baits[,3]))
numsegments_per_job = 100
entry_start = seq(from=1, to=length(baits), by=numsegments_per_job)
entry_end = pmin(entry_start+numsegments_per_job-1, length(baits))

####################################################################################################
# a. Loading the table of putative mutations from each patient
mutations = NULL
for (h in 1:length(entry_start)) {
	if(file.exists(sprintf("./shearwater_temp_all_targeted_pancreas/mismatches_%s_%s.txt", entry_start[h], entry_end[h]))) {
		cat("Going for file=",sprintf("./shearwater_temp_all_targeted_pancreas/mismatches_%s_%s.txt", entry_start[h], entry_end[h]),"\n");
	    m = read.table(file=sprintf("./shearwater_temp_all_targeted_pancreas/mismatches_%s_%s.txt", entry_start[h], entry_end[h]), header=1, sep="\t", stringsAsFactors=F)
		#cat("Going for file=",sprintf("shearwater_temp/mismatches_%s_%s.txt", entry_start[h], entry_end[h]),"\n");
	    #m = read.table  (file=sprintf("shearwater_temp/mismatches_%s_%s.txt", entry_start[h], entry_end[h]), header=1, sep="\t", stringsAsFactors=F)
    	#m = read.table  (file=sprintf("shearwater_%s/shearwater_temp_shearwater_%s/mismatches_%s_%s.txt", dataset_name, dataset_name, entry_start[h], entry_end[h]), header=1, sep="\t", stringsAsFactors=F)
	    mutations = rbind(mutations,m)
	}
}
#mutations <- read.table(mutations_file,header=1, sep="\t", stringsAsFactors=F)
indels_f <- length(grep("[-I]",mutations[,"mut"]));
subs_f   <- nrow(mutations)-indels_f;
cat("#INITIAL_NUMBER_OF_MUTATIONS\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");
mutations$mut_site <- paste(mutations$chr,mutations$pos,mutations$ref,mutations$mut);

# Make the matrix shorter? (if needed)
# mutations <- mutations[mutations$pval < 1e-04,]

L = sum(end(reduce(baits))-start(reduce(baits))+1) # Bait footprint
#s = sample_table$sampleID[substr(sample_table$sampleID,1,7)==dataset_name] # List of samples from this patient
#CHANGED BY FEDE:
s = unique(sample_table$sampleID) # List of samples from this patient

mutations$nu <- (mutations$xfw + mutations$xbw) / (mutations$xfw + mutations$xbw + mutations$nfw + mutations$nbw)
mutations$indiv <- substr(mutations$sampleID,1,7)

# Get samples information:
#sinfo <- read.table("/lustre/scratch116/casm/cgp/users/fa8/EXOMES_PANCREAS/Sample_info_vF.4.tsv",header=T,sep="\t",row.names=1)
#mutations$sampleID2 <- gsub(".sample.dupmarked","",mutations$sampleID)
#mutations$sample_type <- sinfo[mutations$sampleID2,"feature"]


# Let's remove individuals with less than 10 samples:
kk <- unique(mutations[,c("sampleID","indiv")])
kk <- table(kk$indiv)
good_samples <- names(kk[which(kk>=10)])
mutations <- mutations[which(mutations$indiv %in% good_samples),]


###############################################################################################################
############ fa8: Measuring potential contamination by looking into ref-reads for alternative homozygous variants
###########include_REF_REF = FALSE; # whether A>A, C>C mutations are to be included. I decided not to include them
###########                         # for some samples because they were clustering in certain regions of the genome.
###########						 # They didnt look as real mutations but mapping artefacts
###########all <- mutations[which(mutations$pval<1e-05),]
###########all <- unique(all[,c("chr","pos","ref","mut","tum_globalvaf")]);
###########all <- all[grep("[-I]",all[,"mut"],invert=T),]
############> all[which(all$tum_globalvaf>0.8),]  << very few!!!!!
############      chr       pos ref mut tum_globalvaf
############9860   12  56493822   A   C     0.9978358
############14222  16  68857441   T   C     0.9981881
############19925   2 166892788   C   C     0.9977858
############20453   2 166903445   T   T     0.9982624
############20538   2 166897864   A   A     0.9980687
############23286   2 228883721   T   C     0.9982315
############25256  20  41306600   A   G     0.9978034
############What’s there in chr2:166892788-166897864?
############SCN1a, I don’t see TRF or seg dups...
###########
###########all$mut_site <- paste(all$chr,all$pos,all$ref,all$mut);
###########mutations$vaf = (mutations$xfw+mutations$xbw)/(mutations$nfw+mutations$nbw)
############mutations[which(mutations$mut_site %in% all[which(all$tum_globalvaf>0.8),"mut_site"]),]
###########if(include_REF_REF) {
###########	alt_homoz <- all[which(all$tum_globalvaf>0.8),];
###########} else {
###########	alt_homoz <- all[which(all$tum_globalvaf>0.8 & all$mut != all$ref),];
###########}
###########if(nrow(alt_homoz) == 0) {
###########	cat("No alternative homozygotes to measure contamination (this typically happens under match-normal settings)\n");
###########} else {
###########	samples <- unique(mutations$sampleID);
###########	matches_by_sample           <- vector(length=length(samples));
###########	mismatches_by_sample        <- vector(length=length(samples));
###########	names(matches_by_sample   ) <- samples;
###########	names(mismatches_by_sample) <- samples;
###########	letters <- c("A","C","G","T","a","c","g","t");
###########	for(alt_hom in c(1:nrow(alt_homoz))) {
###########	    for(s in c(1:length(samples))) {
###########			s <- samples[s];
###########			#cat("Going to read: ", s, "\n");
###########			f = sprintf("%s/%s.bam",bams_path,s)
###########			f <- gsub(" ","",f);
###########	   		jorl <- bam2R(f, alt_homoz[alt_hom,"chr"], alt_homoz[alt_hom,"pos"], alt_homoz[alt_hom,"pos"], q=30, mask=3844, mq=10)[1,]
###########	   		match <- c(toupper(alt_homoz[alt_hom,"mut"]),tolower(alt_homoz[alt_hom,"mut"]));
###########	   		if(include_REF_REF) {
###########		   		mism  <- letters[!(letters %in% match)];
###########		   	} else {
###########		   		mism  <- c(toupper(alt_homoz[alt_hom,"ref"]),tolower(alt_homoz[alt_hom,"ref"]));
###########		   	}
###########	   		cat(alt_homoz[alt_hom,"chr"],":",alt_homoz[alt_hom,"pos"],":",alt_homoz[alt_hom,"ref"],"-",alt_homoz[alt_hom,"mut"],"\t",s,"\t",sum(jorl[match]),"\t",sum(jorl[mism]),"\n",sep="");
###########	   		matches_by_sample[s]    <- matches_by_sample[s]+sum(jorl[match]);
###########	   		mismatches_by_sample[s] <- mismatches_by_sample[s]+sum(jorl[mism ]);
###########	   	}
###########	   	cat("\n");
###########	} 
###########	mismatches_by_sample/(matches_by_sample+mismatches_by_sample);
###########}

####################################################################################################
# fa8: Let's merge neighbor indels and check they have consistent VAFs
rownames(mutations) <- paste(mutations$sampleID,mutations$chr,mutations$pos,mutations$mut,sep="|")
indels         <- mutations[which(mutations$mut=="-"),];
indels         <- indels[order(indels$sampleID, indels$chr, indels$pos),];
#rownames(indels) <- paste(indels$sampleID,indels$chr,indels$pos,indels$mut,sep="|")
i              <- 1;
group_counter  <- 1;
# fa8: rarely happening bugs in the following while loop fixed by al28
while(i <= nrow(indels)) {
  indel_from  <- indels[i,"pos"];
  indel_to    <- indel_from;
  from_index  <- i;
  to_index    <- i;
  deleted_seq <- indels[i,"ref"];
  sample_f    <- indels[i,"sampleID"];
  i <- i+1;
  for(j in c(i:nrow(indels))) {
    if(j>nrow(indels)) {
      break; #finish
    }
    else if(indels[j,"sampleID"] != indels[(j-1),"sampleID"]){
      i <- j;
      break;
    } else {
      if(indels[j,"pos"] != (indels[(j-1),"pos"]+1)) {
        i <- j; #start a new indel
        break;
      } else {
        if(indels[j,"chr"] != indels[(j-1),"chr"]) {
          i <- j; #start a new indel
          break;
        } else {
          #This is a candidate, but check their VAFs are compatible with a Fishers exact test:
          mat <- matrix(nrow=2,ncol=2,0);
          #mat[1,] <- c(indels[max(i,j-1),  "xfw"]+indels[max(i,j-1),  "xbw"], indels[max(i,j-1),  "nfw"]+indels[max(i,j-1),  "nbw"])
          mat[1,] <- c(indels[j-1,  "xfw"]+indels[j-1,  "xbw"], indels[j-1,  "nfw"]+indels[j-1,  "nbw"])
          mat[2,] <- c(indels[j,    "xfw"]+indels[j,    "xbw"], indels[j,    "nfw"]+indels[j,    "nbw"])
          mat[1,2] <- mat[1,2]-mat[1,1];
          mat[2,2] <- mat[2,2]-mat[2,1];
          pvalue <- fisher.test(mat)$p.value;
          if(pvalue < 0.01) {
            cat(" Breaking up indel because VAFs do not match\n");
            cat("             ",j-1, " vs ", j, ": pval=", pvalue, " [",indels[j-1,"pos"],"-",indels[j,"pos"],"]",sep="");
            cat("    (mat=", mat[1,1],",",mat[1,2],",",mat[2,1],",",mat[2,2],")\n",sep="");
            i <- j;
            break;
          }
          indel_to <- indels[j,"pos"];
          to_index <- j;
          deleted_seq <- paste(deleted_seq,indels[j,"ref"],sep="");
        }
      }
    }
  }
  #cat("   [Sample=",sample_f,"] Indel goes from=",indel_from,", to=", indel_to," [",deleted_seq,">-]\n",sep="");
  indels[c(from_index:to_index),"groupID"    ] <- group_counter;
  indels[c(from_index:to_index),"deleted_seq"] <- deleted_seq;
  group_counter <- group_counter + 1;
}


## while(i <= nrow(indels)) {
## 	indel_from  <- indels[i,"pos"];
## 	indel_to    <- indel_from;
## 	from_index  <- i;
## 	to_index    <- i;
## 	deleted_seq <- indels[i,"ref"];
## 	sample_f    <- indels[i,"sampleID"];
## 	if(i>nrow(indels)) {
## 		break;
## 	}	
## 	i <- i+1;
## 	for(j in c(i:nrow(indels))) {
## 		if(j>nrow(indels)) {
## 			i <- j 
## 			break; #finish
## 		}
## 		else if(indels[j,"pos"]>(indels[(j-1),"pos"]+1)) {
## 			i <- j; #start a new indel
## 			break;
## 		} else {
## 			if(indels[j,"chr"] != indels[(j-1),"chr"]) {
## 				i <- j; #start a new indel
## 				break;
## 			} else {
## 				#This is a candidate, but check their VAFs are compatible with a Fishers exact test:
## 				mat <- matrix(nrow=2,ncol=2,0);
## 				#mat[1,] <- c(indels[max(i,j-1),  "xfw"]+indels[max(i,j-1),  "xbw"], indels[max(i,j-1),  "nfw"]+indels[max(i,j-1),  "nbw"])
## 				mat[1,] <- c(indels[j-1,  "xfw"]+indels[j-1,  "xbw"], indels[j-1,  "nfw"]+indels[j-1,  "nbw"])
## 				mat[2,] <- c(indels[j,    "xfw"]+indels[j,    "xbw"], indels[j,    "nfw"]+indels[j,    "nbw"])
## 				mat[1,2] <- mat[1,2]-mat[1,1];
## 				mat[2,2] <- mat[2,2]-mat[2,1];
## 				pvalue <- fisher.test(mat)$p.value;
## 				if(pvalue < 0.01) {
## 					cat(" Breaking up indel because VAFs do not match\n");
## 					cat("             ",j-1, " vs ", j, ": pval=", pvalue, " [",indels[j-1,"pos"],"-",indels[j,"pos"],"]",sep="");
## 					cat("    (mat=", mat[1,1],",",mat[1,2],",",mat[2,1],",",mat[2,2],")\n",sep="");
## 					cat(i, " out of ", nrow(indels), "\n")
## 					i <- j;
## 					break;
## 				}		
## 				indel_to <- indels[j,"pos"];
## 				to_index <- j;
## 				i<-j
## 				deleted_seq <- paste(deleted_seq,indels[j,"ref"],sep="");
## 			}
## 		}
## 	}
## 	#cat("   [Sample=",sample_f,"] Indel goes from=",indel_from,", to=", indel_to," [",deleted_seq,">-]\n",sep="");
## 	indels[c(from_index:to_index),"groupID"    ] <- group_counter;
## 	indels[c(from_index:to_index),"deleted_seq"] <- deleted_seq;
## 	group_counter <- group_counter + 1;
## }
mutations$indel_group <- NA;
mutations$deleted_seq <- NA;
save.image("SESSION.1.Rdat")
#for(j in c(1:nrow(indels))) {
#	mutations[which(mutations$sampleID==indels[j,"sampleID"] & mutations$chr==indels[j,"chr"] 
#	                                                         & mutations$pos==indels[j,"pos"] 
#	                                                         & mutations$mut==indels[j,"mut"]),c("indel_group","deleted_seq")] = indels[j,c("groupID","deleted_seq")];
#}
mutations[rownames(indels),"indel_group"] <- indels$groupID
mutations[rownames(indels),"deleted_seq"] <- indels$deleted_seq
save.image("SESSION.2.Rdat")

#I may prefer make groups of positions... Then label each group as ok or not depending on whether they have FDR-good ones
#To get all of them, require that sample, chr, pos, mut, ref match!
#Learn how to work with lists in R!!!!!
#	new_index <- length(grouped_indels[[sample_f]])+1;
#	grouped_indels[[sample_f]][[new_index]] <- list();
#	grouped_indels[[sample_f]][[new_index]][["sampleID"]] <- sample_f;
#	grouped_indels[[sample_f]][[new_index]][["from_to"]] <- paste(indel_from,"-",indel_to,sep="");
#	grouped_indels[[sample_f]][[new_index]][["mut"]] <- deleted_seq;
#	grouped_indels[[sample_f]][[new_index]][["positions"]] <- c(indel_from:indel_to);

####################################################################################################
# b. Identifying putative germline or somatic indels to flag mutations near them
putative_indelsites = mutations[mutations$mut %in% c("-","INS"),]
#s = sample_table$sampleID[grep(dataset_name,sample_table$sampleID)] # List of samples from this patient
s = unique(sample_table$sampleID) # List of samples from this patient

putative_indelsites$qval = p.adjust(putative_indelsites$pval, method="BH", n=L*length(s)*2)
putative_indelsites = unique(putative_indelsites[putative_indelsites$qval<0.20, c("sampleID","chr","pos")])
indel_flank = 10
putative_indelsites_gr = GRanges(putative_indelsites$chr, IRanges(putative_indelsites$pos-indel_flank, putative_indelsites$pos+indel_flank))

####################################################################################################
# c. Removing germline SNPs:
#    - Any mutation present in >20% of all reads across samples (a low cutoff to remove germline indels too, as they present lower VAFs) 
# fa8: This filter is not appropriate for all cases. For example, when there is just one bladder disk in a patient.
# fa8: The removal of germline SNPs is the piece that needs more adjustments from project to project
mutations <- unique(mutations);
#mutations$indiv <- substr(mutations$sampleID,1,7)
mutations$label <- "";

#mutations <- mutations[which(mutations$tum_globalvaf<0.2),]

cat("GERMLINE: I will filter later using dbSNP\n")
# Germline filter: remove mutations present in >10% of an individual samples
kk <- unique(mutations[,c("sampleID","indiv")])
kk <- table(kk$indiv)
good_samples <- names(kk[which(kk>=10)])
mutations <- mutations[which(mutations$indiv %in% good_samples),]
for(indiv in good_samples) {
	num_samples <- kk[indiv]
	maximum <- num_samples/10
	cat(num_samples, " samples for ", indiv, "; Removing when #muts>", maximum, "\n")
	mm <- mutations[which(mutations$indiv == indiv),]
	mut_counts <- table(mm$mut_site)
	remove_muts <- names(mut_counts[which(mut_counts>maximum)])
	cat("  Removing ",length(remove_muts), " mutations for individual ", indiv, "\n")
	mutations <- mutations[-which(mutations$indiv == indiv & mutations$mut_site %in% remove_muts),]
	cat("  Remaining mutations: ", nrow(mutations), "\n")
}
save.image("SESSION.3a.Rdat")


##Get the list of mutations with nu >= 20%
#mm <- mutations[which( mutations$nu >= 0.20 ),] # if we choose  nu >= 0.10 much more are removed
#reps <- table(mm$mut_site)
#germline <- names(reps[which(reps>1)])
#mutations <- mutations[-which(mutations$mut_site %in% germline),]
#
##mutations[which(mutations$tum_globalvaf>0.2),"label"] = "germline";
###remove the germline:
##if(length(which(mutations$label == "germline")) > 0) {
##	mutations <- mutations[-which(mutations$label == "germline"),]
##}
#indels_f <- length(grep("[-I]",mutations[which(mutations$label == ""),"mut"]));
#subs_f   <- nrow(mutations[which(mutations$label == ""),])-indels_f;
#cat("#AFTER_GLOBAL_VAF\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");

#######################################################################################################
########### e. Removing mutations in muscle ("b")
########### fa8: This bit is very specific of the bladder project, for which there are muscle samples
##########if(dataset_name == "PD31007") {
##########	muts_in_muscle <- mutations[grep("PD31007r", mutations$sampleID),];
##########} else {
##########	muts_in_muscle <- mutations[grep("b.bam.out", mutations$sampleID),];
##########}
##########mutations[which(mutations$mut_site %in% muts_in_muscle$mut_site),"label"] <- paste(mutations[which(mutations$mut_site %in% muts_in_muscle$mut_site),"label"],"in_muscle;",sep="");
##########if(dataset_name == "PD31007") {
##########	mutations <- mutations[grep("PD31007r",mutations$sampleID,invert=T),]
##########} else {
##########	mutations <- mutations[grep("b.bam.out",mutations$sampleID,invert=T),]
##########}
###########For the FDR to work properly I should remove these variants before:
##########mutations <- mutations[-which(mutations$label == "in_muscle;"),]
##########
##########indels_f <- length(grep("[-I]",mutations[which(mutations$label == ""),"mut"]));
##########subs_f   <- nrow(mutations[which(mutations$label == ""),])-indels_f;
##########cat("#AFTER_MUSCLE\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");

####################################################################################################
# d. FDR calculation: significant mutations (after removing SNPs, to avoid inflating the FDR adjustment) 
L = sum(end(reduce(baits))-start(reduce(baits))+1) # Bait footprint
#s = sample_table$sampleID[grep(dataset_name,sample_table$sampleID)] # List of samples from this patient
s = unique(sample_table$sampleID) # List of samples from this patient

mutations$qval = p.adjust(mutations$pval, method="BH", n=L*length(s)*5)
mutations[which(mutations$qval>=0.01),"label"] = "no-fdr;";
prefdr.mutations <- mutations;                        # This will be the matrix used for the rescuing
mutations <- mutations[which(mutations$qval < 0.1),]; # To make the matrix smaller 
mutations = mutations[order(mutations$chr,mutations$pos),]
mutations$vaf = (mutations$xfw+mutations$xbw)/(mutations$nfw+mutations$nbw)
indels_f <- length(grep("[-I]",mutations[which(mutations$label == ""),"mut"]));
subs_f   <- nrow(mutations[which(mutations$label == ""),])-indels_f;
cat("#AFTER_FDR\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");

#rm(mm)
#rm(reps)
save.image("SESSION.3b.Rdat")


####################################################################################################
# Coverage filter: remove muts with too low cov... analyse cov stats
mutations$coverage <- apply(mutations[,c("nfw","nbw")],1,sum)
#mean_cov <- mean(mutations$coverage)
#sdev     <- sd(mutations$coverage)
#max_cov <- mean_cov + 2*sdev
#min_cov <- mean_cov - 2*sdev
#mutations <- mutations[which(mutations$coverage > 30 & mutations$coverage < 200),]
mutations <- mutations[which(mutations$coverage > 30),]
#mutations             <- mutations[which(mutations$base_counts > 30),]
indels_f <- length(grep("[-I]",mutations[which(mutations$label==""),"mut"]));
subs_f   <- nrow(mutations[which(mutations$label==""),])-indels_f;
cat("#AFTER_COVERAGE_FILTER\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");


###############################################################################################################
########## e. Removing mutations in >25% of samples
########## fa8: This bit is very specific of the oesophagus project
#########mutstr = paste(mutations$chr,mutations$pos,mutations$ref,mutations$mut)
#########mutfreq = table(mutstr)/length(s)
#########filt1 = mutations[which(mutstr %in% names(mutfreq[mutfreq>=0.25])),]; 
##########CHANGED BY FEDE: added IF
#########if(nrow(filt1) > 0) {
#########	filt1$filter = "Over_25prc_samples";
#########}
#########mutations[which(mutstr %in% names(mutfreq[which(mutfreq>=0.25)])),"label"] = paste(mutations[which(mutstr %in% names(mutfreq[which(mutfreq>=0.25)])),"label"],"over_25%;",sep="");
#########indels_f <- length(grep("[-I]",mutations[which(mutations$label == ""),"mut"]));
#########subs_f   <- nrow(mutations[which(mutations$label == ""),])-indels_f;
#########cat("#AFTER_OVER_25SAMPLES\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");

######################################################################################################
# e. FILTERING OF SNPS, there are many out there (e.g. chr1:9795726-9795726)
# Primero plotear las VAFs para ver cómo van...	
# P.e. filter variants present in 3 or more samples with at least a frequency of ....

#Algunas muestras tienen réplicas o whatever, e.g. PD29267c/d/e,  PD29596,  PD29560



####################################################################################################
# f. Requesting at least 1 supporting read from both strands and annotating substitutions near indels (somatic or germline)
#rmpos = (mutations$mut %in% c("-","INS")) & (mutations$xfw==0 | mutations$xbw==0) # Asking for 1 supporting read in both strands only for indels
rmpos = (mutations$xfw==0 | mutations$xbw==0) # Asking for 1 supporting read in both strands for all mutations
filt2 = mutations[rmpos,]; 
if(nrow(filt2)>0) {
	filt2$filter = "Strandness"
}
mutations[rmpos,"label"] = paste(mutations[rmpos,"label"],"strandness;",sep="");
#####mutations = mutations[!rmpos,]
samples = unique(mutations$sampleID)
rmpos = NULL
for (h in 1:length(samples)) {
    m = mutations[mutations$sampleID==samples[h] & !(mutations$mut %in% c("-","INS")),]
    m_gr = GRanges(m$chr, IRanges(m$pos,m$pos))
    i_gr = putative_indelsites_gr[putative_indelsites$sampleID==samples[h]]
    ol = as.matrix(suppressWarnings(findOverlaps(m_gr, i_gr, type="any", select="all")))
    rmpos = c(rmpos, rownames(m)[unique(ol[,1])])
}
filt3 = mutations[rmpos,]; 
if(nrow(filt3) > 0) {
	filt3$filter = "Near_indel"
}
mutations[rmpos,"label"] = paste(mutations[rmpos,"label"],"near_indel;",sep="");
###########mutations = mutations[!(rownames(mutations) %in% rmpos),]
write.table(mutations, file=sprintf("%s/mutations_q01.txt", outdir), col.names=T, row.names=F, sep="\t", quote=F)
#write.table(rbind(filt1,filt2,filt3), file=sprintf("%s/filteredout_q01_%s.txt", outdir, dataset_name), col.names=T, row.names=F, sep="\t", quote=F)
indels_f <- length(grep("[-I]",mutations[which(mutations$label == ""),"mut"]));
subs_f   <- nrow(mutations[which(mutations$label == ""),])-indels_f;
cat("#AFTER_BOTH_STRANDS_FILT\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");
save.image("SESSION.4.Rdat")

#######################################################################################################
#### fa8: f2. Now try to rescue mutations...
#### For each mutation seen in mutations and labelled as OK, go to the rescue of others:
#### Now this will be done re-running FDR correction but only for those sites, muts and samples of interest
####  hence reducing a lot the number of hypotheses to test (--> increasing statistical power)
###good_muts             <- unique(mutations[which(mutations$label == ""),"mut_site"]);
###mutations_tmp         <- prefdr.mutations;
###new_muts              <- mutations[0,];
###to_rescue             <- mutations_tmp[which(mutations_tmp$mut_site %in% good_muts),]
###to_rescue             <- to_rescue[-which(to_rescue$label == ""),]
###cat("Attempting to rescue ", nrow(to_rescue), " potential mutations\n");
####s = sample_table$sampleID[grep(dataset_name,sample_table$sampleID)] # List of samples from this patient
####num_tests             <- length(good_muts)*length(s) - length(which(mutations$label == "")) 
####to_rescue$qval.new    <- p.adjust(to_rescue$pval, method="BH", n=num_tests)
####to_rescue$qval.new2   <- p.adjust(to_rescue$pval, method="BH", n=L*length(s)*5)
###
####*************{{{ TODO }}}*****************
####THE MUTATIONS COME FROM THE "prefdr" DATAFRAME.
####IF THEY ARE IN THE MUTATIONS MATRIX ALREADY, I SHOULD REANNOTATE THEM IF RESCUED
####IF THEY ARE NOT IN THE MUTATIONS MATRIX, I SHOULD ADD THEM AND NOTE THAT SOME OTHER FILTERS
####   HAVE NOT BEEN RUN ON THEM
####*************{{{ TODO }}}*****************
####
###for(good in c(1:length(good_muts))) {
###	site <- good_muts[good];
###	rescue_these <- which(mutations$mut_site==site);
###	for(r in c(1:length(rescue_these))) {
###		if(mutations[rescue_these[r],"label"] != "") {
###			#if(length(grep("rescued",mutations[rescue_these[r],"label"] )) > 0) {
###			#	next;
###			#}
###			mutations[rescue_these[r],"label"] = paste(mutations[rescue_these[r],"label"],"OK-rescued;",sep="");
###		}
###	}
###}

####################################################################################################
# fa8: Write table before merging consecutive subs/indels 
#      Interesting to save this because after the merging all error labels will be lost (only OK and OK-rescued 
#      will pass)
mutations[which(mutations$label == ""),"label"] <- "OK;";
write.table(mutations, file=sprintf("%s/mutations_including_failed_ones.txt", outdir), col.names=T, row.names=F, sep="\t", quote=F)

# Because rescueing was disabled:
mutations <- mutations[grep("OK",mutations$label),         ]

####################################################################################################
# g. Annotating possible dinucleotides or runs of changes [WARNING: This annotation assumes that consecutive changes in a sample belong to a complex event. This may not be the case in 100% of the cases]
# fa8: We could compare VAFs and make sure they are consistent (Fisher's exact test) [not done yet, only for deletions at the beginning]
#    : Now I process subs, del, ins separately (they were being mixed sometimes
#    : For indels I use the group information defined at the beginning of the script

#mutations[which(mutations$label == ""),"label"] <- "OK";

ok_muts <- mutations[grep("OK",mutations$label),         ]
sub     <- ok_muts[grep("[-I]",ok_muts[,"mut"],invert=T),]
ins     <- ok_muts[grep("I",   ok_muts[,"mut"]),         ]
#del    <- ok_muts[grep("-",   ok_muts[,"mut"]),         ]  # CORRECTION SUGGESTED BY ANDREW
del     <- mutations[grep("-",   mutations[,"mut"]),     ]
sub     <- sub[order(sub$sampleID, sub$chr, sub$pos),    ]
ins     <- ins[order(ins$sampleID, ins$chr, ins$pos),    ]
del     <- del[order(del$sampleID, del$chr, del$pos),    ]

#To store the new data
mutations     <- mutations[grep("OK;",mutations$label),] # WITH THIS WE DISABLE THE RESCUEING...
mutations$vaf <- 1
mutations$qval<- 1
new_mutations <- mutations[0,]

# Deletions (defined in mutations$indel_group):
# For every "OK" deletion, get its del-groupID and find all the other deletions belonging
# to that group. Merge them and create a new entry in mutations: combine pvalues, vaf, etc
# For every "OK" deletion first check it hasn't been already merged
for(j in 1:nrow(del)) {
	indel_group       <- del[j,"indel_group"]
	if(nrow(new_mutations[which(new_mutations$indel_group==indel_group),]) > 0) {
		next; #we already have one from the group of indels
	}
	indels_from_group                                  <- mutations[which(mutations$indel_group==indel_group),]
	new_mutations                                      <- rbind(new_mutations,indels_from_group[1,])
	new_mutations[nrow(new_mutations),"pos"          ] <- min  (indels_from_group$pos              )
	new_mutations[nrow(new_mutations),"vaf"          ] <- mean (indels_from_group$vaf              )
	new_mutations[nrow(new_mutations),"tum_globalvaf"] <- mean (indels_from_group$tum_globalvaf    )
	new_mutations[nrow(new_mutations),"pval"         ] <- min  (indels_from_group$pval             )
	new_mutations[nrow(new_mutations),"qval"         ] <- min  (indels_from_group$qval             )
	new_mutations[nrow(new_mutations),"label"        ] <- paste(indels_from_group$label,collapse="")
	new_mutations[nrow(new_mutations),"ref"          ] <- indels_from_group[1,"deleted_seq"]
}

# Insertions. No need to look for consecutive INS. Just add them to new_mutations
new_mutations <- rbind(new_mutations,ins);


# Substitutions: merge consecutive... [using Iñigo's code]
d = sub$pos-(1:nrow(sub))
runs = rle(d)
rmpos = rep(0,nrow(sub))
runstarts = cumsum(runs$length)-runs$length+1
for (h in 1:length(runs$length)) {
    if (runs$length[h]>1) { # Adjacent mutations
        mutcluster                         = runstarts[h]:(runstarts[h]+runs$lengths[h]-1)
        rmpos[mutcluster[-1]             ] = 1 # Removing all the affected rows except the first one (which we will edit to capture the complex event)
        sub[mutcluster[1],"ref"          ] = paste(sub[mutcluster,"ref"          ],collapse="")
        sub[mutcluster[1],"mut"          ] = paste(sub[mutcluster,"mut"          ],collapse="")
        sub[mutcluster[1],"mu"           ] = mean (sub[mutcluster,"mu"           ]            )
        sub[mutcluster[1],"tum_globalvaf"] = mean (sub[mutcluster,"tum_globalvaf"]            )
        sub[mutcluster[1],"vaf"          ] = mean (sub[mutcluster,"vaf"          ]            )
        sub[mutcluster[1],"pval"         ] = min  (sub[mutcluster,"pval"         ]            )
        sub[mutcluster[1],"qval"         ] = min  (sub[mutcluster,"qval"         ]            )
    }
}
sub = sub[!rmpos,]        
new_mutations <- rbind(new_mutations,sub);

mutations.old <- mutations
mutations     <- new_mutations
mutations[which(mutations$label == ""),"label"] <- "OK;";
mutations <- mutations[grep("OK",mutations$label),      ]
write.table(mutations, file=sprintf("%s/mutations_q01.annotruns.1.txt", outdir), col.names=T, row.names=F, sep="\t", quote=F)
indels_f <- length(grep("[-I]",mutations[,"mut"]));
subs_f   <- nrow(mutations)-indels_f;
cat("#AFTER_MERGING_RUNS\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");

####################################################################################################
# h. Filtering out mutations if reads carry indels within 10 bp of the mutation (allowing for mapQ<30)
#    Using mapQ>=30 in the calling of indels tends to remove reads with genuine indels and
#    so misses potential false positives associated with the indel

flank = 15 # bp around the mutation where we will look for indels
rmpos = NULL
for (h in 1:length(samples)) {
	cat("sample=",h," out of=", length(samples),"\n")
    m = mutations[mutations$sampleID==samples[h] & !(mutations$mut %in% c("-","INS")),]
    s <- samples[h]
    if (nrow(m)>0) { 
        for (r in 1:nrow(m)) {
        	f <- NULL;
        	#if(length(grep("out",m$sampleID[r]))>=1) {
        		f = sprintf("/lustre/scratch116/casm/cgp/users/fa8/TARGETED_PANCREAS/BAMS/%s.bam",m$sampleID[r])
			#	f = sprintf("%s/%s.bam",bams_path,s)
				#cat("Going to read: ", f, "\n");
        	#} else {
	        #    f = sprintf("/lustre/scratch116/casm/cgp/users/fa8/BLADDER_SHEARWATER/BAMS/%s.bam",m$sampleID[r])
	        #}
			f <- gsub(" ","",f);
			#cat("     ", m$chr[r], m$pos[r]-flank, m$pos[r]+flank," f=",f,"\n",sep=" : ")
            n_indels = sum(bam2R(f, m$chr[r], m$pos[r]-flank, m$pos[r]+flank, q=30, mask=3844, mq=10)[,c("-","INS","_","ins")]) # Number of reads with an indel around the mutation
            cat("n_indels=",n_indels," for mutation=",r, " of sample=",s,", xfw+xbw=",m$xfw[r]+m$xbw[r],"\n")
            #if (n_indels > 5*(m$xfw[r]+m$xbw[r])) { 
            if (n_indels > 3*(m$xfw[r]+m$xbw[r])) {  # Modified by Fede
                rmpos = c(rmpos, rownames(m)[r])
            }
        }
    }
}

if (length(rmpos)>0) {
    filt4 = mutations[rmpos,]; filt4$filter = "Near_hidden_indel"
	mutations[which(mutations$label == "OK;"),"label"] <- "";
	mutations[rmpos,"label"] = paste(mutations[rmpos,"label"],"near_hidden_indel;",sep="");
	mutations[which(mutations$label == ""),"label"] <- "OK;";
    write.table(mutations, file=sprintf("%s/mutations_q01.txt", outdir), col.names=T, row.names=F, sep="\t", quote=F)
    #write.table(rbind(filt1,filt2,filt3,filt4), file=sprintf("%s/filteredout_q01_%s.txt", outdir, dataset_name), col.names=T, row.names=F, sep="\t", quote=F)
}
write.table(mutations, file=sprintf("%s/mutations_q01.annotruns.2.txt", outdir), col.names=T, row.names=F, sep="\t", quote=F)
indels_f <- length(grep("[-I]",mutations[grep("OK",mutations$label),"mut"]));
subs_f   <- nrow(mutations[grep("OK",mutations$label),])-indels_f;
cat("#AFTER_NEARBY_INDELS\t",nrow(mutations),"\t",indels_f,"\t",subs_f,"\t",indels_f/(indels_f+subs_f),"\n",sep="");

write.table(mutations[grep("OK;",mutations$label),], file=sprintf("%s/mutations_q01.FINAL-onlyOK.txt", outdir), col.names=T, row.names=F, sep="\t", quote=F)

# An specific filter:
mutations <- mutations[which(mutations$sampleID != "PD37586b_lo0005.sample.dupmarked"),]
write.table(mutations[grep("OK;",mutations$label),], file=sprintf("%s/mutations_q01.FINAL-onlyOK.txt", outdir), col.names=T, row.names=F, sep="\t", quote=F)


#########################################################################################################################
# In the laptop:
#########################################################################################################################
load("~/Desktop/TARGETED_PANCREAS/SESSION.3b.Rdat")
#muts <- read.table("~/Desktop/Pancreas_Islets_Exomes/mutations_q01.FINAL-onlyOK.txt",header=T,sep="\t")
#muts$sampleID3 <- paste(muts$sampleID2,muts$sample_type,sep="-")
muts <- mutations
muts$mut_id <- paste(muts$chr,muts$pos,muts$ref,muts$mut,sep=":")
tt          <- t(table(muts[,c("sampleID","mut_id")]))
counts      <- apply(tt,1,sum)
tt_shared   <- tt[which(counts>=2),] # shared by at least two samples

#heatmap(tt_shared,scale="none")

for(sample in colnames(tt_shared)) {
	for(mut_id in rownames(tt_shared)) {
		if(tt_shared[mut_id,sample]>0) {
			tt_shared[mut_id,sample] <- muts[which(muts$sampleID==sample&muts$mut_id==mut_id),"vaf"]
		}
	}
}
heatmap(tt_shared,scale="none",col=colorRampPalette(c("aliceblue","lightblue","darkblue","black"))(10))


for(indiv in good_samples) {
	muts <- mutations[which(mutations$indiv==indiv),]
	muts$mut_id <- paste(muts$chr,muts$pos,muts$ref,muts$mut,sep=":")
	tt          <- t(table(muts[,c("sampleID","mut_id")]))
	counts      <- apply(tt,1,sum)
	tt_shared   <- tt[which(counts>=2),] # shared by at least two samples
	
	if(length(which(counts>=2))>2) {
		
		#for(sample in colnames(tt_shared)) {
		#	for(mut_id in rownames(tt_shared)) {
		#		if(tt_shared[mut_id,sample]>0) {
		#			tt_shared[mut_id,sample] <- muts[which(muts$sampleID==sample&muts$mut_id==mut_id),"vaf"]
		#		}
		#	}
		#}
		pdf(paste("~/Desktop/TARGETED_PANCREAS/",indiv,".heatmap-shared-novafs2.pdf",sep=""))
		par(mar=c(10,10,10,10))
		colnames(tt_shared) <- gsub(".sample.dupmarked","",colnames(tt_shared))
		heatmap(tt_shared,scale="none",col=colorRampPalette(c("aliceblue","lightblue","darkblue","black"))(20),cexRow=0.6,cexCol=0.3)
		dev.off();
	}
	
}



remove_these   <- rownames(tt[which(counts>=4),])
muts <- muts[-which(muts$mut_id %in% remove_these),]

write.table(muts,"~/Desktop/Pancreas_Islets_Exomes/mutations_q01.FINAL-onlyOK-inlessthan4samples.txt",sep="\t",quote=F)


par(mar=c(13,3,3,3))
boxplot(muts$vaf~muts$sampleID3,las=2)

par(mar=c(13,3,3,3))
counts <- sort(table(muts$sampleID3))
samples <- names(sort(table(muts$sampleID3)))
boxplot(NULL,xlim=c(0,115),ylim=c(0,0.7))
for(i in c(1:length(samples))) {
	sample <- samples[i]
	boxplot(muts[which(muts$sampleID3==sample),"vaf"],add=T,at=i)
	text(i,0.7,counts[sample],cex=.6)
}
axis(side=1,at=c(1:length(samples)),labels=samples,las=2,cex=.7)

barplot(sort(table(muts$sampleID3)),las=2,cex.names=.7)



############################################################################################################
# Let's run dndscv
mutations <- read.table("~/Desktop/TARGETED_PANCREAS/Mutation_calls/mutations_q01.FINAL-onlyOK.txt",sep="\t",header=T)
target_genes <- read.table("~/Desktop/TARGETED_PANCREAS/target_genes.txt",header=F)[,1]
library(dndscv)
ddall<-dndscv(mutations,gene_list=target_genes)

############################################################################################################
# Let's run dndscv including mutations failing the strandness filter
load("~/Desktop/TARGETED_PANCREAS/SESSION.3c.Rdat")
ddall2<-dndscv(mutations,gene_list=target_genes)



############################################################################################################################
## # To plot the putative germline calls:
## load("~/Desktop/Pancreas_Islets_Exomes/putative_germline.Rdat")
## muts <- putative_germline
## 
## muts$sampleID3 <- paste(muts$sampleID2,muts$sample_type,sep="-")
## muts$mut_id <- paste(muts$chr,muts$pos,muts$ref,muts$mut,sep=":")
## tt          <- t(table(muts[,c("sampleID3","mut_id")]))
## counts      <- apply(tt,1,sum)
## tt_shared   <- tt[which(counts>=2),] # shared by at least two samples
## 
## #heatmap(tt_shared,scale="none")
## 
## for(sample in colnames(tt_shared)) {
## 	for(mut_id in rownames(tt_shared)) {
## 		if(tt_shared[mut_id,sample]>0) {
## 			tt_shared[mut_id,sample] <- muts[which(muts$sampleID3==sample&muts$mut_id==mut_id),"vaf"]
## 		}
## 	}
## }
## heatmap(tt_shared,scale="none",col=colorRampPalette(c("aliceblue","lightblue","darkblue","black"))(10))
## 

#########################################################################################################################

######
# Back to the farm
library("GenomicRanges")
library("rtracklayer")
library("deepSNV")
library("dndscv")
load("CHECK_THIS_FROM_HOME2.Rdat")
load("sup_conco.Rdat")
target_genes <- read.table("target_genes.txt",header=F)[,1]
sups         <- intersect(suppressors,target_genes)
oncs         <- intersect(oncogenes,target_genes)
dsup         <- dndscv(mutations,gene_list=sups,max_muts_per_gene_per_sample=12)
donc         <- dndscv(mutations,gene_list=oncs,max_muts_per_gene_per_sample=12)
dall         <- dndscv(mutations,gene_list=unique(sups,oncs),max_muts_per_gene_per_sample=12)
ddother      <- dndscv(mutations,gene_list=setdiff(target_genes,unique(sups,oncs)),max_muts_per_gene_per_sample=12)
save(ddother,dsup,donc,dall,file="dndscv_results.prefdr.Rdat")


